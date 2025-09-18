#include "basic_counting.h"
#include "distribution.h"
#include "treap.h"
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>


std::pair<std::vector<double>, double> compute_partial_triangle_weights(const Graph &g, int v, int i, bool inc) {
    const Graph::vertex_descriptor root = boost::vertex(v, g);
    auto weight_map = get(&edge_info::weight, g);
    auto noise_map = get(&edge_info::noise, g);

    std::vector<int> neighbors;
    for (auto [nbrIt, nbrEnd] = boost::adjacent_vertices(root, g); nbrIt != nbrEnd; ++nbrIt) {
        neighbors.push_back(*nbrIt);
    }

    int u = neighbors[i];

    // BOOST_LOG_TRIVIAL(debug) << "v: " << v << " u: " << u;

    std::vector<double> partial_triangle_weights;
    auto [e_uv, exists_uv] = boost::edge(u, v, g);

    for (size_t j = 0; j < neighbors.size(); ++j) {
        if (i == j)
            continue;

        int w = neighbors[j];
        if (w < v)
            continue;

        auto [e_vw, exists_vw] = boost::edge(v, w, g);
        auto [e_uw, exists_uw] = boost::edge(u, w, g);

        if (exists_uw) {
            double partial_triangle_weight = get(noise_map, e_uw) + get(weight_map, e_uw) + get(weight_map, e_vw);
            // BOOST_LOG_TRIVIAL(debug) << "u: " << u << " v: " << v << " w: " << w << " partial_triangle_weight: " << partial_triangle_weight;
            partial_triangle_weights.push_back(partial_triangle_weight);
        }
    }

    std::sort(partial_triangle_weights.begin(), partial_triangle_weights.end());

    if (inc) {
        return {partial_triangle_weights, -1 - get(weight_map, e_uv)};
    } else {
        return {partial_triangle_weights, 1 - get(weight_map, e_uv)};
    }
}

std::pair<double, int> handle_boundary_case(TreapNode *left,
                                            TreapNode *right,
                                            double offL,
                                            double offR,
                                            int k,
                                            int m,
                                            double opt_dis,
                                            int opt_l) {
    TreapNode *primary = (m == 0 ? right : left);
    TreapNode *secondary = (m == 0 ? left : right);
    double offP = (m == 0 ? offR : offL);
    double offS = (m == 0 ? offL : offR);

    double a = kth(primary, k)->key + offP;
    int a_l = (m == 0 ? 0 : k);
    double b = std::max(kth(primary, k - 1)->key + offP, kth(secondary, 1)->key + offS);
    int b_l = (m == 0 ? 1 : k - 1);

    if (a < b && a < opt_dis) {
        return {a, a_l};
    }

    if (b < opt_dis && b < a) {
        return {b, b_l};
    }

    return {opt_dis, opt_l};
}

double handle_nullptr_case(TreapNode *root, double off, int k) {
    auto out = kth(root, k);

    if (out == nullptr) {
        throw std::invalid_argument("find_k_smallest_distance: no k values in left and right treap");
    }

    return out->key + off;
}

std::pair<double, int> find_k_smallest_distance(TreapNode *left, TreapNode *right, int k, double offL, double offR) {
    if (left == nullptr) {
        return {handle_nullptr_case(right, offR, k), 0};
    }
    if (right == nullptr) {
        // std::cout << "find_k_smallest_distance: right is nullptr." << std::endl;
        return {handle_nullptr_case(left, offL, k), k};
    }
    if (k == 1) {
        double l_value = kth(left, 1)->key + offL;
        double r_value = kth(right, 1)->key + offR;

        if (l_value <= r_value) {
            return {l_value, 1};
        }

        return {r_value, 0};
    }

    int left_size = get_size(left);
    int right_size = get_size(right);

    int bound_l = std::max(k - right_size, 0);
    int bound_r = std::min(left_size, k);

    double opt_dis = std::numeric_limits<double>::max();
    int opt_l = 0;

    while (bound_l <= bound_r) {
        const int m = std::min((bound_l + bound_r) / 2, left_size);

        if (m == 0 || m == k) {
            return handle_boundary_case(left, right, offL, offR, k, m, opt_dis, opt_l);
        }

        double l_value = kth(left, m)->key + offL;
        double r_value = kth(right, k - m)->key + offR;
        double n_dis = std::max(l_value, r_value);

        if (n_dis < opt_dis) {
            opt_dis = n_dis;
            opt_l = m;
        }

        if (l_value == r_value) {
            return {opt_dis, m};
        }

        if (l_value < r_value) {
            bound_l = m + 1;
        }

        if (l_value > r_value) {
            bound_r = m - 1;
        }
    }

    return {opt_dis, opt_l};
}


double compute_sum(TreapNode *left, TreapNode *right, double offL, double offR, int l, int r) {
    return sum_first_k(left, l) + offL * l + sum_first_k(right, r) + offR * r;
}

std::pair<bool, int> shift_set_can_be_increased(int t, double beta, TreapNode *left, TreapNode *right, double offL,
                                                double offR) {
    auto [next_dis, i] = find_k_smallest_distance(left, right, t + 1, offL, offR);
    // std::cout << "t: " << t << ", next_dis: " << next_dis << std::endl;
    return {static_cast<double>(t) / (t + 1) < std::exp(-beta * next_dis), i};
}

std::pair<int, int> optimal_shift_set_size(TreapNode *left, TreapNode *right, const double offL,
                                           const double offR, const double beta) {
    int bound_l = 1;
    int bound_r = get_size(left) + get_size(right);

    int opt_l = 1;
    bool can_be_increased;

    while (bound_l < bound_r) {
        int m = (bound_l + bound_r) / 2;
        int l;
        std::tie(can_be_increased, l) = shift_set_can_be_increased(m, beta, left, right, offL, offR);
        // std::cout << "optimal_shift_set_size: can_be_increased: " << can_be_increased << ", m: " << m <<  ", l: " << l << ", opt_l: " <<
        //         opt_l << std::endl;
        if (can_be_increased) {
            bound_l = m + 1;
            opt_l = l;
        } else {
            bound_r = m;
        }
    }

    return {bound_l, opt_l};
}

double fixed_edge_objective(TreapNode *left,
                            TreapNode *right,
                            double offL,
                            double offR,
                            int t,
                            int l,
                            double center_dis,
                            double beta) {
    double sum_dis = compute_sum(left, right, offL, offR, l, t - l) + center_dis;
    double obj = t * std::exp(-beta * sum_dis);
    return obj;
}

std::pair<TreapNode *, TreapNode *> advance_center(TreapNode *left,
                                                   TreapNode *right,
                                                   double &offL,
                                                   double &offR,
                                                   const std::vector<double> &centers,
                                                   int j) {
    if (j >= centers.size() - 1) {
        return {left, right};
    }

    double next_center = centers[j + 1];
    double delta = abs(centers[j] - next_center);

    offL = offL + delta;
    offR = offR - delta;

    if (contains(right, abs(centers[0] - next_center))) {
        right = erase(right, abs(centers[0] - next_center));
        left = insert(left, -offL);
    }

    return {left, right};
}

std::vector<double> build_targets(const std::vector<double> &centers, double target) {
    std::vector<double> result = centers;
    auto it = std::lower_bound(result.begin(), result.end(), target);
    result.insert(it, target);
    return result;
}

double fixed_edge_sensitivity(const Graph &g, int v, int i, double beta, bool inc) {
    // BOOST_LOG_TRIVIAL(info) << "Compute maximum objective value for fixed edge: " << i;
    double opt = 0;

    auto [partial_triangle_weights, target] = compute_partial_triangle_weights(g, v, i, inc);

    if (partial_triangle_weights.size() == 0) {
        return -1;
    }

    // setup left and right tree
    TreapNode *left = nullptr;
    TreapNode *right = nullptr;
    double offL = 0;
    double offR = 0;

    // TODO: If all partial triangle weights are larger than the target we cannot insert 0 as the first key
    if (partial_triangle_weights[0] <= target) {
        left = insert(left, 0);

        for (int j = 1; j < partial_triangle_weights.size(); j++) {
            right = insert(right, abs(partial_triangle_weights[0] - partial_triangle_weights[j]));
        }
    } else {
        // edge case where all partial triangle weights are to the right of the target
        // -> left tree must be empty
        for (int j = 0; j < partial_triangle_weights.size(); j++) {
            right = insert(right, abs(target - partial_triangle_weights[j]));
        }
    }

    std::vector<double> centers = build_targets(partial_triangle_weights, target);
    // BOOST_LOG_TRIVIAL(info) << "Centers:";

    // print_left_right_treap(left, right);
    for (int j = 0; j < centers.size(); j++) {
        auto [t, opt_l] = optimal_shift_set_size(left, right, offL, offR, beta);
        double n = fixed_edge_objective(left, right, offL, offR, t, opt_l,
                                        abs(centers[j] - target), beta);

        // BOOST_LOG_TRIVIAL(info) << "Center: " << centers[j]
        //                         << ", t: " << t
        //                         << ", obj: " << n;

        opt = std::max(opt, n);

        std::tie(left, right) = advance_center(left, right, offL, offR, centers, j);
        // print_left_right_treap(left, right);
    }

    return opt;
}

double smooth_sensitivity(const Graph &g, int v, double beta) {
    double sens = 0;

    for (int i = 0; i < boost::degree(v, g); ++i) {
        int u = boost::adjacent_vertices(v, g).first[i];
        if (v > u) continue;

        double sens_inc = fixed_edge_sensitivity(g, v, i, beta, true);
        double sens_dec = fixed_edge_sensitivity(g, v, i, beta, false);

        double fixed_edge_sens = std::max(sens_inc, sens_dec);
        sens = std::max(fixed_edge_sens, sens);
    }

    return sens;
}

void apply_smooth_sensitivity(const Graph &g, const PrivateCountingConfig &cfg, std::vector<double> &counts) {
    const double gamma = 4;
    const double beta = cfg.count_eps / (2 * (gamma - 1));
    const double sens_mult = cfg.use_unbiased_estimator
                                 ? (1 + 2 * (std::exp(cfg.weight_eps) / std::pow(1 - std::exp(cfg.weight_eps), 2)))
                                 : 1;
    auto rv = PolynomialTailRV(gamma);

    // BOOST_LOG_TRIVIAL(debug) << "Smooth sensitivity multiplier: " << sens_mult << " for use_unbiased_estimator: " << cfg.use_unbiased_estimator;

    for (int i = 0; i < counts.size(); i++) {
        double sens = smooth_sensitivity(g, i, beta);
        double sample = rv.sample();
//         BOOST_LOG_TRIVIAL(debug) << "Smooth sensitivity: " << sens << ", laplace b: " << sens * sens_mult / cfg.
// count_eps << ", noise sample: " << sample << std::endl;
        double noise = 2 * (gamma - 1) / cfg.count_eps * sens_mult * sens * sample;
        counts[i] += noise;
    }
}
