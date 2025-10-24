#include "counting.h"
#include "distribution.h"
#include <boost/log/utility/setup/console.hpp>
#include "treap.h"


std::pair<std::vector<double>, double> compute_partial_triangle_weights(Graph &g,
                                                                        const Edge &e,
                                                                        const int lamba,
									const std::list<int> &triangle_index_list,
                                                                        const std::vector<Triangle> &triangles,
									bool inc) {
    int w_e = Triangle::get_edge_weight(g, e);
    std::vector<double> partial_triangle_weights;
        
    for (int t_id : triangle_index_list) {
	auto t = triangles[t_id];
        if (!t.contains_edge(g, e))
            continue;

        auto [noise_e1, w_e1, w_e2, w_e3] = t.get_triangle_weights(g);

        double partial_triangle_weight = noise_e1 + w_e1 + w_e2 + w_e3 - w_e;
        partial_triangle_weights.push_back(partial_triangle_weight);
    }

    std::sort(partial_triangle_weights.begin(), partial_triangle_weights.end());

    double target;
    if (inc) {
        target = lamba - 1 - w_e;
    } else {
        target = lamba - w_e;
    }

    return {partial_triangle_weights, target};
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

double fixed_edge_sensitivity(Graph &g,
                              Edge e,
                              const int lambda, 
			      const std::list<int> &triangle_index_list,
			      const std::vector<Triangle> &triangles,
			      double beta,
                              bool inc) {
    double opt = 0;

    auto [partial_triangle_weights, target] = compute_partial_triangle_weights(g, e, lambda, triangle_index_list, triangles, inc);

    if (partial_triangle_weights.empty()) {
        return -1;
    }

    // setup left and right tree
    TreapNode *left = nullptr;
    TreapNode *right = nullptr;
    double offL = 0;
    double offR = 0;

    if (partial_triangle_weights[0] <= target) {
        left = insert(left, 0);
        for (int i = 1; i < partial_triangle_weights.size(); i++) {
            right = insert(right, abs(partial_triangle_weights[0] - partial_triangle_weights[i]));
        }
    } else {
        // edge case where all partial triangle weights are to the right of the target
        // -> left tree must be empty
        for (double partial_triangle_weight: partial_triangle_weights) {
            right = insert(right, abs(target - partial_triangle_weight));
        }
    }

    std::vector<double> centers = build_targets(partial_triangle_weights, target);
    for (int j = 0; j < centers.size(); j++) {
        auto [t, opt_l] = optimal_shift_set_size(left, right, offL, offR, beta);
        double n = fixed_edge_objective(left, right, offL, offR, t, opt_l,
                                        abs(centers[j] - target), beta);
        opt = std::max(opt, n);
        std::tie(left, right) = advance_center(left, right, offL, offR, centers, j);
    }

    return opt;
}

double smooth_sensitivity(Graph &g, Node v, const int lambda, const std::list<int> &triangles_index_list, const std::vector<Triangle> &triangles, double beta) {
    double sens = 0;

    #pragma omp parallel for reduction(max:sens)
    for (int i = 0; i < boost::degree(v, g); ++i) {
        // Fix the edge e = (v,u) and compute the maximum smooth sensitivity achieved by increasing / decreasing e
        Node u = boost::adjacent_vertices(v, g).first[i];
        Edge e = boost::edge(u, v, g).first;

        double sens_inc = fixed_edge_sensitivity(g, e, lambda, triangles_index_list, triangles, beta, true);
        double sens_dec = fixed_edge_sensitivity(g, e, lambda, triangles_index_list, triangles, beta, false);

        // std::cout << "Edge " << e << ", sens_inc: " << sens_inc << ", sens_dec: " << sens_dec << std::endl;

        double fixed_edge_sens = std::max(sens_inc, sens_dec);
        sens = std::max(fixed_edge_sens, sens);
    }

    return sens;
}

void apply_smooth_sensitivity(Graph &g,
                              const PrivateCountingConfig &cfg,
                              std::vector<TriangleCount> &counts,
                              std::vector<std::list<int>> &node_triangle_index_map,
			      const std::vector<Triangle> &triangles) {
    const double beta = cfg.count_eps / (2 * (cfg.gamma - 1));
    const double p = std::exp(-cfg.weight_eps);
    const double unbiased_sens_mult = 1 + 2 * (p / std::pow(1 - p, 2));
    const double smooth_sens_mult = 2*std::pow(cfg.gamma - 1, (cfg.gamma - 1) / cfg.gamma);
    auto rv = PolynomialTailRV(cfg.gamma);

    #pragma omp parallel for
    for (int i = 0; i < counts.size(); i++) {
        auto &triangle_index_list = node_triangle_index_map[i];

        if (triangle_index_list.empty()) {
            continue;
        }

        double sens = smooth_sensitivity(g, i, cfg.lambda, triangle_index_list, triangles, beta);

        counts[i].unbiased += smooth_sens_mult / cfg.count_eps * unbiased_sens_mult * sens * rv.sample();
        counts[i].biased += smooth_sens_mult / cfg.count_eps * sens * rv.sample();
    }
}
