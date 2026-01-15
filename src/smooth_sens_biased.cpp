#include "counting.h"
#include "distribution.h"

#include "distance_retrieval.h"
#include "smooth_sens_unbiased.h"
#include "smooth_sens_utils.h"

bool shift_set_can_be_increased(SingleTargetDistance &std, int k, double beta) {
    auto [next_dis, _] = std.k_th_closest_distance(k + 1);
    return static_cast<double>(k) / (k + 1) < std::exp(-beta * next_dis);
}

int shift_count(SingleTargetDistance &std, double beta) {
    int lo = 1;
    int hi = std.size();

    while (lo < hi) {
        int mid = (lo + hi) / 2;
        bool can_be_increased = shift_set_can_be_increased(std, mid, beta);

        if (can_be_increased) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }

    return lo;
}

double obj(SingleTargetDistance &std, int k, double beta) {
    long long dis_sum = std.sum_k_closest_distances(k);
    return k * std::exp(-beta * dis_sum);
}

double fixed_edge_sensitivity(Graph &g,
                              Edge e,
                              const int lambda,
                              const std::list<int> &triangle_index_list,
                              const std::vector<Triangle> &triangles,
                              double beta,
                              bool inc) {
    std::vector<int> partial_triangle_weights = setup_partial_triangle_weights(g, e, triangle_index_list, triangles);

    if (partial_triangle_weights.empty()) {
        return 0;
    }

    auto [targets, zero_target] = setup_targets_biased(g, e, lambda, partial_triangle_weights, inc);
    SingleTargetDistance std = SingleTargetDistance(partial_triangle_weights, targets[0]);

    double opt = 0;

    for (int i = 0; i < targets.size(); ++i) {
        // for every target compute the optimal shift set
        // distance from current target to zero target (z_i)
        int k = shift_count(std, beta);
        double opt_t = obj(std, k, beta) * std::exp(-beta * std::abs(targets[i] - zero_target));

        // check if new opt was found
        opt = std::max(opt, opt_t);

        // advance target
        if (i < targets.size() - 1) {
            std.update_target(targets[i + 1]);
        }
    }

    return opt;
}

double smooth_sensitivity(Graph &g, Node v, const int lambda, const std::list<int> &triangles_index_list,
                          const std::vector<Triangle> &triangles, double beta) {
    double sens = 0;

    // #pragma omp parallel for reduction(max:sens)
    for (int i = 0; i < boost::degree(v, g); ++i) {
        // Fix the edge e = (v,u) and compute the maximum smooth sensitivity achieved by increasing / decreasing e
        Node u = boost::adjacent_vertices(v, g).first[i];
        Edge e = boost::edge(u, v, g).first;

        double sens_inc = fixed_edge_sensitivity(g, e, lambda, triangles_index_list, triangles, beta, true);
        double sens_dec = fixed_edge_sensitivity(g, e, lambda, triangles_index_list, triangles, beta, false);

        double fixed_edge_sens = std::max(sens_inc, sens_dec);
        sens = std::max(fixed_edge_sens, sens);
    }

    return sens;
}

void apply_smooth_sensitivity(Graph &g,
                              const PrivateCountingConfig &cfg,
                              std::vector<TriangleCount> &counts,
                              std::vector<std::list<int> > &node_triangle_index_map,
                              const std::vector<Triangle> &triangles) {
    const double beta = cfg.count_eps / (2 * (cfg.gamma - 1));
    const double p = std::exp(-cfg.weight_eps);
    const double smooth_sens_mult = 2 * std::pow(cfg.gamma - 1, (cfg.gamma - 1) / cfg.gamma);
    auto rv = PolynomialTailRV(cfg.gamma);

#pragma omp parallel for
    for (int i = 0; i < counts.size(); i++) {
        auto &triangle_index_list = node_triangle_index_map[i];

        if (triangle_index_list.empty()) {
            continue;
        }

        double sens = smooth_sensitivity(g, i, cfg.lambda, triangle_index_list, triangles, beta);
        double unbiased_sens = smooth_sensitivity_unb(g, i, cfg.lambda, triangle_index_list, triangles, beta, p);

        counts[i].unbiased += smooth_sens_mult / cfg.count_eps * unbiased_sens * rv.sample();
        counts[i].biased += smooth_sens_mult / cfg.count_eps * sens * rv.sample();
    }
}
