#include "../include/smooth_sens_unbiased.h"

#include "distance_retrieval.h"
#include "smooth_sens_utils.h"
#include "triangle.h"
#include "utils.h"

bool shift_set_can_be_increased(TargetDistance &td, int k, int m, int n, double beta) {
    auto [next_dis, _] = td.k_th_closest_distance(k + m + n + 1);
    return static_cast<double>(k + m + n) / (k + m + n + 1) < std::exp(-beta * next_dis);
}

/**
 * Returns objective function (without fixed factor e^(-beta * |z_i|)) during the shifting of values from O
 */
double obj_out_target(
    TargetDistance &dtd,
    int k,
    int m,
    int n,
    double f,
    double beta
) {
    long long sum_k = dtd.sum_k_closest_distances(k + m + n);
    double linear = f * (k + n + m);
    double exponent = -beta * sum_k;
    return linear * std::exp(exponent);
}

/**
 *  Computes the optimal number of values shifted from the outside set (O) to the target via binary search.
 */
int opt_outside_shift_count(TargetDistance &td, int m, int n, double beta) {
    int lo = 0;
    int hi = td.outside_target_count();

    while (lo < hi) {
        int mid = (lo + hi) / 2;
        bool can_be_increased = shift_set_can_be_increased(td, mid, m, n, beta);
        if (can_be_increased) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }

    return lo;
}

/**
 * Returns objective function (without fixed factor e^(-beta * |z_i|)) during the shifting of values from N
 */
double obj_on_target_pos(int k, int m, int n, double x, double beta) {
    // objective function during the step of assigning on target values for positive sum contribution
    double linear = ((m + k) * x) + ((n - k) * (-1.0 - 2.0 * x));
    return linear * std::exp(-beta * k);
}

/**
 * Return optimal number of values shifted from N
 */
int opt_on_target_shift_count_double(int m, int n, double x, double beta) {
    double y = m * x + (-1.0 - 2.0 * x) * n;
    double k_star = 1.0 / beta + y / (1.0 + 3.0 * x);

    // Candidate integer values (clamped to [0, n])
    int k0 = std::clamp((int) std::floor(k_star), 0, n);
    int k1 = std::clamp(k0 + 1, 0, n);

    // Evaluate objective
    double v0 = obj_on_target_pos(k0, m, n, x, beta);
    double v1 = obj_on_target_pos(k1, m, n, x, beta);

    return (v1 > v0) ? k1 : k0;
}

double opt_double_target_obj(DoubleTargetDistance dtd, double x, double beta) {
    // Compute |M|, |N|
    int m = dtd.adjacent_to_target_count();
    int n = dtd.on_target_count();

    // Compute optimal number of values connected from N
    int on_t_k_star = opt_on_target_shift_count_double(m, n, x, beta);

    // if k < N done
    if (on_t_k_star < n) {
        return obj_on_target_pos(on_t_k_star, m, n, x, beta);
    }

    // Compute optimal number of values shifted from O
    auto outside_t_k_star = opt_outside_shift_count(dtd, m, n, beta);

    return obj_out_target(dtd, outside_t_k_star, m, n, x, beta);
}

/**
 * Computes the smooth sensitivity contribution for a fixed edge when the objective function contributes positively.
 * Solves Problem (4) in the Appendix.
 */
double fix_edge_sens_pos(Graph &g,
                         Edge e,
                         const int lambda,
                         const std::list<int> &triangle_index_list,
                         const std::vector<Triangle> &triangles,
                         double beta,
                         double p,
                         bool inc) {
    // Computes the optimal value for a positive sum contribution
    double x = p / (std::pow(1.0 - p, 2));

    // build triangle values and target array
    std::vector<int> partial_triangle_weights = setup_partial_triangle_weights(g, e, triangle_index_list, triangles);

    if (partial_triangle_weights.empty()) {
        return 0;
    }

    auto [targets, zero_target] = setup_targets_unbiased(g, e, lambda, partial_triangle_weights, inc);

    // setup double target data structure
    DoubleTargetDistance dtd = DoubleTargetDistance(partial_triangle_weights, targets[0]);

    double opt = 0;

    for (int i = 0; i < targets.size(); ++i) {
        // for every target compute the optimal shift set
        double opt_t = opt_double_target_obj(dtd, x, beta) * std::exp(-beta * std::abs(targets[i] - zero_target));
        // check if new opt was found
        opt = std::max(opt, opt_t);
        // advance target
        if (i < targets.size() - 1) {
            dtd.update_target(targets[i + 1]);
        }
    }

    return opt;
}

// <-------- Negative Contribution To Sum ---------->

/**
 * Returns objective function (without fixed factor e^(-beta * |z_i|)) during the shifting of values from M
 */
double obj_on_target_neg(int k, int m, int n, double x, double beta) {
    // objective function during the step of assigning on target values for negative sum contribution
    double linear = (-(m - k) * x) + ((n + k) * (1.0 + 2.0 * x));
    return linear * std::exp(-beta * k);
}

/**
 * Computes the optimal number of values shifted from M
 */
int opt_on_target_shift_count_single(int m, int n, double x, double beta) {
    double k_star = 1.0 / beta + (-m * x + (1.0 + 2.0 * x) * n) / (1.0 + 3.0 * x);

    // Candidate integer values (clamped to [0, n])
    int k0 = std::clamp((int) std::floor(k_star), 0, m);
    int k1 = std::clamp(k0 + 1, 0, m);

    // Evaluate objective
    double v0 = obj_on_target_neg(k0, m, n, x, beta);
    double v1 = obj_on_target_neg(k1, m, n, x, beta);

    return (v1 > v0) ? k1 : k0;
}


/**
 * Computes the optimal objective value for a fixed single target
 */
double opt_single_target_obj(SingleTargetDistance dtd, double x, double beta) {
    // Compute |M|, |N|
    int m = dtd.adjacent_to_target_count();
    int n = dtd.on_target_count();

    // Compute optimal number of values connected from M
    int on_t_k_star = opt_on_target_shift_count_single(m, n, x, beta);

    // if k < M done
    if (on_t_k_star < m) {
        return obj_on_target_neg(on_t_k_star, m, n, x, beta);
    }

    // Compute optimal number of values shifted from O
    auto outside_t_k_star = opt_outside_shift_count(dtd, m, n, beta);

    return obj_out_target(dtd, outside_t_k_star, m, n, 1 + 2 * x, beta);
}

/**
 * Computes the smooth sensitivity contribution for a fixed edge when the objective function contributes negatively.
 * Solves Problem (5) in the Appendix.
 */
double fix_edge_sens_neg(Graph &g,
                         Edge e,
                         const int lambda,
                         const std::list<int> &triangle_index_list,
                         const std::vector<Triangle> &triangles,
                         double beta,
                         double p,
                         bool inc) {
    // Computes the optimal value for a negative sum contribution
    double x = p / std::pow(1.0 - p, 2);

    // build triangle values and target array
    std::vector<int> partial_triangle_weights = setup_partial_triangle_weights(g, e, triangle_index_list, triangles);

    if (partial_triangle_weights.empty()) {
        return 0;
    }

    auto [targets, zero_target] = setup_targets_unbiased(g, e, lambda, partial_triangle_weights, inc);

    // setup double target data structure
    SingleTargetDistance std = SingleTargetDistance(partial_triangle_weights, targets[0]);

    double opt = 0;

    for (int i = 0; i < targets.size(); ++i) {
        // for every target compute the optimal shift set
        double opt_t = opt_single_target_obj(std, x, beta) * std::exp(-beta * std::abs(targets[i] - zero_target));
        // check if new opt was found
        opt = std::max(opt, opt_t);
        // advance target
        if (i < targets.size() - 1) {
            std.update_target(targets[i + 1]);
        }
    }

    return opt;
}

double smooth_sensitivity_unb(Graph &g,
                              Node v,
                              const int lambda,
                              const std::list<int> &triangles_index_list,
                              const std::vector<Triangle> &triangles,
                              double beta,
                              double p) {
    double sens = 0;

#pragma omp parallel for reduction(max:sens)
    for (int i = 0; i < boost::degree(v, g); ++i) {
        // Fix the edge e = (v,u) and compute the maximum smooth sensitivity achieved by increasing / decreasing e
        Node u = boost::adjacent_vertices(v, g).first[i];
        Edge e = boost::edge(u, v, g).first;

        // Case b_i=1:
        double sens_inc_pos = fix_edge_sens_pos(g, e, lambda, triangles_index_list, triangles, beta, p, true);
        double sens_inc_neg = fix_edge_sens_neg(g, e, lambda, triangles_index_list, triangles, beta, p, true);

        // Case b_i=-1:
        double sens_dec_pos = fix_edge_sens_neg(g, e, lambda, triangles_index_list, triangles, beta, p, false);
        double sens_dec_neg = fix_edge_sens_pos(g, e, lambda, triangles_index_list, triangles, beta, p, false);

        double iter_max = std::max({
            sens_inc_pos,
            sens_inc_neg,
            sens_dec_pos,
            sens_dec_neg
        });

        sens = std::max(sens, iter_max);
    }

    return sens;
}
