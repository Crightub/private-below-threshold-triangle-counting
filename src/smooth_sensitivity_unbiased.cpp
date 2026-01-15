#include "../include/smooth_sensitivity_unbiased.h"

#include "distance_retrieval.h"
#include "smooth_sens_utils.h"
#include "triangle.h"
#include "utils.h"


std::pair<bool, int> shift_set_can_be_increased(TargetDistance &td, int k, int m, int n, double beta) {
    auto [next_dis, l_idx] = td.k_th_closest_distance(k + 1);
    return {static_cast<double>(k + m + n) / (k + m + n + 1) < std::exp(-beta * next_dis), l_idx};
}

/**
 *  @brief Computes the optimal number of values shifted from the outside set (O) to the target
 *  This function performs a binary search over the number of elements
 *  shifted outside the target set. For a given shift count `k`, it checks
 *  whether the shift can still be increased based on the objective value.
 *
 *  @complexity O(log^2 d) where d is the number of outside targets
 */
std::pair<int, int> opt_outside_shift_count(TargetDistance &td, int m, int n, double beta) {
    int lo = 0;
    int hi = td.outside_target_count();

    int opt_l = 0;
    bool can_be_increased;

    while (lo < hi) {
        int mid = (lo + hi) / 2;
        int l_idx;
        std::tie(can_be_increased, l_idx) = shift_set_can_be_increased(td, mid, m, n, beta);
        if (can_be_increased) {
            lo = mid + 1;
            opt_l = l_idx;
        } else {
            hi = mid;
        }
    }

    return {lo, opt_l};
}

double obj_on_target_pos(int k, int m, int n, double x, double beta) {
    // objective function during the step of assigning on target values for positive sum contribution
    double linear = ((m + k) * x) + ((n - k) * (-1.0 - 2.0 * x));
    return linear * std::exp(-beta * k);
}

double obj_on_target_neg(int k, int m, int n, double x, double beta) {
    // objective function during the step of assigning on target values for negative sum contribution
    double linear = (-(m - k) * x) + ((n + k) * (1.0 + 2.0 * x));
    return linear * std::exp(-beta * k);
}


double obj_out_target(
    TargetDistance &dtd,
    int k,
    int l_idx,
    int m,
    int n,
    double f,
    double beta
) {
    double sum_k = dtd.sum_k_closest_distances(k, l_idx);
    double linear = f * (k + n + m);
    double exponent = -beta * (n + sum_k);
    return linear * std::exp(exponent);
}

int opt_on_target_shift_count(int m, int n, double x, double beta, double pos) {
    double y = pos ? m * x + (-1.0 - 2.0 * x) * n : -m * x + (1.0 + 2.0 * x) * n;
    double k_star = 1.0 / beta + y / (1.0 + 3.0 * x);

    // Candidate integer values (clamped to [0, n])
    int k0 = std::clamp((int) std::floor(k_star), 0, n);
    int k1 = std::clamp(k0 + 1, 0, n);

    // Evaluate objective
    double v0 = obj_on_target_pos(k0, m, n, x, beta);
    double v1 = obj_on_target_pos(k1, m, n, x, beta);

    return (v1 > v0) ? k1 : k0;
}

double opt_fixed_target_obj(DoubleTargetDistance dtd, double x, double beta) {
    // Compute |M|, |N|
    int m = dtd.adjacent_to_target_count();
    int n = dtd.on_target_count();

    // Compute optimal number of values connected from N
    int on_t_k_star = opt_on_target_shift_count(m, n, x, beta, true);

    // if k < N done
    if (on_t_k_star < n) {
        return obj_on_target_pos(on_t_k_star, m, n, x, beta);
    }

    // Compute optimal number of values shifted from O
    auto [outside_t_k_star, l_idx] = opt_outside_shift_count(dtd, m, n, beta);

    return obj_out_target(dtd, outside_t_k_star, l_idx, m, n, x, beta);
}


/**
 * @brief Computes the smooth sensitivity contribution for a fixed edge
 *        when the objective function contributes positively.
 *
 * This function evaluates the smooth sensitivity of the triangle-counting
 * statistic with respect to a fixed edge `e`, restricted to configurations
 * where the contribution of the objective function is positive.
 *
 * @complexity O(d + d * log^2 d)
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
    double x = p / std::pow(1.0 - p, 2);

    // build triangle values and target array
    std::vector<int> partial_triangle_weights = setup_partial_triangle_weights(g, e, triangle_index_list, triangles);
    std::vector<int> targets = setup_targets_unbiased(g, e, lambda, partial_triangle_weights, inc);

    // setup double target data structure
    DoubleTargetDistance dtd = DoubleTargetDistance(partial_triangle_weights, targets[0]);

    double opt;

    for (int i = 0; i < targets.size(); ++i) {
        // for every target compute the optimal shift set
        double opt_t = opt_fixed_target_obj(dtd, x, beta);
        // check if new opt was found
        opt = std::max(opt, opt_t);
        // advance target
        if (i < targets.size() - 1) {
            dtd.update_target(targets[i + 1]);
        }
    }

    return opt;
}

double obj_out_target_single(
    DoubleTargetDistance &dtd,
    int k,
    int l_idx,
    int m,
    int n,
    double x,
    double beta
) {
    double sum_k = dtd.sum_k_closest_distances(k, l_idx);
    double linear = (1 + 2 + x) * (k + n + m);
    double exponent = -beta * (n + sum_k);
    return linear * std::exp(exponent);
}

double opt_fixed_target_single_obj(SingleTargetDistance dtd, double x, double beta) {
    // Compute |M|, |N|
    int m = dtd.adjacent_to_target_count();
    int n = dtd.on_target_count();

    // Compute optimal number of values connected from N
    int on_t_k_star = opt_on_target_shift_count(m, n, x, beta, false);

    // if k < N done
    if (on_t_k_star < n) {
        return obj_on_target_neg(on_t_k_star, m, n, x, beta);
    }

    // Compute optimal number of values shifted from O
    auto [outside_t_k_star, l_idx] = opt_outside_shift_count(dtd, m, n, beta);

    return obj_out_target(dtd, outside_t_k_star, l_idx, m, n, 1 + 2 * x, beta);
}

double fixed_edge_sens_neg(Graph &g,
                           Edge e,
                           const int lambda,
                           const std::list<int> &triangle_index_list,
                           const std::vector<Triangle> &triangles,
                           double beta,
                           double eps,
                           bool inc) {
    // Computes the optimal value for a negative sum contribution
    double p = std::exp(-eps);
    double x = p / std::pow(1.0 - p, 2);

    // build triangle values and target array
    std::vector<int> partial_triangle_weights = setup_partial_triangle_weights(g, e, triangle_index_list, triangles);
    std::vector<int> targets = setup_targets_unbiased(g, e, lambda, partial_triangle_weights, inc);

    // setup double target data structure
    SingleTargetDistance std = SingleTargetDistance(partial_triangle_weights, targets[0]);

    double opt;

    for (int i = 0; i < targets.size(); ++i) {
        // for every target compute the optimal shift set
        double opt_t = opt_fixed_target_single_obj(std, x, beta);
        // check if new opt was found
        opt = std::max(opt, opt_t);
        // advance target
        if (i < targets.size() - 1) {
            std.update_target(targets[i + 1]);
        }
    }

    return opt;
}

double smooth_sensitivity(Graph &g,
                          Node v,
                          const int lambda,
                          const std::list<int> &triangles_index_list,
                          const std::vector<Triangle> &triangles,
                          double beta,
                          double p) {
    double sens = 0;

    // #pragma omp parallel for reduction(max:sens)
    for (int i = 0; i < boost::degree(v, g); ++i) {
        // Fix the edge e = (v,u) and compute the maximum smooth sensitivity achieved by increasing / decreasing e
        Node u = boost::adjacent_vertices(v, g).first[i];
        Edge e = boost::edge(u, v, g).first;

        double sens_inc_pos = fix_edge_sens_pos(g, e, lambda, triangles_index_list, triangles, beta, p, true);
        double sens_inc_neg = fix_edge_sens_pos(g, e, lambda, triangles_index_list, triangles, beta, p, true);

        std::cout << "Edge e: sens_inc_pos: " << sens_inc_pos << ", sens_inc_neg: " << sens_inc_neg << std::endl;
        // std::cout << "Edge " << e << ", sens_inc: " << sens_inc << ", sens_dec: " << sens_dec << std::endl;

        double fixed_edge_sens = std::max(sens_inc_pos, sens_inc_neg);
        sens = std::max(fixed_edge_sens, sens);
    }

    return sens;
}
