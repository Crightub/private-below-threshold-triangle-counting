#include "../include/smooth_sens_utils.h"

#include <utility>
#include <vector>

#include "triangle.h"
#include "utils.h"

/**
 * Builds the vector c from Subsection A.1
 */
std::vector<int> setup_partial_triangle_weights(Graph &g,
                                                const Edge &e,
                                                const std::list<int> &triangle_index_list,
                                                const std::vector<Triangle> &triangles
) {
    int w_e = Triangle::get_edge_weight(g, e);
    std::vector<int> partial_triangle_weights;

    for (int t_id: triangle_index_list) {
        const Triangle t = triangles[t_id];
        if (!t.contains_edge(g, e))
            continue;

        auto [noise_e1, w_e1, w_e2, w_e3] = t.get_triangle_weights(g);
        int partial_triangle_weight = noise_e1 + w_e1 + w_e2 + w_e3 - w_e;
        partial_triangle_weights.push_back(partial_triangle_weight);
    }

    std::sort(partial_triangle_weights.begin(), partial_triangle_weights.end());
    return partial_triangle_weights;
}

std::pair<std::vector<int>, int> setup_targets_biased(Graph &g,
                                                      const Edge &e,
                                                      int lambda,
                                                      const std::vector<int> &centers,
                                                      bool inc) {
    std::vector<int> targets;

    int w_e = Triangle::get_edge_weight(g, e);

    // Add targets derived from centers
    targets.reserve(centers.size());
    for (int c: centers) {
        targets.push_back(c);
    }

    // Add lambda-based targets (z_i = 0, 1, -1)
    int t = inc ? lambda - w_e - 1 : lambda - w_e;
    targets.push_back(t);

    // Sort and remove duplicates
    std::sort(targets.begin(), targets.end());
    targets.erase(std::unique(targets.begin(), targets.end()), targets.end());

    return {targets, t};
}

/**
 * Setup set of possible targets T = {c_j + d :  c_j \in c, d \in {-1,0,1}} union {lambda - 1 - w_e}.
 */
std::pair<std::vector<int>, int> setup_targets_unbiased(Graph &g,
                                                        const Edge &e,
                                                        int lambda,
                                                        const std::vector<int> &centers,
                                                        bool inc) {
    std::vector<int> targets;

    int w_e = Triangle::get_edge_weight(g, e);

    // Add targets derived from centers
    for (int c: centers) {
        targets.push_back(c - 1);
        targets.push_back(c);
        targets.push_back(c + 1);
    }

    // Add lambda-based targets (z_i = 0, 1, -1)
    int t = inc ? lambda - w_e - 1 : lambda - w_e;
    targets.push_back(t - 1);
    targets.push_back(t);
    targets.push_back(t + 1);

    // Sort and remove duplicates
    std::sort(targets.begin(), targets.end());
    targets.erase(std::unique(targets.begin(), targets.end()), targets.end());

    return {targets, t};
}
