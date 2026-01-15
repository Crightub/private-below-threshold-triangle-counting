#include "../include/smooth_sens_utils.h"

#include <utility>
#include <vector>

#include "triangle.h"
#include "utils.h"

std::vector<int> setup_partial_triangle_weights(Graph &g,
                                                const Edge &e,
                                                const std::list<int> &triangle_index_list,
                                                const std::vector<Triangle> &triangles
) {
    int w_e = Triangle::get_edge_weight(g, e);
    std::vector<int> partial_triangle_weights;

    for (int t_id: triangle_index_list) {
        auto t = triangles[t_id];
        if (!t.contains_edge(g, e))
            continue;

        auto [noise_e1, w_e1, w_e2, w_e3] = t.get_triangle_weights(g);

        int partial_triangle_weight = noise_e1 + w_e1 + w_e2 + w_e3 - w_e;
        partial_triangle_weights.push_back(partial_triangle_weight);
    }

    std::sort(partial_triangle_weights.begin(), partial_triangle_weights.end());

    return partial_triangle_weights;
}

int setup_target(Graph &g,
                 const Edge &e,
                 const int lamba,
                 bool inc) {
    int w_e = Triangle::get_edge_weight(g, e);

    if (inc) {
        return lamba - 1 - w_e;
    } else {
        return lamba - w_e;
    }
}

std::vector<int> setup_targets_unbiased(Graph &g,
                                        const Edge &e,
                                        int lambda,
                                        const std::vector<int> &centers,
                                        bool inc) {
    std::vector<int> targets;

    int w_e = Triangle::get_edge_weight(g, e);

    // Add targets derived from centers
    for (int c: centers) {
        targets.push_back(c - 1.0);
        targets.push_back(c);
        targets.push_back(c + 1.0);
    }

    // Add lambda-based targets (z_i = 0, 1, -1)
    int t = inc ? lambda - w_e - 1 : lambda - w_e;
    targets.push_back(t - 1);
    targets.push_back(t);
    targets.push_back(t + 1);

    // Sort and remove duplicates
    std::sort(targets.begin(), targets.end());
    targets.erase(std::unique(targets.begin(), targets.end()), targets.end());

    return targets;
}
