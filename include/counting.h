#pragma once

#include "utils.h"
#include "triangle.h"
#include <nlohmann/json.hpp>

using nlohmann::json;

struct PrivateCountingResult {
    int opt;
    double global_unbiased;
    double global_biased;
    double smooth_unbiased;
    double smooth_biased;
};

inline void to_json(json &j, const PrivateCountingResult &res) {
    j = {
        {"opt", res.opt},
        {"global_unbiased", res.global_unbiased},
        {"global_biased", res.global_biased},
        {"smooth_unbiased", res.smooth_unbiased},
        {"smooth_biased", res.smooth_biased}
    };
}

struct PrivateCountingConfig {
    double weight_eps;
    double count_eps;
    int lambda;
    bool use_load_balancing;
    bool use_smooth_sensitivity;
};

PrivateCountingResult private_counting(Graph &g, PrivateCountingConfig &base_cfg,
                                       const std::list<Triangle> *ptr_triangles);

void count_local_negative_triangles(Graph &g,
                                    const PrivateCountingConfig &base_cfg,
                                    std::vector<TriangleCount> &counts,
                                    std::vector<std::list<Triangle> > &node_triangle_map,
                                    std::list<Triangle> &triangles
);