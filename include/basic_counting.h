#pragma once

#include "utils.h"

struct PrivateCountingResult {
    double negative_triangle_count;
    double weight_noise;
    double count_noise;
};


struct PrivateCountingConfig {
    bool use_unbiased_estimator;
    double weight_eps;
    double count_eps;
    bool use_smooth_sensitivity;
};

int count_negative_triangles(const Graph &g);
PrivateCountingResult randomized_private_counting(const Graph &g, const PrivateCountingConfig &cfg);
void setup_graph_matrices(const Graph &g,
                          std::vector<std::vector<bool> > &A,
                          std::vector<std::vector<int> > &W,
                          std::vector<std::vector<int> > &N);