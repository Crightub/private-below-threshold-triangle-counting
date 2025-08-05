#pragma once

#include "utils.hpp"

struct PrivateCountingResult {
    double negative_triangle_count;
    double weight_noise;
    double count_noise;
};

int count_negative_triangles(const Graph &g);
PrivateCountingResult randomized_private_counting(const Graph &g, double eps, double eps2, bool use_weight_noise, bool use_count_noise,bool use_biased_estimator);