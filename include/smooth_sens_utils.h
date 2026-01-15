#ifndef NEGATIVE_TRIANGLE_COUNTING_SMOOTH_SENS_UTILS_H
#define NEGATIVE_TRIANGLE_COUNTING_SMOOTH_SENS_UTILS_H
#include <vector>

#include "triangle.h"
#include "utils.h"

std::vector<int> setup_partial_triangle_weights(Graph &g,
                                                const Edge &e,
                                                const std::list<int> &triangle_index_list,
                                                const std::vector<Triangle> &triangles);

int setup_target(Graph &g,
                 const Edge &e,
                 int lamba,
                 bool inc);

std::vector<int> setup_targets_unbiased(Graph &g,
                                        const Edge &e,
                                        int lambda,
                                        const std::vector<int> &centers,
                                        bool inc);

#endif //NEGATIVE_TRIANGLE_COUNTING_SMOOTH_SENS_UTILS_H
