#ifndef NEGATIVE_TRIANGLE_COUNTING_SMOOTH_SENSITIVITY_UNBIASED_H
#define NEGATIVE_TRIANGLE_COUNTING_SMOOTH_SENSITIVITY_UNBIASED_H
#include "triangle.h"
#include "utils.h"


double smooth_sensitivity_unb(Graph &g,
                              Node v,
                              int lambda,
                              const std::list<int> &triangles_index_list,
                              const std::vector<Triangle> &triangles,
                              double beta,
                              double p);

double fix_edge_sens_pos(Graph &g,
                         Edge e,
                         int lambda,
                         const std::list<int> &triangle_index_list,
                         const std::vector<Triangle> &triangles,
                         double beta,
                         double p,
                         bool inc);


double fix_edge_sens_neg(Graph &g,
                         Edge e,
                         int lambda,
                         const std::list<int> &triangle_index_list,
                         const std::vector<Triangle> &triangles,
                         double beta,
                         double p,
                         bool inc);

#endif //NEGATIVE_TRIANGLE_COUNTING_SMOOTH_SENSITIVITY_UNBIASED_H
