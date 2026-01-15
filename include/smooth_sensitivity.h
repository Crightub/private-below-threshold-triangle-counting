#pragma once

#include "treap.h"
#include "triangle.h"

struct PrivateCountingConfig;

double smooth_sensitivity(Graph &g,
                          Node v,
                          int lambda,
                          const std::list<int> &triangles_index_list,
                          const std::vector<Triangle> &triangles,
                          double beta);

void apply_smooth_sensitivity(Graph &g,
                              const PrivateCountingConfig &cfg,
                              std::vector<TriangleCount> &counts,
                              std::vector<std::list<int> > &node_triangle_map,
                              const std::vector<Triangle> &triangles);
