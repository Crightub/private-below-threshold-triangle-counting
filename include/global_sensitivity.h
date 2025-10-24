#pragma once

#include <vector>

#include "triangle.h"
#include "utils.h"

struct PrivateCountingConfig;

int global_sensitivity(const Graph &g, int i, const std::list<int> &node_triangle_index_list, const std::vector<Triangle> &triangles);

void apply_global_sensitivity(const Graph &g,
                              const PrivateCountingConfig &cfg,
                              std::vector<TriangleCount> &counts,
                              const std::vector<std::list<int>> &node_triangle_index_map,
			      const std::vector<Triangle> &triangles);
