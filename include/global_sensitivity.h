#pragma once

#include <vector>

#include "triangle.h"
#include "utils.h"

struct PrivateCountingConfig;

int global_sensitivity(const Graph &g, int i, const std::list<Triangle> &triangles);

void apply_global_sensitivity(const Graph &g,
                              const PrivateCountingConfig &cfg,
                              std::vector<TriangleCount> &counts,
                              const std::vector<std::list<Triangle> > &node2TriangleMap);
