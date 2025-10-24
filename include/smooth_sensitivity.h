#pragma once

#include "treap.h"
#include "triangle.h"

struct PrivateCountingConfig;

double find_k_smallest_distance(TreapNode *left, TreapNode *right, int k, double offL, double offR);

int optimal_shift_set_size(TreapNode *left, TreapNode *right, double offL,
                           double offR, double beta);

double smooth_sensitivity(Graph &g, Node v, int lambda, std::list<int> &triangle_index_list, const std::vector<Triangle> &triangle, double beta);

void apply_smooth_sensitivity(Graph &g,
                              const PrivateCountingConfig &cfg,
                              std::vector<TriangleCount> &counts,
                              std::vector<std::list<int>> &node_triangle_map,
			      const std::vector<Triangle> &triangles);

double compute_sum(TreapNode *left, TreapNode *right, double offL, double offR, int l, int r);

std::pair<TreapNode *, TreapNode *> advance_center(TreapNode *left,
                                                   TreapNode *right,
                                                   double &offL,
                                                   double &offR,
                                                   const std::vector<double> &centers,
                                                   int j);
