#ifndef NEGATIVE_TRIANGLE_COUNTING_SMOOTH_SENSITIVITY_H
#define NEGATIVE_TRIANGLE_COUNTING_SMOOTH_SENSITIVITY_H
#include "treap.h"

struct PrivateCountingConfig;

std::list<double> compute_partial_triangle_weights(Graph g, int v, int i);

double fixed_edge_sensitivity(const Graph &g, int v, int i, double beta, bool inc);

double find_k_smallest_distance(TreapNode *left, TreapNode *right, int k, double offL, double offR);

int optimal_shift_set_size(TreapNode *left, TreapNode *right, double offL,
                           double offR, double beta);

double smooth_sensitivity(const Graph &g, int v, double beta);

void apply_smooth_sensitivity(const Graph &g, const PrivateCountingConfig &cfg, std::vector<double> &counts);

double compute_sum(TreapNode *left, TreapNode *right, double offL, double offR, int l, int r);

std::pair<TreapNode *, TreapNode *> advance_center(TreapNode *left,
                                                   TreapNode *right,
                                                   double &offL,
                                                   double &offR,
                                                   const std::vector<double> &centers,
                                                   int j);

#endif //NEGATIVE_TRIANGLE_COUNTING_SMOOTH_SENSITIVITY_H
