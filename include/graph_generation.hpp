#pragma once

#include "utils.hpp" 

Graph generateERGraph(size_t n, double p);
Graph assignGaussianWeights(Graph g, double mu, double std);
Graph addDiscreteLaplaceNoise(Graph g, double eps);
Graph generateGraph(int num_vertices, double edge_probability, double weight_mu, double weight_std, double weight_eps);