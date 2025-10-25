#pragma once

#include "utils.h"

Graph generate_er_graph(size_t n, double p);

void assign_gaussian_weights(Graph &g, const double mu, const double std);

void add_discrete_laplace_noise(Graph &g, const double weight_eps);

Graph generate_graph(int num_vertices, double edge_probability, double weight_mu, double weight_std, double weight_eps);

Graph load_traffic_graph();

Graph load_telecom_graph();

Graph load_telecom_625_graph();

Graph load_telecom_400_graph();

Graph load_telecom_278_graph();