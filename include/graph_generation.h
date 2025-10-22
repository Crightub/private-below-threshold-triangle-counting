#pragma once

#include "utils.h"

Graph generate_er_graph(size_t n, double p);

void assign_gaussian_weights(Graph &g, const double mu, const double std);

void add_discrete_laplace_noise(Graph &g, const double weight_eps);

Graph generate_graph(int num_vertices, double edge_probability, double weight_mu, double weight_std, double weight_eps);

Graph load_snap_bitcoin_graph(const std::string &filename);

Graph load_snap_epinions_graph(const std::string &file_name);

Graph load_snap_wikipedia_graph(const std::string& file_name);

Graph load_snap_slashdot_graph(const std::string& file_name);

Graph load_traffic_graph(const std::string& file_name);

Graph load_telecomm_graph(const std::string& file_name);
