#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "distribution.h"
#include "utils.h"
#include "graph_generation.h"
#include <fstream>

Graph generate_er_graph(size_t n, double p)
{
    Graph g(ERGen(rng, n, p), ERGen(), n);
    return g;
}

Graph assign_gaussian_weights(Graph g, double mu, double std)
{
    boost::normal_distribution<> nd(mu, std);
    boost::variate_generator<boost::minstd_rand, boost::normal_distribution<>> var_n(rng, nd);
    boost::graph_traits<Graph>::edge_iterator ei, ei_end;

    for (boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei)
    {
        int weight = static_cast<int>(std::round(var_n()));
        put(&edge_info::weight, g, *ei, weight);
    }

    return g;
}

Graph add_discrete_laplace_noise(Graph g, double eps)
{   
    double p = std::exp(-eps);

    double unbias_control = 0;
    boost::graph_traits<Graph>::edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei)
    {
        int noise = sample_discrete_laplace(p);
        put(&edge_info::noise, g, *ei, noise);
        unbias_control += unbiased_estimator(noise, p, 0);
    }

    return g;
}

Graph generate_graph(int num_vertices, double edge_probability, double weight_mu, double weight_std, double weight_eps){
    Graph g = generate_er_graph(num_vertices, edge_probability);
    g = assign_gaussian_weights(g, weight_mu, weight_std);
    Graph noisyGraph = add_discrete_laplace_noise(g, weight_eps);

    return noisyGraph;
}

Graph load_snap_bitcoin_graph(const std::string &file_name) {
    Graph g;

    std::ifstream file(file_name);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file " + file_name);
    }

    std::string line;

    while (std::getline(file, line)) {
        if (line.empty()) continue;

        std::istringstream iss(line);
        std::string token;

        int source, target, rating;
        double time;

        std::getline(iss, token, ',');
        source = std::stoi(token);
        std::getline(iss, token, ',');
        target = std::stoi(token);
        std::getline(iss, token, ',');
        rating = std::stoi(token);
        std::getline(iss, token, ',');
        time = std::stod(token);

        if (std::max(source, target) >= boost::num_vertices(g)) {
            while (boost::num_vertices(g) <= std::max(source, target)) {
                boost::add_vertex(g);
            }
        }

        auto e = boost::add_edge(source, target, g);
        if (e.second) {
            g[e.first].weight = rating;
            g[e.first].load = 0;
            g[e.first].noise = 0;
        }
    }

    return g;
}

Graph load_snap_epinions_graph(const std::string& file_name) {
    Graph g;

    std::ifstream file(file_name);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file " + file_name);
    }

    std::string line;

    while (std::getline(file, line)) {
        if (line.empty()) continue;

        std::istringstream iss(line);
        std::string token;

        int source, target, rating;

        std::getline(iss, token, '\t');
        source = std::stoi(token);
        std::getline(iss, token, '\t');
        target = std::stoi(token);
        std::getline(iss, token, '\t');
        rating = std::stoi(token);

        if (std::max(source, target) >= boost::num_vertices(g)) {
            while (boost::num_vertices(g) <= std::max(source, target)) {
                boost::add_vertex(g);
            }
        }
        auto e = boost::add_edge(source, target, g);
        if (e.second) {
            g[e.first].weight = rating;
            g[e.first].load = 0;
            g[e.first].noise = 0;
        }
    }

    return g;
}

Graph load_snap_wikipedia_graph(const std::string& file_name) {
    // TODO: Implement
    return Graph(0);
}

Graph load_snap_slashdot_graph(const std::string& file_name) {
    // TODO: Implement
    return Graph(0);
}