#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "distribution.h"
#include "utils.h"
#include "graph_generation.h"
#include <fstream>

Graph generate_er_graph(size_t n, double p) {
    Graph g(ERGen(rng, n, p), ERGen(), n);
    return g;
}

void assign_gaussian_weights(Graph &g, const double mu, const double std) {
    boost::normal_distribution<> nd(mu, std);
    boost::variate_generator<boost::minstd_rand, boost::normal_distribution<> > var_n(rng, nd);

    for (auto [ei, ei_end] = boost::edges(g); ei != ei_end; ++ei) {
        int weight = static_cast<int>(std::round(var_n()));
        put(&edge_info::weight, g, *ei, weight);
    }
}

void add_discrete_laplace_noise(Graph &g, const double weight_eps, const double count_eps) {
    double p = std::exp(-weight_eps);

    double triv_p = std::exp(-(weight_eps + count_eps));

    for (auto [ei, ei_end] = boost::edges(g); ei != ei_end; ++ei) {
        int noise = sample_discrete_laplace(p);
        put(&edge_info::noise, g, *ei, noise);

        int triv_noise = sample_discrete_laplace(triv_p);
        put(&edge_info::triv_noise, g, *ei, triv_noise);
    }
}

//
// Graph load_snap_bitcoin_graph(const std::string &file_name) {
//     Graph g;
//
//     std::ifstream file(file_name);
//     if (!file.is_open()) {
//         throw std::runtime_error("Could not open file " + file_name);
//     }
//
//     std::string line;
//
//     std::map<std::pair<int, int>, int> edge_ratings;
//
//     while (std::getline(file, line)) {
//         if (line.empty()) continue;
//
//         std::istringstream iss(line);
//         std::string token;
//
//         int source, target, rating;
//         double time;
//
//         std::getline(iss, token, ',');
//         source = std::stoi(token);
//         std::getline(iss, token, ',');
//         target = std::stoi(token);
//         std::getline(iss, token, ',');
//         rating = std::stoi(token);
//         std::getline(iss, token, ',');
//         time = std::stod(token);
//
//         if (std::max(source, target) >= boost::num_vertices(g)) {
//             while (boost::num_vertices(g) <= std::max(source, target)) {
//                 boost::add_vertex(g);
//             }
//         }
//
//         auto key = std::make_pair(std::min(source, target), std::max(target, source));
//
//         if (edge_ratings.find(key) != edge_ratings.end()) {
//             int reverse_rating = edge_ratings[key];
//             int avg_rating = static_cast<int>(std::round((rating + reverse_rating) / 2.0));
//
//             auto e1 = boost::edge(std::max(source, target), std::min(target, source), g);
//             if (e1.second) {
//                 g[e1.first].weight = avg_rating;
//                 g[e1.first].load = 0;
//                 g[e1.first].noise = 0;
//             }
//         } else {
//             edge_ratings[key] = rating;
//
//             auto e = boost::add_edge(std::max(source, target), std::min(target, source), g);
//             if (e.second) {
//                 g[e.first].weight = rating;
//                 g[e.first].load = 0;
//                 g[e.first].noise = 0;
//             }
//         }
//     }
//
//     return g;
// }

void ensure_vertices(Graph& g, int u, int v) {
    while (boost::num_vertices(g) <= std::max(u, v)) {
        boost::add_vertex(g);
    }
}

void add_edge_with_weight(Graph& g, int u, int v, int weight) {
    auto e = boost::add_edge(u, v, g);
    if (e.second) {
        g[e.first].weight = weight;
        g[e.first].load = 0;
        g[e.first].noise = 0;
    }
}

Graph load_weighted_graph(const std::string& file_name) {
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
        int values[3];  // source, target, weight

        for (int i = 0; i < 2; ++i) {
            std::getline(iss, token, ',');
            values[i] = std::stoi(token);
        }

        std::getline(iss, token, ',');
        values[2] = std::stoi(token);

        ensure_vertices(g, values[0], values[1]);
        add_edge_with_weight(g, values[0], values[1], values[2]);
    }

    return g;
}

Graph load_traffic_graph() {
    return load_weighted_graph("../data/istanbul_traffic.csv");
}


Graph load_telecom_graph() {
    return load_weighted_graph("../data/milan_telecom.csv");
}

Graph load_telecom_625_graph() {
    return load_weighted_graph("../data/MItoMI-2013-11-03-G-625.csv");
}

Graph load_telecom_400_graph() {
    return load_weighted_graph("../data/MItoMI-2013-11-03-G-400.csv");
}

Graph load_telecom_278_graph() {
    return load_weighted_graph("../data/MItoMI-2013-11-03-G-278.csv");
}