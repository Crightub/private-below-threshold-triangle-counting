#pragma once

#include <vector>
#include <iostream>
#include <random>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/log/trivial.hpp>
#include <sstream>

// <-- Graph Utils -->

struct edge_info {
    int weight = 0;
    int noise = 0;
    int load = 0;
    int involved_triangles = 0;
};

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, edge_info> Graph;
typedef boost::sorted_erdos_renyi_iterator<boost::minstd_rand, Graph> ERGen;

typedef boost::graph_traits<Graph>::edge_descriptor Edge;
typedef boost::graph_traits<Graph>::vertex_descriptor Node;
typedef boost::graph_traits<Graph>::adjacency_iterator AdjIter;

extern boost::minstd_rand rng;

std::vector<Node> neighbors(const Graph &g, Node &v);

double unbiased_estimator(int x, double p, int lambda);

int biased_estimator(int w1, int w2, int w3, int lambda);

inline bool is_same_edge(const Graph &g, Edge e1, Edge e2) {
    auto s1 = boost::source(e1, g);
    auto t1 = boost::target(e1, g);

    auto s2 = boost::source(e2, g);
    auto t2 = boost::target(e2, g);

    if (s1 == s2 && t1 == t2) {
        return true;
    }

    if (s1 == t2 && t1 == s2) {
        return true;
    }

    return false;
}

struct TriangleCount {
    long naive;
    long opt;
    double unbiased;
    double biased;
};

// <-- Debug function -->
inline long long compute_c4_instances(const Graph &g) {
    long long sum = 0;
    int max_load = 0;
    auto load_map = get(&edge_info::load, g);

    for (auto [e_it, e_it_end] = boost::edges(g); e_it != e_it_end; ++e_it) {
        Edge e = *e_it;
        int l = get(load_map, e);
        if (l >= 2) {
            sum += (static_cast<long long>(l) * (l - 1)) / 2;
        }
        max_load = std::max(max_load, l);
    }

    std::cout << "Max Load: " << max_load << std::endl;
    return sum;
}

// <-- Print functions -->
inline void print_vector(const std::vector<double> &v) {
    std::ostringstream oss;
    oss << "[";
    for (int j = 0; j < v.size(); j++) {
        oss << v[j];
        if (j < v.size() - 1) {
            oss << ", ";
        }
    }
    oss << "]";
    std::cout << oss.str() << std::endl;
}


inline void print_matrix(const std::vector<std::vector<double> > &A) {
    std::ostringstream oss;
    for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A[i].size(); ++j) {
            oss << A[i][j];
            if (j < A[i].size() - 1) {
                oss << "\t";
            }
        }
        if (i < A.size() - 1) {
            oss << "\n";
        }
    }
}
