#ifndef NEGATIVE_TRIANGLE_COUNTING_GRAPH_PRESETS_H
#define NEGATIVE_TRIANGLE_COUNTING_GRAPH_PRESETS_H

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/subgraph.hpp>

#include "utils.h"

inline void assign_weights_and_noise(
    Graph &g,
    const std::map<std::pair<int, int>, int> &weights,
    const std::map<std::pair<int, int>, int> &noise
) {
    boost::graph_traits<Graph>::edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei) {
        int u = boost::source(*ei, g);
        int v = boost::target(*ei, g);

        auto key = std::minmax(u, v);

        if (weights.count(key)) {
            put(&edge_info::weight, g, *ei, weights.at(key));
        }
        if (noise.count(key)) {
            put(&edge_info::noise, g, *ei, noise.at(key));
        }
    }
}


inline Graph empty_graph(int n) {
    return Graph(n);
}

inline Graph complete_graph(int n) {
    Graph g(n);
    for (int u = 0; u < n; ++u) {
        for (int v = u + 1; v < n; ++v) {
            boost::add_edge(u, v, g);
        }
    }

    return g;
}

inline Graph complete_bipartite_graph(int m, int n) {
    Graph g(m + n);
    for (int u = 0; u < m; ++u) {
        for (int v = m; v < m + n; ++v) {
            boost::add_edge(u, v, g);
        }
    }
    return g;
}

inline Graph cycle_graph(int n) {
    Graph g(n);
    if (n < 2) return g;

    for (int u = 0; u < n; ++u) {
        int v = (u + 1) % n;
        boost::add_edge(u, v, g);
    }
    return g;
}

inline Graph path_graph(int n) {
    Graph g(n);
    for (int u = 0; u < n - 1; ++u) {
        boost::add_edge(u, u + 1, g);
    }
    return g;
}

inline Graph star_graph(int n) {
    Graph g(n);
    if (n == 0) return g;

    for (int v = 1; v < n; ++v) {
        boost::add_edge(0, v, g);
    }
    return g;
}

inline Graph star_outline_graph(int n) {
    Graph g(n);
    if (n == 0) return g;

    for (int v = 1; v < n; ++v) {
        boost::add_edge(0, v, g);
    }

    for (int v = 1; v < n - 1; ++v) {
        boost::add_edge(v, v + 1, g);
    }

    boost::add_edge(1, n - 1, g);
    return g;
}

inline Graph triangle_graph(int n) {
    Graph g(n);
    if (n < 2) return g;

    boost::add_edge(0, 1, g);

    for (int v = 2; v < n; ++v) {
        boost::add_edge(0, v, g);
        boost::add_edge(1, v, g);
    }

    return g;
}

inline Graph distinct_triangle_graph(int num_triangles) {
    // non-overlapping triangles all connected to a root
    Graph g(2 * num_triangles + 1);

    for (int i = 1; i < 2 * num_triangles + 1; i += 2) {
        boost::add_edge(0, i, g);
        boost::add_edge(0, i + 1, g);
        boost::add_edge(i, i + 1, g);
    }

    return g;
}

inline Graph example_graph_1() {
    Graph g = triangle_graph(8);

    std::map<std::pair<int, int>, int> weights = {
        {{0, 1}, 0},
        {{0, 2}, -5},
        {{0, 3}, -3},
        {{0, 4}, -1},
        {{0, 5}, 3},
        {{0, 6}, 4},
        {{0, 7}, 5},
        {{1, 2}, 0},
        {{1, 3}, 0},
        {{1, 4}, 0},
        {{1, 5}, 0},
        {{1, 6}, 0},
        {{1, 7}, 0}
    };

    std::map<std::pair<int, int>, int> noise = {
        {{0, 1}, 0},
        {{0, 2}, 0},
        {{0, 3}, 0},
        {{0, 4}, 0},
        {{0, 5}, 0},
        {{0, 6}, 0},
        {{0, 7}, 0},
        {{1, 2}, 0},
        {{1, 3}, 0},
        {{1, 4}, 0},
        {{1, 5}, 0},
        {{1, 6}, 0},
        {{1, 7}, 0}
    };

    assign_weights_and_noise(g, weights, noise);

    return g;
}

inline Graph example_graph_2() {
    Graph g = star_outline_graph(6);

    std::map<std::pair<int, int>, int> weights = {
        {{0, 1}, 0},
        {{0, 2}, 1},
        {{0, 3}, -1},
        {{0, 4}, -1},
        {{0, 5}, -1},
        {{1, 2}, -1},
        {{2, 3}, -1},
        {{3, 4}, -1},
        {{4, 5}, -1},
        {{1, 5}, -1},
    };

    std::map<std::pair<int, int>, int> noise = {
        {{0, 1}, 0},
        {{0, 2}, 0},
        {{0, 3}, 0},
        {{0, 4}, 0},
        {{0, 5}, 0},
        {{1, 2}, 0},
        {{2, 3}, 0},
        {{3, 4}, 0},
        {{4, 5}, 0},
        {{1, 5}, 0},
    };

    assign_weights_and_noise(g, weights, noise);

    return g;
}

inline Graph example_graph_2_scaled() {
    // Same graph structure as example_graph_2 but with weights increased by 5
    Graph g = star_outline_graph(6);

    std::map<std::pair<int, int>, int> weights = {
        {{0, 1}, 5},
        {{0, 2}, 6},
        {{0, 3}, 4},
        {{0, 4}, 4},
        {{0, 5}, 4},
        {{1, 2}, 4},
        {{2, 3}, 4},
        {{3, 4}, 4},
        {{4, 5}, 4},
        {{1, 5}, 4},
    };

    std::map<std::pair<int, int>, int> noise = {
        {{0, 1}, 0},
        {{0, 2}, 0},
        {{0, 3}, 0},
        {{0, 4}, 0},
        {{0, 5}, 0},
        {{1, 2}, 0},
        {{2, 3}, 0},
        {{3, 4}, 0},
        {{4, 5}, 0},
        {{1, 5}, 0},
    };

    assign_weights_and_noise(g, weights, noise);

    return g;
}

inline Graph example_graph_3() {
    Graph g = distinct_triangle_graph(3);

    std::map<std::pair<int, int>, int> weights = {
        {{0, 1}, -1},
        {{0, 2}, -1},
        {{0, 3}, -1},
        {{0, 4}, 0},
        {{0, 5}, 1},
        {{0, 6}, -1},
        {{1, 2}, 2},
        {{3, 4}, 0},
        {{5, 6}, -1}
    };

    std::map<std::pair<int, int>, int> noise = {
        {{0, 1}, 0},
        {{0, 2}, 0},
        {{0, 3}, 0},
        {{0, 4}, 0},
        {{0, 5}, 0},
        {{0, 6}, 0},
        {{1, 2}, 0},
        {{3, 4}, 0},
        {{5, 6}, 0}
    };

    assign_weights_and_noise(g, weights, noise);

    return g;
}


inline Graph example_graph_4() {
    Graph g = star_outline_graph(6);

    std::map<std::pair<int, int>, int> weights = {
        {{0, 1}, -100},
        {{0, 2}, -100},
        {{0, 3}, -100},
        {{0, 4}, -100},
        {{0, 5}, -100},
        {{1, 2}, -100},
        {{2, 3}, -100},
        {{3, 4}, -100},
        {{4, 5}, -100},
        {{1, 5}, -100},
    };

    std::map<std::pair<int, int>, int> noise = {
        {{0, 1}, 0},
        {{0, 2}, 0},
        {{0, 3}, 0},
        {{0, 4}, 0},
        {{0, 5}, 0},
        {{1, 2}, 0},
        {{2, 3}, 0},
        {{3, 4}, 0},
        {{4, 5}, 0},
        {{1, 5}, 0},
    };

    assign_weights_and_noise(g, weights, noise);

    return g;
}

inline Graph example_graph_5() {
    Graph g(11);

    boost::add_edge(0, 1, g);
    boost::add_edge(0, 2, g);
    boost::add_edge(0, 3, g);
    boost::add_edge(0, 4, g);
    boost::add_edge(0, 5, g);
    boost::add_edge(0, 6, g);

    boost::add_edge(1, 2, g);
    boost::add_edge(1, 7, g);
    boost::add_edge(1, 8, g);
    boost::add_edge(1, 9, g);
    boost::add_edge(1, 10, g);

    boost::add_edge(3, 6, g);
    boost::add_edge(4, 6, g);
    boost::add_edge(5, 6, g);

    boost::add_edge(7, 10, g);
    boost::add_edge(8, 10, g);
    boost::add_edge(9, 10, g);

    std::map<std::pair<int, int>, int> weights = {
        {{0, 1}, -1},
        {{0, 2}, -1},
        {{0, 3}, -1},
        {{0, 4}, -1},
        {{0, 5}, -1},
        {{0, 6}, -1},
        {{1, 2}, -1},
        {{1, 7}, -1},
        {{1, 8}, -1},
        {{1, 9}, -1},
        {{1, 10}, -1},
        {{3, 6}, -1},
        {{4, 6}, -1},
        {{5, 6}, -1},
        {{7, 10}, -1},
        {{8, 10}, -1},
        {{9, 10}, -1},
    };

    std::map<std::pair<int, int>, int> noise = {
        {{0, 1}, 0},
        {{0, 2}, 0},
        {{0, 3}, 0},
        {{0, 4}, 0},
        {{0, 5}, 0},
        {{0, 6}, 0},
        {{1, 2}, 0},
        {{1, 7}, 0},
        {{1, 8}, 0},
        {{1, 9}, 0},
        {{1, 10}, 0},
        {{3, 6}, 0},
        {{4, 6}, 0},
        {{5, 6}, 0},
        {{7, 10}, 0},
        {{8, 10}, 0},
        {{9, 10}, 0},
    };

    assign_weights_and_noise(g, weights, noise);

    return g;
}

inline Graph example_graph_6() {
    Graph g(6);

    boost::add_edge(0, 1, g);
    boost::add_edge(0, 2, g);
    boost::add_edge(0, 3, g);
    boost::add_edge(0, 4, g);
    boost::add_edge(0, 5, g);

    boost::add_edge(1, 2, g);
    boost::add_edge(1, 3, g);
    boost::add_edge(1, 4, g);
    boost::add_edge(1, 5, g);

    std::map<std::pair<int, int>, int> weights = {
        {{0, 1}, 0},
        {{0, 2}, 0},
        {{0, 3}, 0},
        {{0, 4}, 0},
        {{0, 5}, 0},
        {{1, 2}, 1},
        {{1, 3}, 1},
        {{1, 4}, -3},
        {{1, 5}, -3}
    };

    std::map<std::pair<int, int>, int> noise = {
        {{0, 1}, 0},
        {{0, 2}, 0},
        {{0, 3}, 0},
        {{0, 4}, 0},
        {{0, 5}, 0},
        {{1, 2}, 0},
        {{1, 3}, 0},
        {{1, 4}, 0},
        {{1, 5}, 0},
    };

    assign_weights_and_noise(g, weights, noise);

    return g;
}

#endif //NEGATIVE_TRIANGLE_COUNTING_GRAPH_PRESETS_H
