#pragma once

#include "utils.h"

class Triangle {
public:
    std::vector<Node> nodes = std::vector<Node>(3);
    Node source_node; // source_node is responsible for counting the

    std::vector<Edge> edges = std::vector<Edge>(3);
    Edge noisy_edge; // estimator uses noisy edge

    bool is_assigned;

    Triangle(const Graph &g, const Node v, const Node u, const Node w) : source_node(-1) {
        nodes[0] = v;
        nodes[1] = u;
        nodes[2] = w;

        auto [e_vu, e_vu_exists] = boost::edge(v, u, g);
        edges[0] = e_vu;
        auto [e_vw, e_vw_exists] = boost::edge(v, w, g);
        edges[1] = e_vw;
        auto [e_uw, e_uw_exists] = boost::edge(u, w, g);
        edges[2] = e_uw;

        is_assigned = false;
    }

    Triangle(const Graph &g, const Edge &e1, const Edge &e2, const Edge &e3) : source_node(-1) {
        edges[0] = e1;
        edges[1] = e2;
        edges[2] = e3;

        Node u = boost::source(e1, g);
        Node w = boost::target(e1, g);

        Node v = boost::source(e1, g) == u || boost::source(e1, g) == w ? boost::target(e1, g) : boost::source(e1, g);

        nodes[0] = u;
        nodes[1] = w;
        nodes[2] = v;

        is_assigned = false;
    }

    void assign_triangle(Graph &g, bool use_load_balancing) {
        if (is_assigned) return;

        if (!use_load_balancing) {
            assign_by_index(g);
            return;
        }

        bool was_assigned = assign_by_triangle_count(g);

        if (was_assigned) {
            return;
        }

        assign_by_load(g);
    }

    std::tuple<int, int, int, int>
    get_triangle_weights(Graph &g) {
        auto w_e1 = get_edge_weight(g, noisy_edge);
        auto noise_e1 = get_noise(g);

        int w_e2 = 0;
        int w_e3 = 0;

        for (const auto &e: edges) {
            if (!is_same_edge(g, e, noisy_edge)) {
                if (w_e2 == 0)
                    w_e2 = get_edge_weight(g,e);
                else
                    w_e3 = get_edge_weight(g, e);
            }
        }

        return {noise_e1, w_e1, w_e2, w_e3};
    }


    friend std::ostream &operator<<(std::ostream &os, const Triangle &t);

    bool operator==(const Triangle &other) const {
        return nodes == other.nodes;
    }

    static int get_edge_weight(Graph &g, const Edge &e) {
        auto weight_map = get(&edge_info::weight, g);

        auto [x, exists] = boost::edge(boost::source(e, g), boost::target(e, g), g);
        if (!exists) {
            std::cerr << "ERROR: edge descriptor invalid!" << std::endl;
            return 0;
        }

        return get(weight_map, x);
    }

    int get_noise(Graph &g) {
        auto noise_map = get(&edge_info::noise, g);

        auto [x, exists] = boost::edge(boost::source(noisy_edge, g), boost::target(noisy_edge, g), g);
        if (!exists) {
            std::cerr << "ERROR: edge descriptor invalid!" << std::endl;
            return 0;
        }

        return get(noise_map, x);
    }

    bool contains_edge(Graph &g, const Edge &e) const {
        for (Edge triangle_e : edges) {
            if (is_same_edge(g, triangle_e, e))
                return true;
        }
        return false;
    }

private:
    void assign_by_index(Graph &g) {
        source_node = *std::min_element(nodes.begin(), nodes.end());
        update_noisy_edge(g);
        increase_edge_load(g, noisy_edge);
        is_assigned = true;
    }

    bool assign_by_triangle_count(Graph &g) {
        auto triangle_count_map = get(&edge_info::involved_triangles, g);

        for (Edge &e: edges) {
            int involved_triangle_count = get(triangle_count_map, e);
            if (involved_triangle_count == 1) {
                noisy_edge = e;
                update_noisy_edge(g);
                update_source_node(g);
                increase_edge_load(g, noisy_edge);
                is_assigned = true;
                return true;
            }
        }

        return false;
    }

    void assign_by_load(Graph &g) {
        int l0 = get_edge_load(g, edges[0]);
        int l1 = get_edge_load(g, edges[1]);
        int l2 = get_edge_load(g, edges[2]);

        int min_load = std::min({l0, l1, l2});

        std::vector<Edge> candidates;
        for (auto &p: edges) {
            if (get_edge_load(g, p) == min_load)
                candidates.push_back(p);
        }

        static std::random_device rd;
        static std::mt19937 gen(rd());
        std::uniform_int_distribution<> dist(0, (int) candidates.size() - 1);
        Edge chosen = candidates[dist(gen)];

        auto [ne, _] = boost::edge(boost::source(chosen, g), boost::target(chosen, g), g);
        noisy_edge = ne;

        increase_edge_load(g, noisy_edge);
        update_source_node(g);

        is_assigned = true;
    }

    static void increase_edge_load(Graph &g, Edge &e) {
        auto load_map = get(&edge_info::load, g);

        auto [x, exists] = boost::edge(boost::source(e, g), boost::target(e, g), g);
        if (!exists) {
            std::cerr << "ERROR: edge descriptor invalid!" << std::endl;
            return;
        }

        int load = get(load_map, x);
        put(&edge_info::load, g, x, load + 1);
    }

    static int get_edge_load(Graph &g, Edge &e) {
        auto load_map = get(&edge_info::load, g);

        auto [x, exists] = boost::edge(boost::source(e, g), boost::target(e, g), g);
        if (!exists) {
            std::cerr << "ERROR: edge descriptor invalid!" << std::endl;
            return 0;
        }

        return get(load_map, x);
    }


    void update_noisy_edge(const Graph &g) {
        for (auto &e: edges) {
            Node u = boost::source(e, g);
            Node v = boost::target(e, g);
            if (u != source_node && v != source_node) {
                auto [ne, _] = boost::edge(u, v, g);
                noisy_edge = ne;
                break;
            }
        }
    }

    void update_source_node(const Graph &g) {
        for (Node v: nodes) {
            if (v != boost::source(noisy_edge, g) && v != boost::target(noisy_edge, g)) {
                source_node = v;
            }
        }
    }
};