#include "graph_statistics.h"

std::vector<std::set<Node>> forward_adjacency_lists(const Graph &g, const int n, const std::vector<int> &rank) {
    auto forward = std::vector<std::set<Node>>(n, std::set<Node>());
    for (auto [vi, vi_end] = vertices(g); vi != vi_end; ++vi) {
        const Node v = *vi;
        for (auto [ai, ai_end] = adjacent_vertices(v, g); ai != ai_end; ++ai) {
            if (const Node u = *ai; rank[v] < rank[u]) {
                forward[v].insert(u);
            }
        }
    }

    return forward;
}

void compute_rank(const std::vector<size_t> &degree, const size_t n, std::vector<int> &order, std::vector<int> &rank) {
    for (int i = 0; i < n; ++i) order[i] = i;

    std::sort(order.begin(), order.end(),
              [&](int a, int b) {
                  if (degree[a] != degree[b]) return degree[a] > degree[b];
                  return a < b;
              });

    for (int i = 0; i < n; ++i)
        rank[order[i]] = i;

}

void degree_vector(const Graph &g, std::vector<size_t> &degree) {
    for (auto [vi, vi_end] = vertices(g); vi != vi_end; ++vi) {
        degree[*vi] = boost::degree(*vi, g);
    }
}

std::list<Triangle> find_triangles(Graph &g) {
    std::list<Triangle> triangles;
    auto involved_triangles_map = get(&edge_info::involved_triangles, g);

    const size_t n = boost::num_vertices(g);

    std::vector<size_t> degree(n);
    degree_vector(g, degree);

    // order: sorts the vertex indices by degree in decreasing order
    // rank: ranking of vertex i in according to degree
    std::vector<int> order(n), rank(n);
    compute_rank(degree, n, order, rank);

    auto forward = forward_adjacency_lists(g, n, rank);

    for (int i = 0; i < n; ++i) {
        const Node v = order[i];

        for (const Node u: forward[v]) {
            const auto &Nu = forward[u];
            const auto &Nv = forward[v];

            auto nu_itr = Nu.begin();
            auto nv_itr = Nv.begin();

            while (nu_itr != Nu.end() && nv_itr != Nv.end()) {
                if (*nu_itr == *nv_itr) {
                    Node w = *nu_itr;

                    triangles.emplace_back(g, v, u, w);

                    auto [e1, e1_exists] = edge(u, v, g);
                    auto [e2, e2_exists] = edge(u, w, g);
                    auto [e3, e3_exists] = edge(v, w, g);

                    involved_triangles_map[e1]++;
                    involved_triangles_map[e2]++;
                    involved_triangles_map[e3]++;

                    ++nu_itr;
                    ++nv_itr;
                } else if (*nu_itr < *nv_itr) {
                    ++nu_itr;
                } else if (*nv_itr < *nu_itr) {
                    ++nv_itr;
                }
            }
        }
    }

    return triangles;
}


std::list<Triangle> find_triangles_naive(const Graph &g) {
    std::list<Triangle> triangles;

    size_t n = boost::num_vertices(g);

    for (int v = 0; v < n; ++v) {
        for (int u = v + 1; u < n; ++u) {
            auto [e_vu, e_vu_exists] = boost::edge(u, v, g);
            if (!e_vu_exists) continue;
            for (int w = u + 1; w < n; ++w) {
                auto [e_wu, e_wu_exists] = boost::edge(w, u, g);
                if (!e_wu_exists) continue;

                auto [e_vw, e_vw_exists] = boost::edge(w, v, g);
                if (!e_vw_exists) continue;

                if (v >= u || v >= w || u >= w) {
                    std::cout << "Error" << std::endl;
                }

                triangles.push_back(Triangle(g, u, v, w));
            }
        }
    }

    return triangles;
}
