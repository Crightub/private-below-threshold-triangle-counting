#include "graph_statistics.h"
#include <omp.h>

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

std::vector<Triangle> find_triangles(Graph &g) {
    auto involved_triangles_map = get(&edge_info::involved_triangles, g);

    const size_t n = boost::num_vertices(g);

    std::vector<size_t> degree(n);
    degree_vector(g, degree);

    std::vector<int> order(n), rank(n);
    compute_rank(degree, n, order, rank);

    auto forward = forward_adjacency_lists(g, n, rank);

    const int max_threads = 4;
    omp_set_num_threads(max_threads);
    std::vector<std::vector<Triangle>> local_triangles(max_threads);

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < static_cast<int>(n); ++i) {
        const Node v = order[i];
        const int tid = omp_get_thread_num();
        auto &local_vec = local_triangles[tid];

        for (const Node u : forward[v]) {
            const auto &Nu = forward[u];
            const auto &Nv = forward[v];

            auto nu_itr = Nu.begin();
            auto nv_itr = Nv.begin();

            while (nu_itr != Nu.end() && nv_itr != Nv.end()) {
                if (*nu_itr == *nv_itr) {
                    const Node w = *nu_itr;

                    local_vec.emplace_back(g, v, u, w);

                    auto [e1, e1_exists] = edge(u, v, g);
                    auto [e2, e2_exists] = edge(u, w, g);
                    auto [e3, e3_exists] = edge(v, w, g);

#pragma omp atomic
                    involved_triangles_map[e1]++;
#pragma omp atomic
                    involved_triangles_map[e2]++;
#pragma omp atomic
                    involved_triangles_map[e3]++;

                    ++nu_itr;
                    ++nv_itr;
                } else if (*nu_itr < *nv_itr) {
                    ++nu_itr;
                } else {
                    ++nv_itr;
                }
            }
        }
    }

    // === Merge all local triangle vectors into a single contiguous vector ===
    size_t total_size = 0;
    for (const auto &vec : local_triangles)
        total_size += vec.size();

    std::vector<Triangle> triangles;
    triangles.reserve(total_size);

    for (auto &vec : local_triangles) {
        triangles.insert(triangles.end(),
                         std::make_move_iterator(vec.begin()),
                         std::make_move_iterator(vec.end()));
        vec.clear();
    }

    return triangles;
}
