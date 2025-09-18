#include "graph_generation.h"
#include "basic_counting.h"
#include "distribution.h"
#include "smooth_sensitivity.h"
#include "global_sensitivity.h"

int count_negative_triangles(const Graph &g) {
    using Vertex = typename boost::graph_traits<Graph>::vertex_descriptor;
    using Edge = typename boost::graph_traits<Graph>::edge_descriptor;

    std::size_t n = num_vertices(g);

    std::vector<std::vector<bool> > A(n, std::vector<bool>(n, false));
    std::vector<std::vector<double> > W(n, std::vector<double>(n, 0.0));

    auto index_map = get(boost::vertex_index, g);
    auto weight_map = get(&edge_info::weight, g);

    for (auto e: boost::make_iterator_range(edges(g))) {
        Vertex u = source(e, g);
        Vertex v = target(e, g);

        int ui = get(index_map, u);
        int vi = get(index_map, v);

        A[ui][vi] = true;
        A[vi][ui] = true;

        double w = get(weight_map, e);
        W[ui][vi] = w;
        W[vi][ui] = w;
    }

    int count = 0;

#pragma omp parallel for reduction(+ : count) schedule(dynamic)
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (!A[i][j])
                continue;
            for (int k = j + 1; k < n; ++k) {
                if (A[j][k] && A[k][i]) {
                    double sum = W[i][j] + W[j][k] + W[k][i];
                    if (sum < 0)
                        ++count;
                }
            }
        }
    }

    return count;
}

void setup_graph_matrices(const Graph &g,
                          std::vector<std::vector<bool> > &A,
                          std::vector<std::vector<int> > &W,
                          std::vector<std::vector<int> > &N) {
    using Vertex = typename boost::graph_traits<Graph>::vertex_descriptor;

    auto index_map = get(boost::vertex_index, g);
    auto weight_map = get(&edge_info::weight, g);
    auto noise_map = get(&edge_info::noise, g);

    for (auto e: boost::make_iterator_range(edges(g))) {
        Vertex u = source(e, g);
        Vertex v = target(e, g);

        int ui = get(index_map, u);
        int vi = get(index_map, v);

        A[ui][vi] = true;
        A[vi][ui] = true;

        double w = get(weight_map, e);
        double noise = get(noise_map, e);
        W[ui][vi] = w;
        W[vi][ui] = w;
        N[ui][vi] = w + noise;
        N[vi][ui] = w + noise;
    }
}

double count_triangle(const PrivateCountingConfig &cfg, int w_ij, int w_ik, int n_jk) {
    // std::cout << "count_triangle: w_ij = " << w_ij << ", w_ik = " << w_ik << ", n_jk = " << n_jk << std::endl;
    if (cfg.use_unbiased_estimator) {
        return unbiased_estimator(w_ij + w_ik + n_jk, std::exp(-cfg.weight_eps));
    } else {
        return biased_estimator(w_ij, w_ik, n_jk);
    }
}


void count_local_negative_triangles(const PrivateCountingConfig &cfg,
                                    std::vector<std::vector<bool> > &A,
                                    std::vector<std::vector<int> > &W,
                                    std::vector<std::vector<int> > &N,
                                    std::vector<double> &counts
) {
    int n = counts.size();
    for (std::size_t i = 0; i < n; ++i) {
        double count_i = 0;
#pragma omp parallel for reduction(+ : count_i) schedule(dynamic)
        for (std::size_t j = i + 1; j < n; ++j) {
            if (!A[i][j])
                continue;

            for (std::size_t k = j + 1; k < n; ++k) {
                if (A[j][k] && A[k][i]) {
                    // Triangle counted by i
                    count_i += count_triangle(cfg, W[i][j], W[i][k], N[j][k]);
                }
            }
        }
        counts[i] = count_i;
    }
}

void publish_local_counts(const Graph &g,
                          const std::vector<std::vector<bool> > &A,
                          const PrivateCountingConfig &cfg,
                          std::vector<double> &counts) {
    if (!cfg.use_smooth_sensitivity) {
        apply_global_sensitivity(A, cfg, counts);
    }

    if (cfg.use_smooth_sensitivity) {
        apply_smooth_sensitivity(g, cfg, counts);
    }
}

PrivateCountingResult randomized_private_counting(const Graph &g, const PrivateCountingConfig &cfg) {
    std::size_t n = num_vertices(g);

    std::vector<std::vector<bool> > A(n, std::vector<bool>(n, false));
    std::vector<std::vector<int> > W(n, std::vector<int>(n, 0));
    std::vector<std::vector<int> > N(n, std::vector<int>(n, 0));

    setup_graph_matrices(g, A, W, N);

    std::vector<double> counts(n, 0.0);

    count_local_negative_triangles(cfg, A, W, N, counts);
    publish_local_counts(g, A, cfg, counts);

    double total_count = 0;
    for (double c: counts) {
        total_count += c;
    }

    return PrivateCountingResult{total_count, total_count, total_count};
}
