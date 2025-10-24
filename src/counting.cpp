#include "graph_generation.h"
#include "counting.h"

#include "distribution.h"
#include "smooth_sensitivity.h"
#include "global_sensitivity.h"
#include "graph_statistics.h"
#include "triangle.h"


void count_local_negative_triangles(Graph &g,
                                    const PrivateCountingConfig &base_cfg,
                                    std::vector<TriangleCount> &counts,
                                    std::vector<std::list<Triangle>> &node_triangle_map,
                                    std::vector<Triangle> &triangles
) {
    const double p = std::exp(-base_cfg.weight_eps);

    int n_nodes = counts.size();

    // Thread-local storage
    int n_threads = omp_get_max_threads();
    std::vector<std::vector<TriangleCount>> local_counts(n_threads, std::vector<TriangleCount>(n_nodes));
    std::vector<std::vector<std::list<Triangle>>> local_maps(n_threads, std::vector<std::list<Triangle>>(n_nodes));

    #pragma omp parallel for
    for (size_t i = 0; i < triangles.size(); i++) {
        int tid = omp_get_thread_num();
        Triangle &t = triangles[i];

        if (!t.is_assigned) {
            #pragma omp critical
            std::cerr << "Triangle is not yet assigned!" << std::endl;
        }

        auto [noise_e1, w_e1, w_e2, w_e3] = t.get_triangle_weights(g);
        auto [noise_1, noise_2, noise_3] = t.get_triangle_noise(g);

        int node = t.source_node;
        local_counts[tid][node].naive += biased_estimator(w_e1 + noise_1, w_e2 + noise_2, w_e3 + noise_3, base_cfg.lambda);
        local_counts[tid][node].opt +=  biased_estimator(w_e1, w_e2, w_e3, base_cfg.lambda);
        local_counts[tid][node].unbiased += unbiased_estimator(noise_e1 + w_e1 + w_e2 + w_e3, p, base_cfg.lambda);
        local_counts[tid][node].biased += biased_estimator(noise_e1 + w_e1, w_e2, w_e3, base_cfg.lambda);

        local_maps[tid][node].push_back(t);
    }

    // Merge thread-local results
    for (int tid = 0; tid < n_threads; tid++) {
        for (int node = 0; node < n_nodes; node++) {
            counts[node].naive += local_counts[tid][node].naive;
            counts[node].opt += local_counts[tid][node].opt;
            counts[node].unbiased += local_counts[tid][node].unbiased;
            counts[node].biased += local_counts[tid][node].biased;

            node_triangle_map[node].splice(node_triangle_map[node].end(), local_maps[tid][node]);
        }
    }
}

void publish_local_counts(Graph &g,
                          const PrivateCountingConfig &cfg,
                          std::vector<TriangleCount> &counts,
                          std::vector<std::list<Triangle> > &node_triangle_map) {
    if (!cfg.use_smooth_sensitivity) {
        apply_global_sensitivity(g, cfg, counts, node_triangle_map);
    }

    if (cfg.use_smooth_sensitivity) {
        apply_smooth_sensitivity(g, cfg, counts, node_triangle_map);
    }
}

PrivateCountingResult private_counting(
    Graph &g,
    PrivateCountingConfig &base_cfg,
    const std::vector<Triangle> *ptr_triangles) {
    std::size_t n = num_vertices(g);
    std::vector<TriangleCount> global_counts(n);
    std::vector<std::list<Triangle> > node_triangle_map = std::vector<std::list<Triangle> >(n);
    std::vector<Triangle> triangles;

    if (ptr_triangles == nullptr) {
        triangles = find_triangles(g);
    } else {
        triangles = *ptr_triangles;
    }

    std::cout << "Count local below-threshold triangles." << std::endl;
    count_local_negative_triangles(g, base_cfg, global_counts, node_triangle_map, triangles);

    auto smooth_counts = global_counts;

    std::cout << "Publish counts with global sensitivity." << std::endl;
    auto global_cfg = base_cfg;
    global_cfg.use_smooth_sensitivity = false;
    publish_local_counts(g, global_cfg, global_counts, node_triangle_map);

    std::cout << "Publish counts with smooth sensitivity." << std::endl;
    auto smooth_cfg = base_cfg;
    smooth_cfg.use_smooth_sensitivity = true;
    publish_local_counts(g, smooth_cfg, smooth_counts, node_triangle_map);

    long opt = 0;
    long naive = 0;
    double global_unbiased = 0;
    double global_biased = 0;
    for (TriangleCount c: global_counts) {
        opt += c.opt;
        naive += c.naive;
        global_unbiased += c.unbiased;
        global_biased += c.biased;
    }

    double smooth_unbiased = 0;
    double smooth_biased = 0;
    for (TriangleCount c: smooth_counts) {
        smooth_unbiased += c.unbiased;
        smooth_biased += c.biased;
    }

    return PrivateCountingResult{opt, naive, global_unbiased, global_biased, smooth_unbiased, smooth_biased};
}
