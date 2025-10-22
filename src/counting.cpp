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
                                    std::vector<std::list<Triangle> > &node_triangle_map,
                                    std::vector<Triangle> &triangles
) {
    const double p = std::exp(-base_cfg.weight_eps);

    for (Triangle &t: triangles) {
        t.assign_triangle(g, base_cfg.use_load_balancing);
        auto [noise_e1, w_e1, w_e2, w_e3] = t.get_triangle_weights(g);

        counts[t.source_node].opt +=  biased_estimator(w_e1, w_e2, w_e3, base_cfg.lambda);
        counts[t.source_node].unbiased += unbiased_estimator(noise_e1 + w_e1 + w_e2 + w_e3, p, base_cfg.lambda);;
        counts[t.source_node].biased += biased_estimator(noise_e1 + w_e1, w_e2, w_e3, base_cfg.lambda);
        node_triangle_map[t.source_node].push_back(t);
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

void print_c4_info(const Graph &g, const std::vector<std::list<Triangle> > &node_triangle_map,
                   const bool use_load_balancing) {
    int max_size = 0;
    for (auto triangle_list: node_triangle_map) {
        max_size = std::max(max_size, static_cast<int>(triangle_list.size()));
    }
    std::cout << "Max triangle list sizes: " << max_size << std::endl;

    long long c4_instances = compute_c4_instances(g);
    std::cout << "#C4: " << c4_instances << " with use_load_balancing=" << use_load_balancing << std::endl;
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

    std::cout << "Count local negative triangles." << std::endl;
    count_local_negative_triangles(g, base_cfg, global_counts, node_triangle_map, triangles);

    // print_c4_info(g, node_triangle_map, use_load_balancing);

    auto smooth_counts = global_counts;

    std::cout << "Publish counts with global sensitivity." << std::endl;
    auto global_cfg = base_cfg;
    global_cfg.use_smooth_sensitivity = false;
    publish_local_counts(g, global_cfg, global_counts, node_triangle_map);

    // std::cout << "Publish counts with smooth sensitivity." << std::endl;
    // auto smooth_cfg = base_cfg;
    // smooth_cfg.use_smooth_sensitivity = true;
    // publish_local_counts(g, smooth_cfg, smooth_counts, node_triangle_map);

    int opt = 0;
    double global_unbiased = 0;
    double global_biased = 0;
    for (TriangleCount c: global_counts) {
        opt += c.opt;
        global_unbiased += c.unbiased;
        global_biased += c.biased;
    }

    double smooth_unbiased = 0;
    double smooth_biased = 0;
    for (TriangleCount c: smooth_counts) {
        smooth_unbiased += c.unbiased;
        smooth_biased += c.biased;
    }

    return PrivateCountingResult{opt, global_unbiased, global_biased, smooth_unbiased, smooth_biased};
}
