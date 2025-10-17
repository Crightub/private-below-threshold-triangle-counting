#include "global_sensitivity.h"

#include <boost/graph/subgraph.hpp>

#include "counting.h"
#include "distribution.h"
#include "triangle.h"
#include "utils.h"
struct PrivateCountingConfig;

int global_sensitivity(const Graph &g, int i, const std::list<Triangle> &triangles) {
    int sens = 0;

    std::vector<Edge> incident_edges;
    for (auto [ei, ei_end] = boost::out_edges(i, g); ei != ei_end; ++ei) {
        incident_edges.push_back(*ei);
    }

    #pragma omp parallel for reduction(max:sens)
    for (size_t i = 0; i < incident_edges.size(); ++i) {
        Edge e = incident_edges[i];
        int count = 0;

        for (const auto &t : triangles) {
            for (Edge edge_in_triangle : t.edges) {
                if (is_same_edge(g, edge_in_triangle, e)) {
                    ++count;
                }
            }
        }

        sens = std::max(sens, count);
    }

    return sens;
}


void apply_global_sensitivity(const Graph &g,
                              const PrivateCountingConfig &cfg,
                              std::vector<TriangleCount> &counts,
                              const std::vector<std::list<Triangle> > &node2TriangleMap) {
    const double unbiased_sens_mul = 1 + 2 * (std::exp(cfg.weight_eps) / std::pow(1 - std::exp(cfg.weight_eps), 2));
    // BOOST_LOG_TRIVIAL(info) << "Global sensitivity multiplier: " << sens_mult << " for use_unbiased_estimator: " << cfg.use_unbiased_estimator;

    for (int i = 0; i < counts.size(); ++i) {
        std::list<Triangle> triangles = node2TriangleMap[i];
        double sens = global_sensitivity(g, i, triangles);
        double biased_global_noise = sample_laplace(0, sens / cfg.count_eps);
        double unbiased_global_noise = sample_laplace(0, sens * unbiased_sens_mul / cfg.count_eps);

        counts[i].biased += biased_global_noise;
        counts[i].unbiased += unbiased_global_noise;
    }
}
