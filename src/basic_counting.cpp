#include "graph_generation.hpp"
#include "basic_counting.hpp"
#include "distribution.hpp"

int count_negative_triangles(const Graph &g)
{
    using Vertex = typename boost::graph_traits<Graph>::vertex_descriptor;
    using Edge = typename boost::graph_traits<Graph>::edge_descriptor;

    std::size_t n = num_vertices(g);

    std::vector<std::vector<bool>> A(n, std::vector<bool>(n, false));
    std::vector<std::vector<double>> W(n, std::vector<double>(n, 0.0));

    auto index_map = get(boost::vertex_index, g);
    auto weight_map = get(&edge_info::weight, g);

    for (auto e : boost::make_iterator_range(edges(g)))
    {
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
    for (int i = 0; i < n; ++i)
    {
        for (int j = i + 1; j < n; ++j)
        {
            if (!A[i][j])
                continue;
            for (int k = j + 1; k < n; ++k)
            {
                if (A[j][k] && A[k][i])
                {
                    double sum = W[i][j] + W[j][k] + W[k][i];
                    if (sum < 0)
                        ++count;
                }
            }
        }
    }

    return count;
}

PrivateCountingResult randomized_private_counting(const Graph &g, double eps, double eps2, bool use_weight_noise, bool use_count_noise, bool use_biased_estimator)
{
    using Vertex = typename boost::graph_traits<Graph>::vertex_descriptor;
    using Edge = typename boost::graph_traits<Graph>::edge_descriptor;

    double lambda = std::exp(-eps / 2.0);

    std::size_t n = num_vertices(g);

    std::vector<std::vector<bool>> A(n, std::vector<bool>(n, false));
    std::vector<std::vector<double>> W(n, std::vector<double>(n, 0.0));
    std::vector<std::vector<double>> N(n, std::vector<double>(n, 0.0));

    auto index_map = get(boost::vertex_index, g);
    auto weight_map = get(&edge_info::weight, g);
    auto noise_map = get(&edge_info::noise, g);

    for (auto e : boost::make_iterator_range(edges(g)))
    {
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

    double count = 0;
    double count_noise = 0;
    double weight_noise = 0;

#pragma omp parallel for reduction(+ : count) schedule(dynamic)
    for (std::size_t i = 0; i < n; ++i)
    {
        double count_i = 0;
        double triangle_i = 0;
        for (std::size_t j = i + 1; j < n; ++j)
        {
            if (!A[i][j])
                continue;
            for (std::size_t k = j + 1; k < n; ++k)
            {
                if (A[j][k] && A[k][i])
                {
                    triangle_i += 1;
                    double unknown_weight;

                    if (use_weight_noise){
                        unknown_weight = N[j][k];
                    } else {
                        unknown_weight = W[j][k];
                    }

                    if(use_biased_estimator){
                        count_i += biased_estimator(W[i][j], W[k][i], unknown_weight);
                        weight_noise += std::abs(biased_estimator(W[i][j], W[k][i], N[j][k]) - biased_estimator(W[i][j], W[k][i], W[j][k]));
                    } else {
                        count_i += unbiased_estimator(-W[i][j] - W[k][i] - unknown_weight, lambda);
                        weight_noise += unbiased_estimator(-W[i][j] - W[k][i] - N[j][k], lambda) - biased_estimator(W[i][j], W[k][i], W[j][k]);
                    }
                }
            }
        }

        if (!use_count_noise){
            count += count_i;
            continue;
        }

        if (triangle_i > 0)
        {
            double sensitivity;
            if (use_biased_estimator){
                sensitivity = triangle_i;
            } else {
                sensitivity = triangle_i*(1 + 2*lambda / std::pow((1 - lambda), 2));
            }

            double local_count_noise = sampleLaplace(0, sensitivity / eps2);
            // std::cout << "local_count_noise: " << local_count_noise << std::endl;
            count += count_i + local_count_noise;
            count_noise += local_count_noise;
        }
    }

    return PrivateCountingResult{count, weight_noise, count_noise};
}