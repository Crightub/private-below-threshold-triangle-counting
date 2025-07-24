#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/math/distributions/laplace.hpp>
#include <tuple>
#include <set>
#include <omp.h>
#include <random>

struct edge_info{
    int weight;
    double noise;
};
 
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, edge_info> Graph;
typedef boost::sorted_erdos_renyi_iterator<boost::minstd_rand, Graph> ERGen;
typedef boost::minstd_rand RandomType;

Graph generateERGraph(size_t n, double p){
    RandomType gen;
    Graph g(ERGen(gen, n, p), ERGen(), n);
    return g;
}

Graph assignGaussianWeights(Graph g, double mu, double std){
    boost::normal_distribution<> nd(mu, std);
    RandomType gen;

    boost::variate_generator<RandomType&, boost::normal_distribution<> > var_n(gen, nd);

    boost::graph_traits<Graph>::edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei) {
        int weight = var_n();
        std::cout << "weight: " << weight << std::endl;
        put(&edge_info::weight, g, *ei, weight);
        std::cout << "stored weight: " << get(&edge_info::weight, g)[*ei] << std::endl;
    }

    return g;
}

int sampleDiscreteLaplace(double lambda, RandomType &rng) {
    std::geometric_distribution<int> geom(1 - lambda);
    
    int u_sample = geom(rng);
    int v_sample = geom(rng);

    return u_sample - v_sample;
}

double sampleLaplace(double mu, double b){
    RandomType gen;

    boost::random::uniform_real_distribution<> uniform(0.0, 1.0);
    double u = uniform(gen);

    boost::math::laplace_distribution<> laplace(mu, b);
    return boost::math::quantile(laplace, u);
}

Graph addDiscreteLaplaceNoise(Graph g, double eps) {
    RandomType gen;
    double lambda = std::exp(-eps / 2.0);

    boost::graph_traits<Graph>::edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei) {
        int noise = sampleDiscreteLaplace(lambda, gen);
        put(&edge_info::noise, g, *ei, noise);
        std::cout << "stored noise: " << get(&edge_info::noise, g)[*ei] << std::endl;
    }

    return g;
}


template <typename Graph>
std::size_t count_negative_triangles(const Graph& g) {
    using Vertex = typename boost::graph_traits<Graph>::vertex_descriptor;
    using Edge = typename boost::graph_traits<Graph>::edge_descriptor;

    std::size_t n = num_vertices(g);

    std::vector<std::vector<bool>> A(n, std::vector<bool>(n, false));
    std::vector<std::vector<double>> W(n, std::vector<double>(n, 0.0));

    auto index_map = get(boost::vertex_index, g);
    auto weight_map = get(&edge_info::weight, g);

    for (auto e : boost::make_iterator_range(edges(g))) {
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

    std::size_t count = 0;

#pragma omp parallel for reduction(+:count) schedule(dynamic)
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            if (!A[i][j]) continue;
            for (std::size_t k = j + 1; k < n; ++k) {
                if (A[j][k] && A[k][i]) {
                    double sum = W[i][j] + W[j][k] + W[k][i];
                    if (sum < 0) ++count;
                }
            }
        }
    }

    return count;
}

double unbiased_estimator(int x, double p){
    if (x > 0){
        return 1.0;
    }
    if (x == 0){
        1 + p/std::pow((1-p),2);
    }
    if(x == - 1){
        -p/std::pow((1-p),2);
    }
    if (x < 0){
        return 0;
    }

    return 0;
}

double randomized_private_counting(const Graph &g, double eps, double eps2){
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

    for (auto e : boost::make_iterator_range(edges(g))) {
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

#pragma omp parallel for reduction(+:count) schedule(dynamic)
    for (std::size_t i = 0; i < n; ++i) {
        double count_i = 0;
        double triangle_i = 0;
        for (std::size_t j = i + 1; j < n; ++j) {
            if (!A[i][j]) continue;
            for (std::size_t k = j + 1; k < n; ++k) {
                if (A[j][k] && A[k][i]) {
                    triangle_i += 1;
                    // Let i count the triangle -> uses noisy weight for edge (j,k)
                    count_i += unbiased_estimator(-W[i][j] -  W[k][i] - N[j][k], lambda);
                }
            }
        }
        
        if (triangle_i > 0){
            double count_noise = sampleLaplace(0, triangle_i / eps2);
            std::cout << "Node " << i << ": triangles: " << triangle_i << ", negative: " << count_i << ", count noise: " << count_noise << std::endl;
            count += count_i + count_noise;
        }
    }

    return count;
}


class edge_label_writer {
public:
    edge_label_writer(const Graph& g) : g_(g) {}

    template <class Edge>
    void operator()(std::ostream& out, const Edge& e) const {
        const auto& ei = g_[e];
        out << "[label=\"w=" << ei.weight + ei.noise << "\"]";
    }

private:
    const Graph& g_;
};


void save_graph(Graph g, std::string filename){
    std::ofstream file(filename);
    boost::write_graphviz(file, g,
        boost::default_writer(),
        edge_label_writer(g));
    file.close();
}

int main() {
    size_t num_vertices = 10;
    double edge_probability = 0.4; 
    double eps = 1;
    double eps2 = 10;

    Graph g = generateERGraph(num_vertices, edge_probability);
    g = assignGaussianWeights(g, 10, 20);
    Graph noisyGraph = addDiscreteLaplaceNoise(g, eps);

    std::cout << "Generated graph with " << boost::num_vertices(g) << " vertices and "
              << boost::num_edges(g) << " edges." << std::endl;

    std::size_t negative_triangles = count_negative_triangles(noisyGraph);

    double smallest_index_count = randomized_private_counting(noisyGraph, eps, eps2);

    std::cout << "opt: " << negative_triangles << std::endl; 
    std::cout << "smallest index: " << smallest_index_count << std::endl;

    save_graph(noisyGraph, "graph.dot");
    return 0;
}

