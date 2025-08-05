#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "distribution.hpp"
#include "utils.hpp"
#include "graph_generation.hpp"

Graph generateERGraph(size_t n, double p)
{
    Graph g(ERGen(rng, n, p), ERGen(), n);
    return g;
}

Graph assignGaussianWeights(Graph g, double mu, double std)
{
    boost::normal_distribution<> nd(mu, std);
    boost::variate_generator<boost::minstd_rand, boost::normal_distribution<>> var_n(rng, nd);
    boost::graph_traits<Graph>::edge_iterator ei, ei_end;

    for (boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei)
    {
        int weight = var_n();
        put(&edge_info::weight, g, *ei, weight);
    }

    return g;
}

Graph addDiscreteLaplaceNoise(Graph g, double eps)
{   
    double lambda = std::exp(-eps / 2.0);

    boost::graph_traits<Graph>::edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei)
    {
        int noise = sampleDiscreteLaplace(lambda);
        put(&edge_info::noise, g, *ei, noise);
    }

    return g;
}

Graph generateGraph(int num_vertices, double edge_probability, double weight_mu, double weight_std, double weight_eps){
    Graph g = generateERGraph(num_vertices, edge_probability);
    g = assignGaussianWeights(g, weight_mu, weight_std);
    Graph noisyGraph = addDiscreteLaplaceNoise(g, weight_eps);

    return noisyGraph;
}