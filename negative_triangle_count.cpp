#include <iostream>
#include <boost/graph/graphviz.hpp>
#include <boost/property_map/property_map.hpp>
#include "utils.h"
#include "distribution.h"
#include "graph_generation.h"
#include "qp_counting.h"
#include "basic_counting.h"


class edge_label_writer
{
public:
    edge_label_writer(const Graph &g) : g_(g) {}

    template <class Edge>
    void operator()(std::ostream &out, const Edge &e) const
    {
        const auto &ei = g_[e];
        out << "[label=\"w=" << ei.weight << ", w' = " << ei.weight + ei.noise << "\"]";
    }

private:
    const Graph &g_;
};

void save_graph(Graph g, std::string filename)
{
    std::ofstream file(filename);
    boost::write_graphviz(file, g,
                          boost::default_writer(),
                          edge_label_writer(g));
    file.close();
}

int main(int argc, char *argv[])
{
    if (argc < 7)
    {
        std::cerr << "Usage: " << argv[0] << " <num_vertices> <edge_probability> <weight_mu> <weight_std> <eps> <eps2> [-s]" << std::endl;
        return 1;
    }

    int num_vertices = std::atoi(argv[1]);
    double edge_probability = std::atof(argv[2]);
    double weight_mu = std::atof(argv[3]);
    double weight_std = std::atof(argv[4]);
    double eps = std::atof(argv[5]);
    double eps2 = std::atof(argv[6]);

    bool save_graph_flag = false;
    if (argc >= 8 && std::string(argv[7]) == "-s")
    {
        save_graph_flag = true;
    }

    Graph g = generate_graph(num_vertices, edge_probability, weight_mu, weight_std, eps);

    std::cout << "Generated graph with " << boost::num_vertices(g) << " vertices and "
              << boost::num_edges(g) << " edges." << std::endl;

    int negative_triangles = count_negative_triangles(g);
    std::cout << "opt: " << negative_triangles << std::endl;

    // double smallest_index_count = randomized_private_counting(g, eps, eps2, false, false, false);
    // std::cout << "smallest index: " << smallest_index_count << std::endl;
    //
    // double optimized_count = QPCountNegativeTriangles(g, eps, eps2);
    // std::cout << "optimized count: " << optimized_count << std::endl;

    if (save_graph_flag)
    {
        save_graph(g, "graph.dot");
    }
    return 0;
}
