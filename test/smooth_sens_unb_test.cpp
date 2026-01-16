#define BOOST_TEST_MODULE Smooth_Sens_Unb_Test

#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/context.hpp>

#include "graph_presets.h"
#include "graph_statistics.h"
#include "smooth_sens_unbiased.h"
#include "triangle.h"

std::pair<std::vector<Triangle>, std::vector<std::list<int> > > setup_node_triangle_map(Graph &g, const std::size_t n) {
    std::vector<Triangle> triangles = find_triangles(g);

    std::vector<std::list<int> > node_triangle_map = std::vector<std::list<int> >(n);

    for (int i = 0; i < triangles.size(); i++) {
        auto& t = triangles[i];
        t.assign_triangle(g, false);
        node_triangle_map[t.source_node].push_back(i);
    }

    return {triangles, node_triangle_map};
}


BOOST_AUTO_TEST_CASE(smooth_basic) {
    Graph g = example_graph_7();
    const size_t n = boost::num_vertices(g);
    auto [triangles, assignment] = setup_node_triangle_map(g, n);
    auto [e, e_vu_exists] = boost::edge(0, 1, g);

    double beta = 1;
    double eps = 1;
    double p = std::exp(-eps);
    double x = p/std::pow(1-p, 2);

    double sens_pos = fix_edge_sens_pos(g, e, 0, assignment[0], triangles, beta, p, true);
    BOOST_CHECK_CLOSE(sens_pos, x*std::exp(0), 1e-9);

    double sens_neg = fix_edge_sens_neg(g, e, 0, assignment[0], triangles, beta, p, true);
    BOOST_CHECK_CLOSE(sens_neg, (1+2*x)*std::exp(-1), 1e-9);

    double sens_pos_dec = fix_edge_sens_neg(g, e, 0, assignment[0], triangles, beta, p, false);
    BOOST_CHECK_CLOSE(sens_pos_dec, (1+2*x), 1e-9);

    double sens_neg_dec = fix_edge_sens_pos(g, e, 0, assignment[0], triangles, beta, p, false);
    BOOST_CHECK_CLOSE(sens_neg_dec, x*std::exp(-1), 1e-9);
}

BOOST_AUTO_TEST_CASE(smooth_basic_2) {
    Graph g = example_graph_8();
    const size_t n = boost::num_vertices(g);
    auto [triangles, assignment] = setup_node_triangle_map(g, n);
    auto [e, e_vu_exists] = boost::edge(0, 1, g);

    double beta = 1;
    double eps = 1;
    double p = std::exp(-eps);
    double x = p/std::pow(1-p, 2);

    double sens_pos = fix_edge_sens_pos(g, e, 0, assignment[0], triangles, beta, p, true);
    BOOST_CHECK_CLOSE(sens_pos, 2*x*std::exp(-1), 1e-9);

    double sens_neg = fix_edge_sens_neg(g, e, 0, assignment[0], triangles, beta, p, true);
    BOOST_CHECK_CLOSE(sens_neg, (1+2*x), 1e-9);

    double sens_pos_dec = fix_edge_sens_neg(g, e, 0, assignment[0], triangles, beta, p, false);
    BOOST_CHECK_CLOSE(sens_pos_dec, 1*(1+2*x)*std::exp(-1), 1e-9);

    double sens_neg_dec = fix_edge_sens_pos(g, e, 0, assignment[0], triangles, beta, p, false);
    BOOST_CHECK_CLOSE(sens_neg_dec, 2*x, 1e-9);
}

BOOST_AUTO_TEST_CASE(smooth_basic_3) {
    Graph g = example_graph_9();

    const size_t n = boost::num_vertices(g);
    auto [triangles, assignment] = setup_node_triangle_map(g, n);
    auto [e, e_vu_exists] = boost::edge(0, 1, g);

    double beta = 1;
    double eps = 1;
    double p = std::exp(-eps);
    double x = p/std::pow(1-p, 2);

    double sens_pos = fix_edge_sens_pos(g, e, 0, assignment[0], triangles, beta, p, true);
    BOOST_CHECK_CLOSE(sens_pos, x*std::exp(-2), 1e-9);

    double sens_neg = fix_edge_sens_neg(g, e, 0, assignment[0], triangles, beta, p, true);
    BOOST_CHECK_CLOSE(sens_neg, (1+2*x)*std::exp(-3), 1e-9);

    double sens_pos_dec = fix_edge_sens_neg(g, e, 0, assignment[0], triangles, beta, p, false);
    BOOST_CHECK_CLOSE(sens_pos_dec, (1+2*x)*std::exp(-2), 1e-9);

    double sens_neg_dec = fix_edge_sens_pos(g, e, 0, assignment[0], triangles, beta, p, false);
    BOOST_CHECK_CLOSE(sens_neg_dec, x*std::exp(-1), 1e-9);

    beta = 0.1;
    sens_pos = fix_edge_sens_pos(g, e, 0, assignment[0], triangles, beta, p, true);
    BOOST_CHECK_CLOSE(sens_pos, 3*x*std::exp(-5*beta), 1e-9);

    sens_neg = fix_edge_sens_neg(g, e, 0, assignment[0], triangles, beta, p, true);
    BOOST_CHECK_CLOSE(sens_neg, 3*(1+2*x)*std::exp(-6*beta), 1e-9);

    sens_pos_dec = fix_edge_sens_neg(g, e, 0, assignment[0], triangles, beta, p, false);
    BOOST_CHECK_CLOSE(sens_pos_dec, 3*(1+2*x)*std::exp(-7*beta), 1e-9);

    sens_neg_dec = fix_edge_sens_pos(g, e, 0, assignment[0], triangles, beta, p, false);
    BOOST_CHECK_CLOSE(sens_neg_dec, 3*x*std::exp(-6*beta), 1e-9);
}
