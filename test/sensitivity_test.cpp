#define BOOST_TEST_MODULE SensitivityTest
#include <boost/graph/graph_concepts.hpp>
#include <boost/test/included/unit_test.hpp>

#include "distribution.h"
#include "graph_presets.h"
#include "global_sensitivity.h"
#include "graph_statistics.h"
#include "smooth_sensitivity.h"

std::pair<std::vector<Triangle>, std::vector<std::list<int> > > setup_node_triangle_map(Graph &g, const std::size_t n) {
    std::vector<Triangle> triangles = find_triangles(g);

    std::vector<std::list<int> > node_triangle_map = std::vector<std::list<int> >(n);

    for (int i = 0; i < triangles.size(); i++) {
        auto &t = triangles[i];
        t.assign_triangle(g, false);
        node_triangle_map[t.source_node].push_back(i);
    }

    return {triangles, node_triangle_map};
}

BOOST_AUTO_TEST_CASE(sensitivity_1) {
    Graph g = example_graph_1();
    const size_t n = boost::num_vertices(g);
    auto [triangles, node_triangle_map] = setup_node_triangle_map(g, n);

    BOOST_CHECK_EQUAL(global_sensitivity(g, 0, node_triangle_map[0], triangles), 6);
    for (int i = 1; i < n; i++) {
        BOOST_CHECK_EQUAL(global_sensitivity(g, i, node_triangle_map[i], triangles), 0);
    }

    BOOST_CHECK_EQUAL(smooth_sensitivity(g, 0, 0, node_triangle_map[0], triangles, 1), 1);
    for (int i = 1; i < n; i++) {
        BOOST_CHECK_EQUAL(smooth_sensitivity(g, i, 0, node_triangle_map[i], triangles, 1), 0);
    }
}


BOOST_AUTO_TEST_CASE(sensitivity_2) {
    Graph g = example_graph_2();
    const size_t n = boost::num_vertices(g);
    auto [triangles, node_triangle_map] = setup_node_triangle_map(g, n);

    BOOST_CHECK_EQUAL(global_sensitivity(g, 0, node_triangle_map[0], triangles), 2);
    for (int i = 1; i < n; i++) {
        BOOST_CHECK_EQUAL(global_sensitivity(g, i, node_triangle_map[i], triangles), 0);
    }

    BOOST_CHECK_EQUAL(smooth_sensitivity(g, 0, 0, node_triangle_map[0], triangles, 1), 1);
    for (int i = 1; i < n; i++) {
        BOOST_CHECK_EQUAL(smooth_sensitivity(g, i, 0, node_triangle_map[i],triangles, 1), 0);
    }
}


BOOST_AUTO_TEST_CASE(sensitivity_2_scaled) {
    Graph g = example_graph_2_scaled();
    const size_t n = boost::num_vertices(g);
    auto [triangles, node_triangle_map] = setup_node_triangle_map(g, n);

    BOOST_CHECK_EQUAL(global_sensitivity(g, 0, node_triangle_map[0], triangles), 2);
    for (int i = 1; i < n; i++) {
        BOOST_CHECK_EQUAL(global_sensitivity(g, i, node_triangle_map[i], triangles), 0);
    }

    BOOST_CHECK_CLOSE(smooth_sensitivity(g, 0, 5, node_triangle_map[0], triangles, 1), 2*std::exp(-2), 1e-9);
    for (int i = 1; i < n; i++) {
        BOOST_CHECK_EQUAL(smooth_sensitivity(g, i, 5, node_triangle_map[i], triangles, 1), 0);
    }

    BOOST_CHECK_EQUAL(smooth_sensitivity(g, 0, 10, node_triangle_map[0], triangles, 1), 1);
    for (int i = 1; i < n; i++) {
        BOOST_CHECK_EQUAL(smooth_sensitivity(g, i, 10, node_triangle_map[i], triangles, 1), 0);
    }
}

BOOST_AUTO_TEST_CASE(sensitivity_3) {
    Graph g = example_graph_3();
    const size_t n = boost::num_vertices(g);
    auto [triangles, node_triangle_map] = setup_node_triangle_map(g, n);

    BOOST_CHECK_EQUAL(global_sensitivity(g, 0, node_triangle_map[0], triangles), 1);
    for (int i = 1; i < n; i++) {
        BOOST_CHECK_EQUAL(global_sensitivity(g, i, node_triangle_map[i], triangles), 0);
    }

    BOOST_TEST(smooth_sensitivity(g, 0, 0, node_triangle_map[0],triangles, 1) == 1, boost::test_tools::tolerance(1e-4));
    for (int i = 1; i < n; i++) {
        BOOST_TEST(smooth_sensitivity(g, i, 0, node_triangle_map[i],triangles, 1) == 0,
                   boost::test_tools::tolerance(1e-20));
    }
}

BOOST_AUTO_TEST_CASE(sensitivity_4) {
    Graph g = example_graph_4();
    const size_t n = boost::num_vertices(g);
    auto [triangles, node_triangle_map] = setup_node_triangle_map(g, n);

    BOOST_CHECK_EQUAL(global_sensitivity(g, 0, node_triangle_map[0], triangles), 2);
    for (int i = 1; i < n; i++) {
        BOOST_CHECK_EQUAL(global_sensitivity(g, i, node_triangle_map[i], triangles), 0);
    }

    BOOST_TEST(smooth_sensitivity(g, 0, 0, node_triangle_map[0], triangles,1) == 0,
               boost::test_tools::tolerance(1e-12));
    for (int i = 1; i < n; i++) {
        BOOST_TEST(smooth_sensitivity(g, i, 0, node_triangle_map[0], triangles,1) == 0,
                   boost::test_tools::tolerance(1e-12));
    }
}

BOOST_AUTO_TEST_CASE(sensitvity_5) {
    Graph g = example_graph_5();
    const size_t n = boost::num_vertices(g);
    auto [triangles, node_triangle_map] = setup_node_triangle_map(g, n);

    BOOST_CHECK_EQUAL(global_sensitivity(g, 0, node_triangle_map[0],triangles), 3);
    BOOST_CHECK_EQUAL(global_sensitivity(g, 1, node_triangle_map[1],triangles), 3);

    for (size_t i = 2; i < n; i++) {
        BOOST_CHECK_EQUAL(global_sensitivity(g, i, node_triangle_map[i], triangles), 0);
    }

    BOOST_TEST(smooth_sensitivity(g, 0, 0, node_triangle_map[0],triangles, 1) == 0.40600,
               boost::test_tools::tolerance(1e-4));
    BOOST_TEST(smooth_sensitivity(g, 1, 0, node_triangle_map[1], triangles,1) == 0.40600,
               boost::test_tools::tolerance(1e-4));
    for (int i = 2; i < n; i++) {
        BOOST_TEST(smooth_sensitivity(g, i, 0, node_triangle_map[i], triangles,1) == 0,
                   boost::test_tools::tolerance(1e-12));
    }
}

BOOST_AUTO_TEST_CASE(sensitivity_6) {
    Graph g = example_graph_6();
    const size_t n = boost::num_vertices(g);
    auto [triangles, node_triangle_map] = setup_node_triangle_map(g, n);

    BOOST_TEST(smooth_sensitivity(g, 0, 0, node_triangle_map[0],triangles, 0.1) == 1.81959,
               boost::test_tools::tolerance(1e-4));
    for (int i = 1; i < n; i++) {
        BOOST_TEST(smooth_sensitivity(g, i, 0, node_triangle_map[i], triangles,0.1) == 0,
                   boost::test_tools::tolerance(1e-12));
    }
}
