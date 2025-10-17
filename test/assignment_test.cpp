#define BOOST_TEST_MODULE AssignmentTest

#include <boost/test/included/unit_test.hpp>

#include "counting.h"
#include "graph_generation.h"
#include "utils.h"
#include "graph_presets.h"
#include "graph_statistics.h"

BOOST_AUTO_TEST_CASE(find_triangles_1) {
    Graph g = complete_graph(20);

    auto triangles = find_triangles(g);

    BOOST_CHECK_EQUAL(triangles.size(), 1140);

    for (auto e: boost::make_iterator_range(edges(g))) {
        BOOST_CHECK_EQUAL(g[e].involved_triangles, 18);
    }
}

BOOST_AUTO_TEST_CASE(find_triangles_2) {
    Graph g = distinct_triangle_graph(20000);

    auto triangles = find_triangles(g);

    BOOST_CHECK_EQUAL(triangles.size(), 20000);

    // for (auto e: boost::make_iterator_range(edges(g))) {
    //     BOOST_CHECK_EQUAL(g[e].involved_triangles, 1);
    // }
}

BOOST_AUTO_TEST_CASE(find_triangles_3) {
    Graph g = star_outline_graph(20001);

    auto triangles= find_triangles(g);

    BOOST_CHECK_EQUAL(triangles.size(), 20000);
}


BOOST_AUTO_TEST_CASE(find_triangles_epinions) {
    std::string file_name = "../snap/epinions.csv";
    Graph g = load_snap_epinions_graph(file_name);

    auto triangles = find_triangles(g);

    BOOST_CHECK_EQUAL(triangles.size(), 4910076);
}

namespace std {
    template <>
    struct hash<Triangle> {
        size_t operator()(const Triangle &t) const noexcept {
            size_t h1 = std::hash<Node>()(t.nodes[0]);
            size_t h2 = std::hash<Node>()(t.nodes[1]);
            size_t h3 = std::hash<Node>()(t.nodes[2]);
            // Good hash combination
            return h1 ^ (h2 << 1) ^ (h3 << 2);
        }
    };
}

BOOST_AUTO_TEST_CASE(find_triangles_bitcoinotc) {
    std::string file_name = "../snap/bitcoinotc.csv";
    Graph g = load_snap_bitcoin_graph(file_name);

    std::cout << "Edge count: " << boost::num_edges(g) << std::endl;

    auto triangles = find_triangles(g);

    std::unordered_map<Triangle, int> triangle_counts;

    for (const auto &t : triangles)
        triangle_counts[t]++;

    bool found_duplicates = false;

    for (const auto &[t, count] : triangle_counts) {
        if (count > 1) {
            found_duplicates = true;
            std::cout << "Duplicate triangle (count=" << count << "): {"
                      << t.nodes[0] << ", " << t.nodes[1] << ", " << t.nodes[2]
                      << "}" << std::endl;
        }
    }

    if (!found_duplicates)
        std::cout << "✅ No duplicate triangles found!" << std::endl;


    BOOST_CHECK_EQUAL(triangles.size(), 33493);
}


void benchmark_assignment( Graph &g) {
    auto load_balanced_graph = g;
    auto no_load_balanced_graph = g;
    auto triangles_load_balanced = find_triangles(g);
    auto triangles_no_load_balanced = triangles_load_balanced;

    std::cout << "Load balancing ------------" << std::endl;
    for (auto t : triangles_load_balanced) {
        t.assign_triangle(load_balanced_graph, true);
    }

    auto load_balanced_c4 = compute_c4_instances(load_balanced_graph);
    std::cout << "#C4: " << load_balanced_c4 << std::endl;

    std::cout << "No Load balancing ------------" << std::endl;
    for (auto t : triangles_no_load_balanced) {
        t.assign_triangle(no_load_balanced_graph, false);
    }

    auto no_load_balanced_c4 = compute_c4_instances(no_load_balanced_graph);
    std::cout << "#C4: " << no_load_balanced_c4 << std::endl;
}

BOOST_AUTO_TEST_CASE(snap_bitcoin_assignment) {
    Graph g = load_snap_epinions_graph("../snap/bitcoinotc.csv");
    benchmark_assignment(g);
}

BOOST_AUTO_TEST_CASE(snap_epinions_assignment) {
    Graph g = load_snap_epinions_graph("../snap/epinions.csv");
    benchmark_assignment(g);
}