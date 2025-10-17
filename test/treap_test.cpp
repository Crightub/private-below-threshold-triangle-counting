#define BOOST_TEST_MODULE TreapTest
#include <boost/test/included/unit_test.hpp>

#include "utils.h"
#include "graph_presets.h"
#include "smooth_sensitivity.h"
#include "treap.h"


// <-- Utility Methods for Tests -->
TreapNode *setup_treap(const std::vector<double> &values) {
    TreapNode *tree = nullptr;

    for (int i = 0; i < values.size(); i++) {
        tree = insert(tree, values[i]);
    }

    return tree;
}

BOOST_AUTO_TEST_CASE(treap_1) {
    TreapNode *t = nullptr;
    t = insert(t, -5);
    t = insert(t, -3);
    t = insert(t, -1);
    t = insert(t, 3);
    t = insert(t, 4);
    t = insert(t, 5);

    BOOST_CHECK_EQUAL(get_size(t), 6);
    BOOST_CHECK_EQUAL(kth(t, 3)->key, -1.0);
    BOOST_CHECK_EQUAL(kth(t, 6)->key, 5);
    BOOST_CHECK_EQUAL(sum_first_k(t, 3), -9.0);

    t = erase(t, -1);

    BOOST_CHECK_EQUAL(get_size(t), 5);
    BOOST_CHECK_EQUAL(kth(t, 3)->key, 3.0);
    BOOST_CHECK_EQUAL(sum_first_k(t, 3), -5.0);
}

BOOST_AUTO_TEST_CASE(kth_smallest_distance_1) {
    TreapNode *l = setup_treap(std::vector<double>{0, -2});
    double offL = 2;

    TreapNode *r = setup_treap(std::vector<double>{4, 8, 9, 10});
    double offR = -offL;

    BOOST_CHECK_EQUAL(find_k_smallest_distance(l, r, 1, offL, offR), 0);
    BOOST_CHECK_EQUAL(find_k_smallest_distance(l, r, 2, offL, offR), 2.0);
    BOOST_CHECK_EQUAL(find_k_smallest_distance(l, r, 3, offL, offR), 2.0);
    BOOST_CHECK_EQUAL(find_k_smallest_distance(l, r, 4, offL, offR), 6.0);
    BOOST_CHECK_EQUAL(find_k_smallest_distance(l, r, 5, offL, offR), 7.0);
    BOOST_CHECK_EQUAL(find_k_smallest_distance(l, r, 6, offL, offR), 8.0);
}


BOOST_AUTO_TEST_CASE(kth_smallest_distance_2) {
    TreapNode *l = setup_treap(std::vector<double>{0, -2, -1, -4});
    double offL = 7;

    TreapNode *r = setup_treap(std::vector<double>{7, 8, 9, 10});
    double offR = -offL;

    BOOST_CHECK_EQUAL(find_k_smallest_distance(l, r, 1, offL, offR), 0);
    BOOST_CHECK_EQUAL(find_k_smallest_distance(l, r, 2, offL, offR), 1.0);
    BOOST_CHECK_EQUAL(find_k_smallest_distance(l, r, 3, offL, offR), 2.0);
    BOOST_CHECK_EQUAL(find_k_smallest_distance(l, r, 4, offL, offR), 3.0);
    BOOST_CHECK_EQUAL(find_k_smallest_distance(l, r, 5, offL, offR), 3.0);
    BOOST_CHECK_EQUAL(find_k_smallest_distance(l, r, 6, offL, offR), 5.0);
    BOOST_CHECK_EQUAL(find_k_smallest_distance(l, r, 7, offL, offR), 6.0);
    BOOST_CHECK_EQUAL(find_k_smallest_distance(l, r, 8, offL, offR), 7.0);
}

BOOST_AUTO_TEST_CASE(find_optimal_t_1) {
    TreapNode *l = nullptr;
    TreapNode *r = nullptr;

    l = insert(l, 0);
    l = insert(l, -2);

    r = insert(r, 4);
    r = insert(r, 8);
    r = insert(r, 9);
    r = insert(r, 10);

    double offL = 2;
    double offR = -2;

    std::vector<double> centers = std::vector<double>{-5, -3, -1, 3, 4, 5};

    int opt_t = optimal_shift_set_size(l, r, offL, offR, 1.0 / 3);
    BOOST_CHECK_EQUAL(opt_t, 2);
}


BOOST_AUTO_TEST_CASE(find_optimal_t_2) {
    TreapNode *l = nullptr;
    TreapNode *r = nullptr;

    l = insert(l, 0);
    l = insert(l, -2);
    l = insert(l, -4);

    r = insert(r, 8);
    r = insert(r, 9);
    r = insert(r, 10);

    double offL = 4;
    double offR = -4;


    BOOST_CHECK_EQUAL(optimal_shift_set_size(l, r, offL, offR, 0.25), 2);
    BOOST_CHECK_EQUAL(optimal_shift_set_size(l, r, offL, offR, 1.0/20), 4);
}

BOOST_AUTO_TEST_CASE(find_optimal_t_all_left) {
    TreapNode *l = setup_treap(std::vector<double>{0, -2, -4, -8, -9, -10});
    double offL = 10;
    TreapNode *r = nullptr;
    double offR = -offL;


    BOOST_CHECK_EQUAL(optimal_shift_set_size(l, r, offL, offR, 1.0/20), 3);
}

BOOST_AUTO_TEST_CASE(same_keys) {
    TreapNode *l = nullptr;

    l = insert(l, 0);
    l = insert(l, 0);
    l = insert(l, 0);
    l = insert(l, 0);

    print_treap(l);

    BOOST_CHECK_EQUAL(get_size(l), 4);
}

BOOST_AUTO_TEST_CASE(advance_center_1) {
    TreapNode *l = nullptr;
    TreapNode *r = nullptr;

    l = insert(l, 0);
    l = insert(l, -2);

    double offL = 2;

    r = insert(r, 8);
    r = insert(r, 9);
    r = insert(r, 10);

    double offR = -2;

    std::vector<double> centers = std::vector<double>{-5, -3, -1, 3, 4, 5};

    std::tie(l, r) = advance_center(l, r, offL, offR, centers, 1);

    BOOST_CHECK_EQUAL(offL, 4);
    BOOST_CHECK_EQUAL(offR, -4);

    BOOST_CHECK_EQUAL(get_size(l), 2);
    BOOST_CHECK_EQUAL(get_size(r), 3);
    print_treap(l);

    print_treap(r);
}

BOOST_AUTO_TEST_CASE(advance_center_2) {
    TreapNode *l = nullptr;
    TreapNode *r = nullptr;

    l = insert(l, 0);
    l = insert(l, -2);

    double offL = 2;

    r = insert(r, 8);
    r = insert(r, 9);
    r = insert(r, 10);

    double offR = -2;

    std::vector<double> centers = std::vector<double>{-5, -3, -1, 3, 4, 5};

    std::tie(l, r) = advance_center(l, r, offL, offR, centers, 1);

    BOOST_CHECK_EQUAL(offL, 4);
    BOOST_CHECK_EQUAL(offR, -4);

    BOOST_CHECK_EQUAL(get_size(l), 2);
    BOOST_CHECK_EQUAL(get_size(r), 3);
    print_treap(l);

    print_treap(r);
}


BOOST_AUTO_TEST_CASE(find_sum_distance_1) {
    TreapNode *r = nullptr;
    TreapNode *l = nullptr;

    l = insert(l, 0);
    l = insert(l, -2);

    r = insert(r, 4);
    r = insert(r, 8);
    r = insert(r, 9);
    r = insert(r, 10);

    double offL = 2;
    double offR = -2;

    // BOOST_CHECK_EQUAL(find_k_sum_distance(l, r, 1, offL, offR), 0.0);
    // BOOST_CHECK_EQUAL(find_k_sum_distance(l, r, 2, offL, offR), 2.0);
    // BOOST_CHECK_EQUAL(find_k_sum_distance(l, r, 3, offL, offR), 4.0);
    // BOOST_CHECK_EQUAL(find_k_sum_distance(l, r, 4, offL, offR), 10.0);
    // BOOST_CHECK_EQUAL(find_k_sum_distance(l, r, 5, offL, offR), 17.0);
    // BOOST_CHECK_EQUAL(find_k_sum_distance(l, r, 6, offL, offR), 25.0);
}


BOOST_AUTO_TEST_CASE(find_sum_distance_2) {
    TreapNode *r = nullptr;
    TreapNode *l = nullptr;

    l = insert(l, 0);
    l = insert(l, -2);
    l = insert(l, -4);

    r = insert(r, 8);
    r = insert(r, 9);
    r = insert(r, 10);

    double offL = 4;
    double offR = -4;

    // BOOST_CHECK_EQUAL(find_k_sum_distance(l, r, 1, offL, offR), 0.0);
    // BOOST_CHECK_EQUAL(find_k_sum_distance(l, r, 2, offL, offR), 2.0);
    // BOOST_CHECK_EQUAL(find_k_sum_distance(l, r, 3, offL, offR), 6.0);
    // BOOST_CHECK_EQUAL(find_k_sum_distance(l, r, 4, offL, offR), 10.0);
    // BOOST_CHECK_EQUAL(find_k_sum_distance(l, r, 5, offL, offR), 15.0);
    // BOOST_CHECK_EQUAL(find_k_sum_distance(l, r, 6, offL, offR), 21.0);
}

BOOST_AUTO_TEST_CASE(find_sum_distance_all_left) {
    TreapNode *r = nullptr;
    TreapNode *l = nullptr;

    l = insert(l, 0);
    l = insert(l, -2);
    l = insert(l, -4);
    l = insert(l, -8);
    l = insert(l, -9);
    l = insert(l, -10);

    r = insert(r, 8);
    r = erase(r, 8);

    print_treap(l);

    double offL = 10;
    double offR = -10;

    BOOST_CHECK_EQUAL(compute_sum(l, r, offL, offR, 1, 0), 0);
    BOOST_CHECK_EQUAL(compute_sum(l, r, offL, offR, 2, 0), 1);
    BOOST_CHECK_EQUAL(compute_sum(l, r, offL, offR, 3, 0), 3);
    BOOST_CHECK_EQUAL(compute_sum(l, r, offL, offR, 4, 0), 9);
    BOOST_CHECK_EQUAL(compute_sum(l, r, offL, offR, 5, 0), 17);
    BOOST_CHECK_EQUAL(compute_sum(l, r, offL, offR, 6, 0), 27);
}

BOOST_AUTO_TEST_CASE(fixed_edge_opt_1) {
    const Graph g = example_graph_1();
    double opt;
    // opt = fixed_edge_sensitivity(g, 0, 0, 1.0);
    //
    // BOOST_CHECK_EQUAL(opt, 1);

    opt = fixed_edge_sensitivity(g, 0, 0, 1.0 / 20, false);

    BOOST_TEST(opt == 2.42612, boost::test_tools::tolerance(1e-4));
    std::cout << "opt: " << opt << std::endl;
}


// TODO: Add test that for the optimal value does not change anything

// TODO: Add test that
