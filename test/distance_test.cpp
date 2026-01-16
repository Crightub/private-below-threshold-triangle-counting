#define BOOST_TEST_MODULE Distance_Test

#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/context.hpp>
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>

#include "distance_retrieval.h"
#include "utils.h"

// <-- Naive Implementations for k-th shortest distance and k-th sum -->
static int naive_kth(const std::vector<int> &p, int t, int k) {
    std::vector<int> d;
    for (int x: p)
        d.push_back(std::abs(x - t));
    std::nth_element(d.begin(), d.begin() + (k - 1), d.end());
    return d[k - 1];
}

static long long naive_sum_k(const std::vector<int> &p, int t, int k) {
    std::vector<int> d;
    for (int x: p)
        d.push_back(std::abs(x - t));
    std::nth_element(d.begin(), d.begin() + k, d.end());
    return std::accumulate(d.begin(), d.begin() + k, 0LL);
}

void test_single_distance(std::vector<int> p, SingleTargetDistance dr, int t) {
    for (int k = 1; k <= (int) p.size(); ++k) {
        auto [k_th_dis, l] = dr.k_th_closest_distance(k);
        int exp_kth_dis = naive_kth(p, t, k);

        std::cout << "k: " << k << " l: " << l << " k_th dis: " << k_th_dis << " exp dis: " << exp_kth_dis << std::endl;
        BOOST_CHECK_EQUAL(
            k_th_dis,
            exp_kth_dis
        );

        long long k_th_sum = dr.sum_k_closest_distances(k);
        long long exp_kth_sum = naive_sum_k(p, t, k);

        std::cout << "k: " << k << " k_th sum: " << k_th_sum << " exp sum: " << exp_kth_sum << std::endl;

        BOOST_CHECK_EQUAL(
            k_th_sum,
            exp_kth_sum
        );
    }
}

// <-- Unit Tests -->

BOOST_AUTO_TEST_CASE(basic) {
    std::vector<int> p = {1, 3, 7, 10, 12};
    int t = 6;
    SingleTargetDistance dr(p, t);
    test_single_distance(p, dr, t);
}

BOOST_AUTO_TEST_CASE(negative_values) {
    std::vector<int> p = {-5, -3, 0, 1, 3, 7, 10, 12};
    int t = 1;
    SingleTargetDistance dr(p, t);
    test_single_distance(p, dr, t);
}

BOOST_AUTO_TEST_CASE(duplicate_values) {
    std::vector<int> p = {-10, -5, -5, -5, 1, 4};
    int t = 1;
    SingleTargetDistance dr(p, t);
    test_single_distance(p, dr, t);
}

BOOST_AUTO_TEST_CASE(too_large_k) {
    std::vector<int> p = {-10, -5, -5, -5, 1, 4};
    int t = 1;
    SingleTargetDistance dr(p, t);

    auto [k_th_dis, l] = dr.k_th_closest_distance(10);
}


BOOST_AUTO_TEST_CASE(duplicate_boundary_values) {
    std::vector<int> p = {-98, -74, -73, -73, -71, -57, -26, -17, -9, 2, 7, 15, 20, 20, 27, 28, 33, 89, 89, 89};
    int t = 89;
    SingleTargetDistance dr(p, t);
    test_single_distance(p, dr, t);

    t = 90;
    dr.update_target(t);
    test_single_distance(p, dr, t);
}

BOOST_AUTO_TEST_CASE(negative_target) {
    std::vector<int> p = {-5, -3, 0, 1, 3, 7, 10, 12};
    int t = -4;
    SingleTargetDistance dr(p, t);
    test_single_distance(p, dr, t);
}

BOOST_AUTO_TEST_CASE(boundary_target) {
    std::vector<int> p = {-5, -3, 0, 1, 3, 7, 10, 12, 15, 20, 24};
    int t = -5;

    SingleTargetDistance dr(p, t);

    test_single_distance(p, dr, t);

    t = 24;
    dr.update_target(t);

    test_single_distance(p, dr, t);
}


BOOST_AUTO_TEST_CASE(target_updates) {
    std::vector<int> p = {-20, 0, 2, 20};
    SingleTargetDistance dr(p, -20);

    for (int t = -20; t <= 15; ++t) {
        dr.update_target(t);
        test_single_distance(p, dr, t);
    }
}

BOOST_AUTO_TEST_CASE(outside_range) {
    std::vector<int> p = {10, 20, 30, 40};

    int t = 0;
    SingleTargetDistance dr(p, t);
    test_single_distance(p, dr, t);

    t = 20;
    dr.update_target(t);
    test_single_distance(p, dr, t);

    t = 100;
    dr.update_target(t);
    test_single_distance(p, dr, t);
}

BOOST_AUTO_TEST_CASE(distance_retrieval_random) {
    std::mt19937 rng(42);

    for (int iter = 0; iter < 20; ++iter) {
        int n = 20;
        std::vector<int> p(n);
        std::uniform_int_distribution<int> dist(-100, 100);
        for (int i = 0; i < n; ++i)
            p[i] = dist(rng);

        std::sort(p.begin(), p.end());

        int t = dist(rng);
        SingleTargetDistance dr(p, t);

        test_single_distance(p, dr, t);

        // random monotone updates
        for (int step = 0; step < 10; ++step) {
            t += rng() % 5;
            dr.update_target(t);
            test_single_distance(p, dr, t);
        }
    }
}

// <-- Naive Implementations for double target k-th shortest distance and k-th sum -->
static int naive_double_kth(const std::vector<int> &p, int t, int k) {
    std::vector<int> d;
    for (int x: p) {
        int d1 = std::abs(x - (t - 1));
        int d2 = std::abs(x - (t + 1));
        d.push_back(std::min(d1, d2));
    }
    std::nth_element(d.begin(), d.begin() + (k - 1), d.end());
    return d[k - 1];
}

static long long naive_double_sum_k(const std::vector<int> &p, int t, int k) {
    std::vector<int> d;
    for (int x: p) {
        int d1 = std::abs(x - (t - 1));
        int d2 = std::abs(x - (t + 1));
        d.push_back(std::min(d1, d2));
    }
    std::nth_element(d.begin(), d.begin() + k, d.end());
    return std::accumulate(d.begin(), d.begin() + k, 0LL);
}

void test_double_distance(std::vector<int> p, DoubleTargetDistance dr, int t) {
    BOOST_TEST_CONTEXT(
        "p = " << vector_string(p) << ", t = " << t
    ) {
        for (int k = 1; k <= (int) p.size(); ++k) {
            std::cout << "Testing k: " << k << std::endl;
            auto [k_th_dis, l] = dr.k_th_closest_distance(k);
            int exp_kth_dis = naive_double_kth(p, t, k);

            BOOST_CHECK_EQUAL(k_th_dis, exp_kth_dis);

            long long k_th_sum = dr.sum_k_closest_distances(k);
            long long exp_kth_sum = naive_double_sum_k(p, t, k);

            BOOST_CHECK_EQUAL(k_th_sum, exp_kth_sum);
        }
    }
}


BOOST_AUTO_TEST_CASE(double_basic) {
    std::vector<int> p = {-69, -52, -50, -43, -43, -38, -35, -15, -15, -9, -3, 14, 24, 42, 49, 52, 79, 85, 89, 90};
    int t = 89;
    DoubleTargetDistance dr(p, t);
    test_double_distance(p, dr, t);
}

BOOST_AUTO_TEST_CASE(double_negative_values) {
    std::vector<int> p = {-5, -3, 0, 1, 3, 7, 10, 12, 13};
    int t = -5;
    DoubleTargetDistance dr(p, t);
    test_double_distance(p, dr, t);
}

BOOST_AUTO_TEST_CASE(double_boundary_target) {
    std::vector<int> p = {-5, -5, -3, 0, 1, 3, 7, 10, 12};
    int t = -5;
    DoubleTargetDistance dr(p, t);
    test_double_distance(p, dr, t);

    t = 12;
    dr.update_target(t);
    test_double_distance(p, dr, t);
}


BOOST_AUTO_TEST_CASE(double_target_edge_cases) {
    std::vector<int> p = {-10, -5, 0, 5, 10};

    // target far left
    {
        int t = -100;
        DoubleTargetDistance dr(p, t);
        test_double_distance(p, dr, t);
    }

    // target far right
    {
        int t = 100;
        DoubleTargetDistance dr(p, t);
        test_double_distance(p, dr, t);
    }
}


BOOST_AUTO_TEST_CASE(double_random) {
    std::mt19937 rng(std::chrono::system_clock::now().time_since_epoch().count());

    for (int iter = 0; iter < 20; ++iter) {
        int n = 20;
        std::vector<int> p(n);
        std::uniform_int_distribution<int> dist(-100, 100);
        for (int i = 0; i < n; ++i)
            p[i] = dist(rng);

        std::sort(p.begin(), p.end());

        int t = dist(rng);
        DoubleTargetDistance dr(p, t);

        test_double_distance(p, dr, t);

        // random monotone updates
        for (int step = 0; step < 10; ++step) {
            t += rng() % 5;
            dr.update_target(t);
            test_double_distance(p, dr, t);
        }
    }
}
