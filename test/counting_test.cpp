#define BOOST_TEST_MODULE CountingTest
#include <boost/graph/graph_concepts.hpp>
#include <boost/test/included/unit_test.hpp>

#include "basic_counting.h"
#include "distribution.h"
#include "graph_presets.h"
#include "global_sensitivity.h"
#include "graph_generation.h"
#include "smooth_sensitivity.h"
#include "matplotlibcpp.h"


// <!-- counting test utility functions -->
std::tuple<std::vector<std::vector<bool> >, std::vector<std::vector<int> >, std::vector<std::vector<int> >, int>
setup_test_matrices(const Graph &g) {
    std::size_t n = num_vertices(g);
    std::vector<std::vector<bool> > A(n, std::vector<bool>(n, false));
    std::vector<std::vector<int> > W(n, std::vector<int>(n, 0));
    std::vector<std::vector<int> > N(n, std::vector<int>(n, 0));
    setup_graph_matrices(g, A, W, N);

    return {A, W, N, n};
}

namespace plt = matplotlibcpp;

inline void visualize_polynomial_tail(PolynomialTailRV& rv, int n_samples = 100000) {
    std::vector<double> samples;
    samples.reserve(n_samples);
    for (int i = 0; i < n_samples; ++i) {
        double s = rv.sample();
        // BOOST_LOG_TRIVIAL(debug) << "sample: " << s;
        samples.push_back(s);
    }

    int bins = 1000;
    double min_x = *std::min_element(samples.begin(), samples.end());
    double max_x = *std::max_element(samples.begin(), samples.end());
    double bin_width = (max_x - min_x) / bins;

    std::vector<double> hist_y(bins, 0.0);
    std::vector<double> hist_x(bins);

    for (double s : samples) {
        int bin = std::min(static_cast<int>((s - min_x) / bin_width), bins - 1);
        hist_y[bin] += 1.0;
    }
    for (double& y : hist_y) y /= (n_samples * bin_width);
    for (int i = 0; i < bins; ++i) {
        hist_x[i] = min_x + (i + 0.5) * bin_width;
    }


    // plot only bins in [-10,10]
    std::vector<double> hist_x_to_plot;
    std::vector<double> hist_y_to_plot;
    for (int i = 0; i < bins; ++i) {
        double center = min_x + (i + 0.5) * bin_width;
        if (center >= -10 && center <= 10) {
            hist_x_to_plot.push_back(center);
            hist_y_to_plot.push_back(hist_y[i]);
        }
    }
    plt::plot(hist_x_to_plot, hist_y_to_plot, {{"label","Empirical PDF"},{"color","blue"}});

    std::vector<double> xs, ys;
    double step = (10 - -10) / 500.0;
    for (double x = -10; x <= 10; x += step) {
        xs.push_back(x);
        ys.push_back(rv.pdf(x));
    }
    plt::plot(xs, ys, {{"label", "Analytical PDF"}, {"color", "red"}});

    plt::title("PolynomialTailRV visualization");
    plt::legend();
    plt::show();
}

BOOST_AUTO_TEST_CASE(sensitivity_1) {
    const Graph g = example_graph_1();
    auto [A, W, N, n] = setup_test_matrices(g);

    BOOST_CHECK_EQUAL(global_sensitivity(A, 0), 6);
    for (int i = 1; i < n; i++) {
        BOOST_CHECK_EQUAL(global_sensitivity(A, i), 0);
    }

    BOOST_CHECK_EQUAL(smooth_sensitivity(g, 0, 1), 1);
    for (int i = 1; i < n; i++) {
        BOOST_CHECK_EQUAL(smooth_sensitivity(g, i, 1), 0);
    }
}

BOOST_AUTO_TEST_CASE(sensitivity_2) {
    const Graph g = example_graph_2();
    auto [A, W, N, n] = setup_test_matrices(g);

    BOOST_CHECK_EQUAL(global_sensitivity(A, 0), 2);
    for (int i = 1; i < n; i++) {
        BOOST_CHECK_EQUAL(global_sensitivity(A, i), 0);
    }

    BOOST_CHECK_EQUAL(smooth_sensitivity(g, 0, 1), 1);
    for (int i = 1; i < n; i++) {
        BOOST_CHECK_EQUAL(smooth_sensitivity(g, i, 1), 0);
    }
}

BOOST_AUTO_TEST_CASE(sensitivity_3) {
    const Graph g = example_graph_3();
    auto [A, W, N, n] = setup_test_matrices(g);

    BOOST_CHECK_EQUAL(global_sensitivity(A, 0), 1);
    for (int i = 1; i < n; i++) {
        BOOST_CHECK_EQUAL(global_sensitivity(A, i), 0);
    }

    BOOST_TEST(smooth_sensitivity(g, 0, 1) == 1, boost::test_tools::tolerance(1e-4));
    for (int i = 1; i < n; i++) {
        BOOST_TEST(smooth_sensitivity(g, i, 1) == 0, boost::test_tools::tolerance(1e-20));
    }
}

BOOST_AUTO_TEST_CASE(sensitivity_4) {
    const Graph g = example_graph_4();
    auto [A, W, N, n] = setup_test_matrices(g);

    BOOST_CHECK_EQUAL(global_sensitivity(A, 0), 2);
    for (int i = 1; i < n; i++) {
        BOOST_CHECK_EQUAL(global_sensitivity(A, i), 0);
    }

    BOOST_TEST(smooth_sensitivity(g, 0, 1) == 0, boost::test_tools::tolerance(1e-12));
    for (int i = 1; i < n; i++) {
        BOOST_TEST(smooth_sensitivity(g, i, 1) == 0, boost::test_tools::tolerance(1e-12));
    }
}

BOOST_AUTO_TEST_CASE(sensitvity_5) {
    const Graph g = example_graph_5();
    auto [A, W, N, n] = setup_test_matrices(g);

    BOOST_CHECK_EQUAL(global_sensitivity(A, 0), 3);
    BOOST_CHECK_EQUAL(global_sensitivity(A, 1), 3);

    for (int i = 2; i < n; i++) {
        BOOST_CHECK_EQUAL(global_sensitivity(A, i), 0);
    }

    BOOST_TEST(smooth_sensitivity(g, 0, 1) == 0.40600, boost::test_tools::tolerance(1e-4));
    BOOST_TEST(smooth_sensitivity(g, 1, 1) == 0.40600, boost::test_tools::tolerance(1e-4));
    for (int i = 2; i < n; i++) {
        BOOST_TEST(smooth_sensitivity(g, i, 1) == 0, boost::test_tools::tolerance(1e-12));
    }
}

BOOST_AUTO_TEST_CASE(sensitivity_6) {
    const Graph g = example_graph_6();
    auto [A, W, N, n] = setup_test_matrices(g);

    BOOST_TEST(smooth_sensitivity(g, 0, 0.1) == 2.01096, boost::test_tools::tolerance(1e-4));
    for (int i = 1; i < n; i++) {
        BOOST_TEST(smooth_sensitivity(g, i, 0.1) == 0, boost::test_tools::tolerance(1e-12));
    }
}


BOOST_AUTO_TEST_CASE(apply_global_sensitivity_1) {
    const Graph g = example_graph_1();

    PrivateCountingConfig cfg_unbiased{true, 0.1, 1, false};
    PrivateCountingConfig cfg_biased{false, 0.1, 1, false};

    auto [A, W, N, n] = setup_test_matrices(g);
    std::vector<double> unbiased_noise(n, 0.0);
    std::vector<double> biased_noise(n, 0.0);

    apply_global_sensitivity(A, cfg_unbiased, unbiased_noise);
    apply_global_sensitivity(A, cfg_biased, biased_noise);

    BOOST_LOG_TRIVIAL(debug) << "unbiased noise: ";
    print_vector(unbiased_noise);

    BOOST_TEST(unbiased_noise[0] == -12061.9, boost::test_tools::tolerance(1e-1));

    BOOST_LOG_TRIVIAL(debug) << "biased noise: ";
    print_vector(biased_noise);

    BOOST_TEST(biased_noise[0] == -3.85679, boost::test_tools::tolerance(1e-3));
}

BOOST_AUTO_TEST_CASE(polynomial_tail_rv) {
    const int gamma = 3;
    PolynomialTailRV rv(gamma);

    // 1. Check normalization of PDF
    using boost::math::quadrature::tanh_sinh;
    tanh_sinh<double> integrator;
    auto f = [&](double x) { return rv.pdf(x); };

    double integral = integrator.integrate(f, -INFINITY, INFINITY);
    BOOST_CHECK_CLOSE(integral, 1.0, 1e-6);

    // 2. Check symmetry of PDF
    BOOST_CHECK_CLOSE(rv.pdf(1.23), rv.pdf(-1.23), 1e-12);

    // 3. Generate samples and check empirical mean is close to 0
    const int N = 100000;
    double sum = 0.0;
    for (int i = 0; i < N; ++i) {
        sum += rv.sample();
    }
    double mean = sum / N;
    BOOST_CHECK_SMALL(mean, 0.05);

    // 4. Check that probability mass in tails decreases with |x|
    double center = rv.pdf(0.0);
    double mid = rv.pdf(5.0);
    double tail = rv.pdf(20.0);

    BOOST_CHECK(center > mid);
    BOOST_CHECK(mid > tail);
}

BOOST_AUTO_TEST_CASE(polynomial_tail_rv_visualization) {
    const int gamma = 4;
    PolynomialTailRV rv(gamma);

    visualize_polynomial_tail(rv, 100000);
}

BOOST_AUTO_TEST_CASE(apply_smooth_sensitivity_1) {
    const Graph g = example_graph_1();

    PrivateCountingConfig cfg_unbiased{true, 1, 1, true};
    PrivateCountingConfig cfg_biased{false, 1, 1, true};

    std::vector<double> unbiased_noise(g.vertex_set().size(), 0.0);
    std::vector<double> biased_noise(g.vertex_set().size(), 0.0);

    apply_smooth_sensitivity(g, cfg_unbiased, unbiased_noise);
    apply_smooth_sensitivity(g, cfg_biased, biased_noise);

    BOOST_LOG_TRIVIAL(debug) << "unbiased noise: ";
    print_vector(unbiased_noise);

    // BOOST_TEST(unbiased_noise[0] == -12061.9, boost::test_tools::tolerance(1e-1));

    BOOST_LOG_TRIVIAL(debug) << "biased noise: ";
    print_vector(biased_noise);

    // BOOST_TEST(biased_noise[0] == -3.85679, boost::test_tools::tolerance(1e-3));
}


BOOST_AUTO_TEST_CASE(apply_smooth_sensitivity_2) {
    const Graph g = generate_graph(500, 0.1, 0, 5, 1);

    PrivateCountingConfig cfg_unbiased{true, 1, 100, true};
    PrivateCountingConfig cfg_biased{false, 1, 100, true};

    // std::vector<double> unbiased_noise(g.vertex_set().size(), 0.0);
    std::vector<double> biased_noise(g.vertex_set().size(), 0.0);

    // apply_smooth_sensitivity(g, cfg_unbiased, unbiased_noise);
    apply_smooth_sensitivity(g, cfg_biased, biased_noise);

    // BOOST_LOG_TRIVIAL(debug) << "unbiased noise: ";
    // print_vector(unbiased_noise);

    // BOOST_TEST(unbiased_noise[0] == -12061.9, boost::test_tools::tolerance(1e-1));

    BOOST_LOG_TRIVIAL(debug) << "biased noise: ";
    print_vector(biased_noise);

    // BOOST_TEST(biased_noise[0] == -3.85679, boost::test_tools::tolerance(1e-3));
}

BOOST_AUTO_TEST_CASE(privacy_counting_1) {
    const Graph g = example_graph_2();
    auto res = randomized_private_counting(g, {false, 1, 1, false});
}
