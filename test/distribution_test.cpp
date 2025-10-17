#define BOOST_TEST_MODULE DistributionTest
#include <boost/graph/graph_concepts.hpp>
#include <boost/test/included/unit_test.hpp>

#include "distribution.h"
#include "graph_presets.h"
#include "global_sensitivity.h"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

inline void visualize_polynomial_tail(PolynomialTailRV &rv, int n_samples = 100000) {
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

    for (double s: samples) {
        int bin = std::min(static_cast<int>((s - min_x) / bin_width), bins - 1);
        hist_y[bin] += 1.0;
    }
    for (double &y: hist_y) y /= (n_samples * bin_width);
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
    plt::plot(hist_x_to_plot, hist_y_to_plot, {{"label", "Empirical PDF"}, {"color", "blue"}});

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

inline double discrete_laplace_pmf(int k, double p) {
    double norm = (1 - p) / (1 + p);
    return norm * std::pow(p, std::abs(k));
}

inline void visualize_discrete_laplace(double p, int n_samples = 100000) {
    std::vector<int> samples;
    samples.reserve(n_samples);
    for (int i = 0; i < n_samples; ++i) {
        samples.push_back(sample_discrete_laplace(p));
    }

    int min_x = *std::min_element(samples.begin(), samples.end());
    int max_x = *std::max_element(samples.begin(), samples.end());

    std::vector<int> hist_x;
    std::vector<double> hist_y;
    for (int k = min_x; k <= max_x; ++k) {
        hist_x.push_back(k);
        hist_y.push_back(0.0);
    }

    for (int s: samples) {
        hist_y[s - min_x] += 1.0;
    }

    for (double &y: hist_y)
        y /= n_samples;

    std::vector<double> pmf_x, pmf_y;
    for (int k = min_x; k <= max_x; ++k) {
        pmf_x.push_back(k);
        pmf_y.push_back(discrete_laplace_pmf(k, p));
    }

    plt::plot(hist_x, hist_y);

    plt::title("Discrete Laplace Distribution (λ = " + std::to_string(p) + ")");
    plt::legend();
    plt::show();
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

BOOST_AUTO_TEST_CASE(discrete_laplace_visualization) {
    const double eps = 0.1;
    visualize_discrete_laplace(std::exp(-eps));
}
