#define BOOST_TEST_MODULE CountingTest

#include <boost/test/included/unit_test.hpp>

#include "distribution.h"

BOOST_AUTO_TEST_CASE(test_local_counting) {

}

BOOST_AUTO_TEST_CASE(test_unbiased_estimator) {
    const int n = 20;
    const double eps = 0.5;
    const int lambda = 0;

    const double p = std::exp(-eps);

    double count = 0;
    for (int i = 0; i < n; i++) {
        const int sample_noise = sample_discrete_laplace(p);
        double res = unbiased_estimator(sample_noise, p, lambda);
        count += res;
        std::cout << res << std::endl;
    }

    std::cout << "Average Value: " << count/n << std::endl;
}
