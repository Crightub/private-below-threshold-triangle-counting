#ifndef NEGATIVE_TRIANGLE_COUNTING_DISTRIBUTION_H
#define NEGATIVE_TRIANGLE_COUNTING_DISTRIBUTION_H

#include "boost/random/uniform_real_distribution.hpp"
#include <boost/random/cauchy_distribution.hpp>
#include "random"
#include "utils.h"
#include <boost/math/quadrature/tanh_sinh.hpp>

inline int sample_discrete_laplace(const double lambda)
{
    std::geometric_distribution<int> geom(1 - lambda);

    int u_sample = geom(rng);
    int v_sample = geom(rng);

    return u_sample - v_sample;
}

inline double sample_laplace(const double mu, const double b)
{
    const boost::random::uniform_real_distribution<> uniform(0.0, 1.0);
    double u = uniform(rng);

    return (u < 0.5) ? mu + b * log(2.0 * u) : mu - b * log(2 * (1 - u));
}


class PolynomialTailRV {
private:
    int gamma;
    double leading_coefficient;

    boost::random::uniform_real_distribution<double> uniform;

public:
    explicit PolynomialTailRV(int gamma_)
        : gamma(gamma_), uniform(0.0, 1.0)
    {
        using boost::math::quadrature::tanh_sinh;

        auto f = [&](double y) {
            return 1.0 / (1.0 + std::pow(std::abs(y), gamma));
        };

        tanh_sinh<double> integrator;
        leading_coefficient = integrator.integrate(f, -INFINITY, INFINITY);
    }

    double pdf(double x) const {
        return 1.0 / (leading_coefficient * (1.0 + std::pow(std::abs(x), gamma)));
    }

    double sample() {
        boost::random::cauchy_distribution<double> cauchy(0.0, 1.0);

        auto g = [](double y) {
            return 1.0 / (M_PI * (1.0 + y * y));
        };
        const double M = 2.0;

        while (true) {
            double y = cauchy(rng);
            double u = uniform(rng);

            if (u < pdf(y) / (M * g(y))) {
                return y;
            }
        }
    }
};

#endif //NEGATIVE_TRIANGLE_COUNTING_DISTRIBUTION_H