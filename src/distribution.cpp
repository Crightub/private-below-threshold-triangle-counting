#include "distribution.hpp"

int sampleDiscreteLaplace(double lambda)
{
    std::geometric_distribution<int> geom(1 - lambda);

    int u_sample = geom(rng);
    int v_sample = geom(rng);

    return u_sample - v_sample;
}

double sampleLaplace(double mu, double b)
{
    boost::random::uniform_real_distribution<> uniform(0.0, 1.0);
    double u = uniform(rng);

    return (u < 0.5) ? mu + b * log(2.0 * u) : mu - b * log(2 * (1 - u));
}