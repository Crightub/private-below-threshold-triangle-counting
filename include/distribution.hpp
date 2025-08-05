#pragma once

#include "boost/random/uniform_real_distribution.hpp"
#include <boost/random/linear_congruential.hpp>
#include "random"
#include "utils.hpp"

int sampleDiscreteLaplace(double lambda);
double sampleLaplace(double mu, double b);