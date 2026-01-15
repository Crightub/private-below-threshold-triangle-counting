#include "utils.h"
#include <vector>

boost::minstd_rand rng(static_cast<unsigned int>(std::time(0)));

std::vector<Node> neighbors(const Graph &g, Node &v) {
    std::vector<Node> neighbors;
    for (std::pair<AdjIter, AdjIter> nbr = adjacent_vertices(v, g); nbr.first != nbr.second; ++nbr.first) {
        neighbors.push_back(*nbr.first);
    }

    return neighbors;
}


double unbiased_estimator(const int x, double p, int lambda) {
    if (x > lambda) {
        return 0;
    }
    if (x == lambda) {
        return -p / std::pow(1 - p, 2);
    }
    if (x == lambda - 1) {
        return 1 + p / std::pow(1 - p, 2);
    }
    if (x < lambda - 1) {
        return 1;
    }

    return 1;
}

int biased_estimator(int w1, int w2, int w3, int lambda) {
    const long long sum = w1 + w2 + w3;
    return sum < lambda;
}