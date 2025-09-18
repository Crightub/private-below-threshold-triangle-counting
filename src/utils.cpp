#include "utils.h"
#include <iostream>
#include <vector>

// boost::minstd_rand rng(static_cast<unsigned int>(std::time(0)));
boost::minstd_rand rng(static_cast<unsigned int>(0));

double unbiased_estimator(const int x, double p) {
    if (x > 0) {
        return 0;
    }
    if (x == 0) {
        return -p / std::pow((1 - p), 2);
    }
    if (x == -1) {
        return 1 + p / std::pow((1 - p), 2);
    }
    if (x < -1) {
        return 1;
    }

    return 1;
}

int biased_estimator(int w1, int w2, int w3) {
    return w1 + w2 + w3 < 0;
}

double approximate_variance(double const w1, double const w2, double const w3, double const lambda) {
    return std::pow(lambda, std::abs(-w1 - w2 - w3) + 1);
}

std::ostream &operator<<(std::ostream &out, const Triangle &t) {
    out << t.source_id << " - " << t.lower_id << " - " << t.higher_id;
    return out;
}

std::list<Triangle> find_rooted_triangles(const Graph *g, const int root_id) {
    std::list<Triangle> triangles;

    const Graph::vertex_descriptor root = boost::vertex(root_id, *g);

    std::vector<int> neighbors;
    for (auto [nbrIt, nbrEnd] = boost::adjacent_vertices(root, *g); nbrIt != nbrEnd; ++nbrIt) {
        neighbors.push_back(*nbrIt);
    }

    for (size_t i = 0; i < neighbors.size(); ++i) {
        for (size_t j = i + 1; j < neighbors.size(); ++j) {
            int u = neighbors[i];
            int v = neighbors[j];
            if (boost::edge(u, v, *g).second) {
                triangles.emplace_back(root_id, u, v);
            }
        }
    }

    return triangles;
}
