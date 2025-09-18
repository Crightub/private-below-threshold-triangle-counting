#ifndef NEGATIVE_TRIANGLE_COUNTING_UTILS_H
#define NEGATIVE_TRIANGLE_COUNTING_UTILS_H

#include <vector>
#include <iostream>
#include <Eigen/Sparse>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/log/trivial.hpp>
#include <sstream>

struct edge_info {
    int weight;
    int noise;
};

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, edge_info> Graph;
typedef boost::sorted_erdos_renyi_iterator<boost::minstd_rand, Graph> ERGen;

extern boost::minstd_rand rng;

double unbiased_estimator(int x, double p);

int biased_estimator(int w1, int w2, int w3);

double approximate_variance(double w1, double w2, double w3, double lambda);

struct DecisionVariableInfo {
    int matrixIndex;
    int countingNodeId;
    double indicator;
};

class Triangle {
public:
    int source_id;
    int lower_id;
    int higher_id;

    Triangle(const int sId, const int Id1, const int Id2) : source_id(sId) {
        higher_id = Id1 >= Id2 ? Id1 : Id2;
        lower_id = Id1 < Id2 ? Id1 : Id2;
    }

    friend std::ostream &operator<<(std::ostream &os, const Triangle &t);

    bool operator<(const Triangle &other) const {
        std::array<int, 3> thisIds = {source_id, lower_id, higher_id};
        std::array<int, 3> otherIds = {other.source_id, other.lower_id, other.higher_id};
        std::sort(thisIds.begin(), thisIds.end());
        std::sort(otherIds.begin(), otherIds.end());
        return thisIds < otherIds;
    }
};

std::ostream &operator<<(std::ostream &os, const Triangle &t);

std::list<Triangle> find_rooted_triangles(const Graph *g, const int root_id);

inline void print_vector(const std::vector<double> &v) {
    std::ostringstream oss;
    oss << "[";
    for (int j = 0; j < v.size(); j++) {
        oss << v[j];
        if (j < v.size() - 1) {
            oss << ", ";
        }
    }
    oss << "]";
    std::cout << oss.str() << std::endl;
    // BOOST_LOG_TRIVIAL(info) << oss.str();
}

inline void print_matrix(const Eigen::VectorXd &A) {
    std::ostringstream oss;
    for (int i = 0; i < A.size(); ++i) {
        oss << A[i];
        if (i < A.size() - 1) {
            oss << "\t";
        }
    }
    BOOST_LOG_TRIVIAL(info) << oss.str();
}

inline void print_matrix(const std::vector<std::vector<double> > &A) {
    std::ostringstream oss;
    for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A[i].size(); ++j) {
            oss << A[i][j];
            if (j < A[i].size() - 1) {
                oss << "\t";
            }
        }
        if (i < A.size() - 1) {
            oss << "\n";
        }
    }
    BOOST_LOG_TRIVIAL(info) << "\n" << oss.str();
}

inline void print_matrix(const Eigen::SparseMatrix<double> &A) {
    const int max_size = 15;
    std::ostringstream oss;

    if (A.rows() <= max_size && A.cols() <= max_size) {
        Eigen::MatrixXd denseA(A);
        oss << denseA; // Eigen already has stream operator
    } else {
        for (int k = 0; k < A.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
                oss << "(" << it.row() << ", " << it.col() << ") = " << it.value() << "\n";
            }
        }
    }

    BOOST_LOG_TRIVIAL(info) << "\n" << oss.str();
}

#endif