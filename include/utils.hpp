#pragma once

#include <vector>
#include <iostream>
#include <Eigen/Sparse>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/random/linear_congruential.hpp>

struct edge_info
{
    int weight;
    double noise;
};

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, edge_info> Graph;
typedef boost::sorted_erdos_renyi_iterator<boost::minstd_rand, Graph> ERGen;

extern boost::minstd_rand rng;


void printMatrix(std::vector<std::vector<double>> A);
void printMatrix(Eigen::SparseMatrix<double> A);
void printMatrix(Eigen::VectorXd A);

double unbiased_estimator(int x, double p);

int biased_estimator(int w1, int w2, int w3);

double approximate_variance(double w1, double w2, double w3, double lambda);

struct DecisionVariableInfo
{
    int matrixIndex;
    int countingNodeId;
    double indicator;
};

class Triangle
{
public:
    int sourceId;
    int leftTargetId;
    int rightTargetId;
    double variance;
    double indicator;

    Triangle(int sId, int lId, int rId, double var, double ind) : sourceId(sId), leftTargetId(lId), rightTargetId(rId), variance(var), indicator(ind) {}

    friend std::ostream &operator<<(std::ostream &os, const Triangle &t);

    bool operator<(const Triangle &other) const
    {
        std::array<int, 3> thisIds = {sourceId, leftTargetId, rightTargetId};
        std::array<int, 3> otherIds = {other.sourceId, other.leftTargetId, other.rightTargetId};
        std::sort(thisIds.begin(), thisIds.end());
        std::sort(otherIds.begin(), otherIds.end());
        return thisIds < otherIds;
    }
};

std::ostream &operator<<(std::ostream &os, const Triangle &t);