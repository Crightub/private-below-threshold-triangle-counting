#include "utils.hpp"
#include "qp_counting.hpp"
#include "distribution.hpp"
#include <map>
#include <list>
#include <vector>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <boost/property_map/property_map.hpp>
#include "../osqp-cpp/include/osqp++.h"

using TriangleStats = std::tuple<std::list<Triangle>, std::map<std::pair<int, int>, std::list<Triangle>>>;
using ObjectiveMatrixStats = std::tuple<Eigen::SparseMatrix<double>, std::map<Triangle, std::list<DecisionVariableInfo>>>;

TriangleStats setupTriangleStatsForQP(Graph *g, double eps)
{
    using Vertex = typename boost::graph_traits<Graph>::vertex_descriptor;
    using Edge = typename boost::graph_traits<Graph>::edge_descriptor;

    double lambda = std::exp(-eps / 2.0);

    std::size_t n = num_vertices(*g);

    auto index_map = get(boost::vertex_index, *g);
    auto weight_map = get(&edge_info::weight, *g);
    auto noise_map = get(&edge_info::noise, *g);

    std::map<std::pair<int, int>, std::list<Triangle>> edgeTriangleMap{};
    std::list<Triangle> triangles = {};

    for (auto e : boost::make_iterator_range(edges(*g)))
    {
        Vertex u = source(e, *g);
        Vertex v = target(e, *g);

        int ui = get(index_map, u);
        int vi = get(index_map, v);

        double uv_noisy_weight = get(weight_map, e) + get(noise_map, e);

        for (auto w : boost::make_iterator_range(adjacent_vertices(u, *g)))
        {
            if (w == v)
                continue;

            auto wv_edge = edge(w, v, *g);
            if (!wv_edge.second)
            {
                continue;
            }

            int wi = get(index_map, w);
            auto wu_edge = edge(w, u, *g);

            double wu_weight = get(weight_map, wu_edge.first);
            double wv_weight = get(weight_map, wv_edge.first);

            double wu_noisy_weight = get(noise_map, wu_edge.first);
            double wv_noisy_weight = get(noise_map, wv_edge.first);

            double var = approximate_variance(wu_noisy_weight, wv_noisy_weight, uv_noisy_weight, lambda);
            double ind = unbiased_estimator(-wu_weight - wv_weight - uv_noisy_weight, lambda);
            Triangle t(wi, ui, vi, var, ind);
            triangles.push_back(t);
            edgeTriangleMap[std::pair<int, int>{ui, vi}].push_back(t);
        }
    }

    return TriangleStats{triangles, edgeTriangleMap};
}

ObjectiveMatrixStats setupObjectiveMatrix(std::list<Triangle> triangles, std::map<std::pair<int, int>, std::list<Triangle>> edgeTriangleMap)
{
    int triangleCount = triangles.size() / 3;
    int objectiveMatrixRowCount = 3 * triangleCount;
    int objectiveMatrixColumnCount = 3 * triangleCount;
    std::vector<Eigen::Triplet<double>> tripletsP;
    Eigen::SparseMatrix<double> objectiveMatrix(objectiveMatrixRowCount, objectiveMatrixColumnCount);
    std::map<Triangle, std::list<DecisionVariableInfo>> triangleToVariableIndexMap{};

    int offset = 0;
    for (const auto &[e, ts] : edgeTriangleMap)
    {
        int i = offset;
        for (Triangle t : ts)
        {
            int j = offset;
            triangleToVariableIndexMap[t].push_back(DecisionVariableInfo{i, t.sourceId, t.indicator});
            for (Triangle t2 : ts)
            {
                if (i == j)
                {
                    tripletsP.emplace_back(i, j, t.variance);
                }
                else
                {
                    tripletsP.emplace_back(i, j, std::sqrt(t.variance) * std::sqrt(t2.variance));
                }
                ++j;
            }
            ++i;
        }

        offset += ts.size();
    }

    // for (auto [t, varList] : triangleToVariableIndexMap)
    // {
    //     std::cout << t << " is represented by indizes";
    //     for (DecisionVariableInfo var : varList)
    //     {
    //         std::cout << " " << var.matrixIndex << " - " << var.countingNodeId << " ";
    //     }
    //     std::cout << std::endl;
    // }

    objectiveMatrix.setFromTriplets(tripletsP.begin(), tripletsP.end());
    return ObjectiveMatrixStats{objectiveMatrix, triangleToVariableIndexMap};
}

Eigen::SparseMatrix<double> setupConstraintMatrix(int triangleCount, std::map<Triangle, std::list<DecisionVariableInfo>> triangleToVariableIndexMap)
{
    std::vector<Eigen::Triplet<double>> triplets;

    int constraintMatrixRowCount = 4 * triangleCount;
    int constraintMatrixColumnCount = 3 * triangleCount;

    int i = 0;
    for (auto &[t, infoList] : triangleToVariableIndexMap)
    {
        for (DecisionVariableInfo varInfo : infoList)
        {
            triplets.emplace_back(i, varInfo.matrixIndex, 1);
        }
        ++i;
    }

    for (int i = 0; i < 3 * triangleCount; i++)
    {
        triplets.emplace_back(triangleCount + i, i, 1.0);
    }

    Eigen::SparseMatrix<double> constraintMatrix(constraintMatrixRowCount, constraintMatrixColumnCount);
    constraintMatrix.setFromTriplets(triplets.begin(), triplets.end());

    return constraintMatrix;
}

std::tuple<Eigen::VectorXd, Eigen::VectorXd> setupConstraintBounds(int triangleCount)
{
    int constraintCount = 4 * triangleCount;
    Eigen::VectorXd upperBound = Eigen::VectorXd::Constant(constraintCount, std::numeric_limits<double>::infinity());
    Eigen::VectorXd lowerBound(constraintCount);

    lowerBound.head(triangleCount) = Eigen::VectorXd::Ones(triangleCount);
    lowerBound.tail(constraintCount - triangleCount) = Eigen::VectorXd::Zero(constraintCount - triangleCount);

    return std::tuple<Eigen::VectorXd, Eigen::VectorXd>{lowerBound, upperBound};
}

Eigen::VectorXd runSolver(Eigen::SparseMatrix<double> Q, Eigen::VectorXd p, Eigen::SparseMatrix<double> A, Eigen::VectorXd l, Eigen::VectorXd u)
{
    osqp::OsqpInstance instance;
    instance.objective_matrix = Q;
    instance.objective_vector = p;
    instance.constraint_matrix = A;
    instance.lower_bounds = l;
    instance.upper_bounds = u;

    osqp::OsqpSolver solver;
    osqp::OsqpSettings settings;

    auto status = solver.Init(instance, settings);
    if (!status.ok())
    {
        std::cerr << "OSQP initialization failed\n";
        return Eigen::VectorXd();
    }

    osqp::OsqpExitCode exit_code = solver.Solve();

    if (exit_code == osqp::OsqpExitCode::kOptimal)
    {
        double optimal_objective = solver.objective_value();
        return solver.primal_solution();
    }
    else
    {
        std::cerr << "OSQP did not find an optimal solution\n";
    }

    return Eigen::VectorXd();
}

double extractTriangleCountFromQP(std::list<Triangle> triangles, std::map<Triangle, std::list<DecisionVariableInfo>> triangleToVariableIndexMap, Eigen::VectorXd solution, double eps2, double lambda)
{
    std::map<int, std::list<std::tuple<double, double>>> countingNodeToTriangleMap;

    for (auto [t, varInfos] : triangleToVariableIndexMap)
    {
        for (DecisionVariableInfo varInfo : varInfos)
        {
            //std::cout << varInfo.countingNodeId << " counts " << t << std::endl;
            countingNodeToTriangleMap[varInfo.countingNodeId].push_back(std::tuple<double, double>{varInfo.indicator, solution[varInfo.matrixIndex]});
        }
    }

    double negativeTriangleCount = 0;
    double unnoisyTriangleCount = 0;
    for (auto [key, trianglesToCount] : countingNodeToTriangleMap)
    {

        double localNegativeTriangleCount = 0;
        double indicatorSum = 0;

        for (auto i : trianglesToCount)
        {
            double indicator = std::get<0>(i);
            localNegativeTriangleCount += indicator * std::get<1>(i);
            indicatorSum += indicator;
            
            //std::cout << key << " counts " << std::get<0>(i) << " with indicator " << std::get<1>(i) << std::endl;
        }
        
        // std::cout << "node: " << key << "neg T count" << localNegativeTriangleCount << ", " << trianglesToCount.size() << std::endl;
        negativeTriangleCount += localNegativeTriangleCount + sampleLaplace(0, indicatorSum * ((1 + 2*lambda / std::pow((1 - lambda), 2))) / eps2);
        unnoisyTriangleCount += localNegativeTriangleCount;
    }

    std::cout << "Negative Triangle Count without second Laplace: " << unnoisyTriangleCount << std::endl;

    return negativeTriangleCount;
}

double QPCountNegativeTriangles(Graph g, double eps, double eps2)
{
    double lambda = std::exp(-eps / 2.0);
    TriangleStats tStats = setupTriangleStatsForQP(&g, eps);

    std::list<Triangle> triangles = std::get<0>(tStats);
    int triangleCount = triangles.size() / 3;

    if (triangleCount == 0)
    {
        return 0;
    }

    std::cout << "Finished Triangle Counting: " << triangleCount << std::endl;

    std::map<std::pair<int, int>, std::list<Triangle>> edgeTriangleMap = std::get<1>(tStats);

    ObjectiveMatrixStats objStats = setupObjectiveMatrix(triangles, edgeTriangleMap);
    std::cout << "Finished: Setup Objective Matrix" << std::endl;

    Eigen::SparseMatrix<double> objectiveMatrix = std::get<0>(objStats);
    std::map<Triangle, std::list<DecisionVariableInfo>> triangleToDecisionVariableMap = std::get<1>(objStats);

    // printMatrix(objectiveMatrix);

    Eigen::SparseMatrix<double> constraintMatrix = setupConstraintMatrix(triangleCount, triangleToDecisionVariableMap);
    std::cout << "Finished: Setup Constraint Matrix" << std::endl;

    std::tuple<Eigen::VectorXd, Eigen::VectorXd> constraintBounds = setupConstraintBounds(triangleCount);
    std::cout << "Finished: Setup Constraint Bounds" << std::endl;

    Eigen::VectorXd lowerBound = std::get<0>(constraintBounds);
    Eigen::VectorXd upperBound = std::get<1>(constraintBounds);

    Eigen::VectorXd sol = runSolver(objectiveMatrix, Eigen::VectorXd::Ones(3 * triangleCount), constraintMatrix, lowerBound, upperBound);

    if (sol.size() == 0)
    {
        std::cout << "Solver did not produce a solution." << std::endl;
        return 0;
    }

    // printMatrix(sol);

    double negativeTriangleCount = extractTriangleCountFromQP(triangles, triangleToDecisionVariableMap, sol, eps2, lambda);

    return negativeTriangleCount;
}