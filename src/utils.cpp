#include "utils.hpp"

// boost::minstd_rand rng(static_cast<unsigned int>(std::time(0)));
boost::minstd_rand rng(static_cast<unsigned int>(0));

void printMatrix(Eigen::VectorXd A)
{
    for (int i = 0; i < A.size(); ++i)
    {
       std::cout << A[i];
       if (i < A.size()-1){
        std::cout << "\t";
       }
    }
    std::cout << std::endl;
}

void printMatrix(std::vector<std::vector<double>> A)
{
    for (int i = 0; i < A.size(); ++i)
    {
        for (int j = 0; j < A[i].size(); ++j)
        {
            std::cout << A[i][j];
            if (j < A[i].size() - 1){
                std::cout << "\t";
            }
        }
        if (i < A.size() - 1)
        {
            std::cout << std::endl;
        }
    }
}

void printMatrix(Eigen::SparseMatrix<double> A) {
    const int max_size = 15;

    if (A.rows() <= max_size && A.cols() <= max_size) {
        Eigen::MatrixXd denseA = Eigen::MatrixXd(A);
        std::cout << denseA << "\n";
    } else {
        for (int k = 0; k < A.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
                std::cout << "(" << it.row() << ", " << it.col() << ") = " << it.value() << "\n";
            }
        }
    }
}

double unbiased_estimator(int x, double p)
{
    if (x > 0)
    {
        return 1.0;
    }
    if (x == 0)
    {
        1 + p / std::pow((1 - p), 2);
    }
    if (x == -1)
    {
        -p / std::pow((1 - p), 2);
    }
    if (x < 0)
    {
        return 0;
    }

    return 0;
}


int biased_estimator(int w1, int w2, int w3){
    return w1 + w2 + w3 < 0;
}

double approximate_variance(double w1, double w2, double w3, double lambda)
{
    return std::pow(lambda, std::abs(-w1 - w2 - w3) + 1);
}

std::ostream &operator<<(std::ostream &out, const Triangle &t)
{
    out << t.sourceId << " - " << t.leftTargetId << " - " << t.rightTargetId << ": " << t.variance;
    return out;
}