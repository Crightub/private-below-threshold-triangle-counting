#include "global_sensitivity.h"

#include "basic_counting.h"
#include "distribution.h"
#include "utils.h"
struct PrivateCountingConfig;

double global_sensitivity(const std::vector<std::vector<bool> > &A, int i) {
    int sens = 0;
    // #pragma omp parallel for
    for (int j = i + 1; j < A.size(); ++j) {
        if (!A[i][j])
            continue;
        int sens_j = 0;
        for (int k = i+1; k < A.size(); ++k) {
            if (k == j)
                continue;

            if (A[j][k] && A[k][i]) {
                sens_j++;
            }
        }
        sens = std::max(sens, sens_j);
    }

    return sens;
}

void apply_global_sensitivity(const std::vector<std::vector<bool> > &A,
                              const PrivateCountingConfig &cfg,
                              std::vector<double> &counts) {
    const double sens_mult = cfg.use_unbiased_estimator
                                 ? 1 + 2 * (std::exp(cfg.weight_eps) / std::pow(1 - std::exp(cfg.weight_eps), 2))
                                 : 1;

    // BOOST_LOG_TRIVIAL(info) << "Global sensitivity multiplier: " << sens_mult << " for use_unbiased_estimator: " << cfg.use_unbiased_estimator;

    for (int i = 0; i < counts.size(); ++i) {
        double sens = global_sensitivity(A, i);
        // BOOST_LOG_TRIVIAL(debug) << "Global sensitivity: " << sens << ", laplace b: " << sens * sens_mult / cfg.count_eps;
        counts[i] += sample_laplace(0, sens * sens_mult / cfg.count_eps);
    }
}
