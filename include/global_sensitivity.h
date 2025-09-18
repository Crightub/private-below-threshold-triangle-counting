#ifndef NEGATIVE_TRIANGLE_COUNTING_GLOBAL_SENSITIVITY_H
#define NEGATIVE_TRIANGLE_COUNTING_GLOBAL_SENSITIVITY_H

#include <vector>

struct PrivateCountingConfig;

double global_sensitivity(const std::vector<std::vector<bool> > &A, int i);
void apply_global_sensitivity(const std::vector<std::vector<bool> > &A,
                              const PrivateCountingConfig &cfg,
                              std::vector<double> &counts);


#endif //NEGATIVE_TRIANGLE_COUNTING_GLOBAL_SENSITIVITY_H
