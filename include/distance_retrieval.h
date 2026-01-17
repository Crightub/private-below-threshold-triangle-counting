#ifndef NEGATIVE_TRIANGLE_COUNTING_DISTANCE_RETRIEVAL_H
#define NEGATIVE_TRIANGLE_COUNTING_DISTANCE_RETRIEVAL_H

#include <vector>
#include <utility>

class TargetDistance {
public:
    virtual ~TargetDistance() = default;

    virtual void update_target(int new_t) = 0;

    virtual std::pair<int, int> k_th_closest_distance(int k) = 0;

    virtual long long sum_k_closest_distances(int k) = 0;

    virtual int adjacent_to_target_count() const = 0;

    virtual int on_target_count() const = 0;

    virtual int outside_target_count() const = 0;
};


/**
 * Maintains distance statistics to a single target value.
 *
 * Supports:
 *  - k-th closest distance queries
 *  - sum of k closest distances
 *  - target updates
 *
 * This abstraction is used by the smooth sensitivity computation with negative contribution (Appendix A.3).
 */
class SingleTargetDistance : public TargetDistance {
    const std::vector<int> &p; // sorted
    std::vector<long long> pref;
    int t;
    int j; // j denotes the smallest index such that a[j] >= t

    // [ --- x < t-1, t-1, t, t+1, x > t=1 ---]
    int idx_l0, idx_l1; // [idx_l0, idx_l1) == t-1
    int idx_r0, idx_r1; // [idx_r0, idx_r1) == t+1

    int zero_cnt; // count of distance-0 elements (t-1 and t+1)
    int mid_cnt; // count of x == t

public:
    SingleTargetDistance(const std::vector<int> &p, int t);

    std::pair<int, int> k_th_closest_distance(int k) override;

    long long sum_k_closest_distances(int k) override;

    void update_target(int new_t) override;

    [[nodiscard]] int adjacent_to_target_count() const override {
        return zero_cnt;
    }

    [[nodiscard]] int on_target_count() const override {
        return mid_cnt;
    }

    [[nodiscard]] int outside_target_count() const override {
        return p.size() - zero_cnt - mid_cnt;
    }

    [[nodiscard]] int size() const {
        return p.size();
    }

private:
    [[nodiscard]] int left_size() const;

    [[nodiscard]] int right_size() const;

    void rebuild_indices();
};


/**
 * Maintains distance statistics to a double target value.
 *
 * Supports:
 *  - k-th closest distance queries
 *  - sum of k closest distances
 *  - target updates
 *
 * This abstraction is used by the smooth sensitivity computation with positive contribution (Appendix A.1 & A.2).
 */
class DoubleTargetDistance : public TargetDistance {
    const std::vector<int> &p; // sorted
    int t;
    std::vector<long long> pref;

    // [ --- x < t-1, t-1, t, t+1, x > t=1 ---]
    int idx_l0, idx_l1; // [idx_l0, idx_l1) == t-1
    int idx_r0, idx_r1; // [idx_r0, idx_r1) == t+1

    int zero_cnt; // count of distance-0 elements (t-1 and t+1)
    int mid_cnt; // count of x == t


public:
    DoubleTargetDistance(const std::vector<int> &p, int target);

    void update_target(int new_t);

    std::pair<int, int> k_th_closest_distance(int k);

    long long sum_k_closest_distances(int k);

    int adjacent_to_target_count() const {
        return zero_cnt;
    }

    int on_target_count() const {
        return mid_cnt;
    }

    int outside_target_count() const {
        return p.size() - zero_cnt - mid_cnt;
    }

private:
    void rebuild_indices();
};

#endif //NEGATIVE_TRIANGLE_COUNTING_DISTANCE_RETRIEVAL_H
