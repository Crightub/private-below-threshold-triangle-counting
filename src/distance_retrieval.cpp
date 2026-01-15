#include "../include/distance_retrieval.h"

#include <iostream>
#include <vector>
#include "utils.h"

SingleTargetDistance::SingleTargetDistance(std::vector<int> _p, int _t)
    : p(std::move(_p)), t(_t) {
    int n = p.size();
    pref.resize(n + 1, 0);
    for (int i = 0; i < n; ++i)
        pref[i + 1] = pref[i] + p[i];

    j = std::lower_bound(p.begin(), p.end(), t) - p.begin();
    rebuild_indices();
}


long long SingleTargetDistance::sum_k_closest_distances(int k, int l) {
    int r = k - l;
    long long sum = 0;

    // left contribution
    if (l > 0) {
        long long sum_left_p = pref[j] - pref[j - l];
        sum += (long long) l * t - sum_left_p;
    }

    // right contribution
    if (r > 0) {
        long long sum_right_p = pref[j + r] - pref[j];
        sum += sum_right_p - (long long) r * t;
    }

    return sum;
}

int SingleTargetDistance::right_size() const {
    return std::max(static_cast<int>(p.size()) - j, 0);
}

int SingleTargetDistance::left_size() const {
    return std::min(j, static_cast<int>(p.size()));
}


std::pair<int, int> SingleTargetDistance::k_th_closest_distance(int k) {
    int lo = std::max(0, k - right_size());
    int hi = std::min(k, left_size());

    int opt_dis = std::numeric_limits<int>::max();
    int opt_l = 0;

    while (lo <= hi) {
        const int m = std::min((lo + hi) / 2, left_size());
        int left_dist = m > 0 ? t - p[j - m] : -1;
        int right_dist = k - m > 0 ? p[j + (k - m) - 1] - t : -1;
        int n_dis = std::max(left_dist, right_dist);

        if (n_dis <= opt_dis) {
            opt_dis = n_dis;
            opt_l = m;
        }

        if (left_dist == right_dist || hi == lo) {
            return {opt_dis, opt_l};
        }

        if (left_dist < right_dist) {
            // left distance is smaller -> take more element from the right side
            lo = m + 1;
        }

        if (left_dist > right_dist) {
            // right distance is smaller -> take more elements from the left side
            hi = m;
        }
    }

    int l = lo;
    int left_dist = (l > 0) ? t - p[j - l] : -1;
    int right_dist = (k - l > 0) ? p[j + (k - l) - 1] - t : -1;

    if (l == 0) return {right_dist, 0};
    if (k - l == 0) return {left_dist, k};

    return {opt_dis, opt_l};
}

void SingleTargetDistance::update_target(int new_t) {
    t = new_t;

    while (p[j] < t && j < p.size()) {
        ++j;
    }

    rebuild_indices();
}

void SingleTargetDistance::rebuild_indices() {
    idx_l0 = std::lower_bound(p.begin(), p.end(), t - 1) - p.begin();
    idx_l1 = std::upper_bound(p.begin(), p.end(), t - 1) - p.begin();

    idx_r0 = std::lower_bound(p.begin(), p.end(), t + 1) - p.begin();
    idx_r1 = std::upper_bound(p.begin(), p.end(), t + 1) - p.begin();

    zero_cnt = (idx_l1 - idx_l0) + (idx_r1 - idx_r0);

    mid_cnt =
            std::lower_bound(p.begin(), p.end(), t + 1)
            - std::upper_bound(p.begin(), p.end(), t - 1);
}

DoubleTargetDistance::DoubleTargetDistance(const std::vector<int> &points, int target)
    : p(points), t(target) {
    int n = p.size();

    pref.resize(n + 1, 0);
    for (int i = 0; i < n; ++i)
        pref[i + 1] = pref[i] + p[i];

    rebuild_indices();
}

void DoubleTargetDistance::rebuild_indices() {
    idx_l0 = std::lower_bound(p.begin(), p.end(), t - 1) - p.begin();
    idx_l1 = std::upper_bound(p.begin(), p.end(), t - 1) - p.begin();

    idx_r0 = std::lower_bound(p.begin(), p.end(), t + 1) - p.begin();
    idx_r1 = std::upper_bound(p.begin(), p.end(), t + 1) - p.begin();

    zero_cnt = (idx_l1 - idx_l0) + (idx_r1 - idx_r0);

    mid_cnt =
            std::lower_bound(p.begin(), p.end(), t + 1)
            - std::upper_bound(p.begin(), p.end(), t - 1);
}

void DoubleTargetDistance::update_target(int new_t) {
    if (new_t <= t)
        return;

    t = new_t;
    rebuild_indices();
}


std::pair<int, int> DoubleTargetDistance::k_th_closest_distance(int k) {
    // distance 0 elements first
    if (k <= zero_cnt)
        return {0, 0};

    k -= zero_cnt;

    // distance 1 elements (x == t)
    if (k <= mid_cnt)
        return {1, 0};

    k -= mid_cnt;

    // remaining:
    // left side:  x < t-1
    // right side: x > t+1

    int left_cnt = idx_l0;
    int right_start = idx_r1;
    int right_cnt = p.size() - right_start;

    // maximum and minimum number of elements that take be taken from the left side
    int lo = std::max(0, k - right_cnt);
    int hi = std::min(k, left_cnt);

    int opt_dis = std::numeric_limits<int>::max();
    int opt_l = 0;

    while (lo <= hi) {
        int m = (lo + hi) / 2;

        int left_dist = (m > 0)
                            ? (t - 1) - p[idx_l0 - m]
                            : -1;

        int right_dist = (k - m > 0)
                             ? p[right_start + (k - m) - 1] - (t + 1)
                             : -1;

        int n_dis = std::max(left_dist, right_dist);

        std::cout << "m: " << m << " n_dis: " << n_dis << " left_dis = " << left_dist << " right_dis = "
                << right_dist << std::endl;

        if (n_dis <= opt_dis) {
            opt_dis = n_dis;
            opt_l = m;
        }

        if (left_dist == right_dist || lo == hi) {
            return {opt_dis, opt_l};
        }

        if (left_dist < right_dist) {
            // left distance is smaller -> take more element from the left side
            lo = m + 1;
        }

        if (left_dist > right_dist) {
            // right distance is smaller -> take more elements from the right side
            hi = m;
        }
    }

    return {opt_dis, opt_l};
}


long long DoubleTargetDistance::sum_k_closest_distances(int k, int l) {
    long long sum = 0;

    // zero-distance
    int use0 = std::min(k, zero_cnt);
    k -= use0;

    // distance-1 (x == t)
    int use_mid = std::min(k, mid_cnt);
    sum += use_mid;
    k -= use_mid;

    int r = k - l;

    // left contribution
    if (l > 0) {
        long long s = pref[idx_l0] - pref[idx_l0 - l];
        sum += (long long) l * (t - 1) - s;
    }

    // right contribution
    if (r > 0) {
        long long s = pref[idx_r1 + r] - pref[idx_r1];
        sum += s - (long long) r * (t + 1);
    }

    return sum;
}
