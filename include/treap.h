#pragma once

#include <iostream>
#include <random>

// Treap node definition
struct TreapNode {
    double key;
    int priority; // for balancing
    TreapNode* left;
    TreapNode* right;
    double subtree_sum;
    int subtree_size;

    TreapNode(double k)
        : key(k), left(nullptr), right(nullptr),
          subtree_sum(k), subtree_size(1)
    {
        static std::mt19937 rng(std::random_device{}());
        std::uniform_int_distribution<int> dist(INT_MIN, INT_MAX);
        priority = dist(rng);
    }
};

inline int get_size(TreapNode* t) {
    return t ? t->subtree_size : 0;
}

inline double get_sum(TreapNode* t) {
    return t ? t->subtree_sum : 0.0;
}

inline bool contains(TreapNode* t, double key) {
    if (!t) return false;
    if (t->key == key) return true;
    if (key < t->key)
        return contains(t->left, key);
    else
        return contains(t->right, key);
}

inline void update(TreapNode* t) {
    if (!t) return;
    t->subtree_size = 1 + get_size(t->left) + get_size(t->right);
    t->subtree_sum = t->key + get_sum(t->left) + get_sum(t->right);
}

inline void split(TreapNode* t, double key, TreapNode*& left, TreapNode*& right) {
    if (!t) {
        left = right = nullptr;
    } else if (key <= t->key) {
        split(t->left, key, left, t->left);
        right = t;
        update(t);
    } else {
        split(t->right, key, t->right, right);
        left = t;
        update(t);
    }
}

inline TreapNode* merge(TreapNode* left, TreapNode* right) {
    if (!left || !right) return left ? left : right;
    if (left->priority > right->priority) {
        left->right = merge(left->right, right);
        update(left);
        return left;
    } else {
        right->left = merge(left, right->left);
        update(right);
        return right;
    }
}

inline TreapNode* insert(TreapNode* t, double key) {
    TreapNode* new_node = new TreapNode(key);
    TreapNode *L, *R;
    split(t, key, L, R);
    return merge(merge(L, new_node), R);
}

inline TreapNode* erase(TreapNode* t, double key) {
    if (!t) return nullptr;
    if (t->key == key) {
        TreapNode* res = merge(t->left, t->right);
        delete t;
        return res;
    } else if (key < t->key) {
        t->left = erase(t->left, key);
    } else {
        t->right = erase(t->right, key);
    }
    update(t);
    return t;
}

inline TreapNode* kth(TreapNode* t, int k) {
    if (!t || k <= 0 || k > get_size(t)) return nullptr;
    int left_size = get_size(t->left);
    if (k == left_size + 1) return t;
    if (k <= left_size) return kth(t->left, k);
    return kth(t->right, k - left_size - 1);
}

inline double sum_first_k(TreapNode* t, int k) {
    if (!t || k <= 0) return 0.0;
    int left_size = get_size(t->left);
    if (k <= left_size) return sum_first_k(t->left, k);
    if (k == left_size + 1) return get_sum(t->left) + t->key;
    return get_sum(t->left) + t->key + sum_first_k(t->right, k - left_size - 1);
}

inline void print_inorder(TreapNode* t) {
    if (!t) return;
    print_inorder(t->left);
    std::cout << t->key << " ";
    print_inorder(t->right);
}

inline void print_treap(TreapNode* t, int depth = 0) {
    if (!t) return;

    print_treap(t->right, depth + 1);

    for (int i = 0; i < depth; i++) std::cout << "    ";

    std::cout << "[key=" << t->key
              << ", sum=" << t->subtree_sum
              << ", size=" << t->subtree_size
              << "]" << std::endl;

    print_treap(t->left, depth + 1);
}

inline void print_left_right_treap(TreapNode *left, TreapNode *right) {
    std::cout << "Left Treap: " << std::endl;
    print_treap(left);
    std::cout << "Right Treap: " << std::endl;
    print_treap(right);
}