#pragma once
#include <sdsl/suffix_trees.hpp>

namespace noLZSS {

using cst_t = sdsl::cst_sada<>;

/**
 * @brief Computes the longest common prefix between two suffixes.
 *
 * Uses the suffix tree's LCA (Lowest Common Ancestor) to efficiently
 * compute the length of the longest common prefix between suffixes
 * starting at positions i and j.
 *
 * @param cst The compressed suffix tree
 * @param i Starting position of first suffix
 * @param j Starting position of second suffix
 * @return Length of the longest common prefix
 */
inline size_t lcp(const cst_t& cst, size_t i, size_t j) {
    if (i == j) return cst.csa.size() - cst.csa[i];
    auto lca = cst.lca(cst.select_leaf(cst.csa.isa[i]+1), cst.select_leaf(cst.csa.isa[j]+1));
    return cst.depth(lca);
}

/**
 * @brief Advances a leaf node by a specified number of positions.
 *
 * Moves from the current leaf node forward by 'iterations' positions
 * in the suffix array order. This is used to advance the current
 * factorization position.
 *
 * @param cst The compressed suffix tree
 * @param lambda Current leaf node
 * @param iterations Number of positions to advance (default: 1)
 * @return The leaf node at the new position
 */
inline cst_t::node_type next_leaf(const cst_t& cst, cst_t::node_type lambda, size_t iterations = 1) {
    auto lambda_rank = cst.lb(lambda);
    for (size_t i = 0; i < iterations; i++) {
        lambda_rank = cst.csa.psi[lambda_rank];
    }
    return cst.select_leaf(lambda_rank + 1);
}

} // namespace noLZSS
