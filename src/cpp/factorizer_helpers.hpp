#pragma once
#include <sdsl/suffix_trees.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/wt_algorithm.hpp>
#include <algorithm>
#include <cstdint>
#include <tuple>

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

/**
 * @brief Range-predecessor index over an array of text positions.
 *
 * Answers "largest value <= limit among A[lb..rb]" in O(log n) time using a wavelet tree over A.
 * It is the structure behind RefPolicy::Last: the candidates for a factor of length L at position i are
 * the suffix-array interval of the factor's locus node, and the closest non-overlapping one is the largest
 * start <= i - L (forward) or the largest end <= i - 1 (reverse complement) in that interval.
 *
 * Values >= cap are treated as absent (they are stored as cap and never returned), which keeps the tree
 * at ceil(log2(cap + 1)) levels even when the source array uses a large INF marker.
 *
 * Queries are read-only and thread-safe.
 */
class RangePredecessor {
public:
    using wt_type = sdsl::wt_int<sdsl::bit_vector,
                                 sdsl::rank_support_v5<1>,
                                 sdsl::select_support_scan<1>,
                                 sdsl::select_support_scan<0>>;

    RangePredecessor() = default;

    /**
     * @param values Array of positions aligned to the suffix array (one entry per SA rank)
     * @param cap    Smallest value that counts as "absent" (e.g. the text length)
     */
    RangePredecessor(const sdsl::int_vector<64>& values, uint64_t cap) : m_cap(cap) {
        sdsl::int_vector<64> clamped(values.size(), 0);
        for (size_t k = 0; k < values.size(); ++k) {
            uint64_t v = values[k];
            clamped[k] = v < cap ? v : cap;
        }
        m_wt = wt_type(clamped.begin(), clamped.end());
    }

    /**
     * @brief Largest value <= limit in A[lb..rb] (inclusive bounds).
     * @return true and sets out when such a value exists, false otherwise
     */
    bool predecessor(size_t lb, size_t rb, uint64_t limit, uint64_t& out) const {
        if (m_cap == 0 || lb > rb) return false;
        if (limit >= m_cap) limit = m_cap - 1;             // values equal to cap are absent
        auto counts = m_wt.lex_count(lb, rb + 1, limit + 1); // (rank, #smaller than limit+1, #greater)
        size_t cnt = std::get<1>(counts);                  // == number of values <= limit
        if (cnt == 0) return false;
        out = sdsl::quantile_freq(m_wt, lb, rb, cnt - 1).first;
        return true;
    }

    bool empty() const { return m_cap == 0; }

private:
    wt_type  m_wt;
    uint64_t m_cap = 0;
};

} // namespace noLZSS
