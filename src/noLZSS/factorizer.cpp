#include "factorizer.hpp"
#include <sdsl/suffix_trees.hpp>
#include <sdsl/rmq_support.hpp>
#include <sdsl/bit_vectors.hpp>
#include <cassert>
#include <string>

namespace noLZSS {
using cst_t = sdsl::cst_sada<>;

static size_t lcp(cst_t& cst, size_t i, size_t j) {
    if (i == j) return cst.csa.size() - cst.csa[i];
    auto lca = cst.lca(cst.select_leaf(cst.csa.isa[i]+1), cst.select_leaf(cst.csa.isa[j]+1));
    return cst.depth(lca);
}
static cst_t::node_type next_leaf(cst_t& cst, cst_t::node_type lambda, size_t iterations = 1) {
    assert(cst.is_leaf(lambda));
    auto lambda_rank = cst.lb(lambda);
    for (size_t i = 0; i < iterations; i++) lambda_rank = cst.csa.psi[lambda_rank];
    return cst.select_leaf(lambda_rank + 1);
}
static std::vector<Factor> lzss(cst_t& cst) {
    sdsl::rmq_succinct_sct<> rmq(&cst.csa);
    size_t str_len = cst.size() - 1;
    auto lambda = cst.select_leaf(cst.csa.isa[0] + 1);
    size_t lambda_node_depth = cst.node_depth(lambda);
    size_t lambda_sufnum = 0;
    cst_t::node_type v; size_t v_min_leaf_sufnum = 0; size_t u_min_leaf_sufnum = 0;
    std::vector<Factor> factors; factors.reserve(str_len/4 + 4);
    while (lambda_sufnum < str_len) {
        size_t d = 1; size_t l = 1;
        while (true) {
            v = cst.bp_support.level_anc(lambda, lambda_node_depth - d);
            v_min_leaf_sufnum = cst.csa[rmq(cst.lb(v), cst.rb(v))];
            l = cst.depth(v);
            if (v_min_leaf_sufnum + l - 1 < lambda_sufnum) { u_min_leaf_sufnum = v_min_leaf_sufnum; d++; continue; }
            auto u = cst.parent(v); auto u_depth = cst.depth(u);
            if (v_min_leaf_sufnum == lambda_sufnum) {
                if (u == cst.root()) { l = 1; factors.push_back({(unsigned)lambda_sufnum,(unsigned)l}); break; }
                else { l = u_depth; factors.push_back({(unsigned)lambda_sufnum,(unsigned)l}); break; }
            }
            l = std::min(lcp(cst, lambda_sufnum, v_min_leaf_sufnum), (lambda_sufnum - v_min_leaf_sufnum));
            if (l <= u_depth) { l = u_depth; factors.push_back({(unsigned)lambda_sufnum,(unsigned)l}); break; }
            else { factors.push_back({(unsigned)lambda_sufnum,(unsigned)l}); break; }
        }
        lambda = next_leaf(cst, lambda, factors.back().length);
        lambda_node_depth = cst.node_depth(lambda);
        lambda_sufnum = cst.sn(lambda);
    }
    if (!factors.empty()) factors.pop_back();
    return factors;
}
std::vector<Factor> factorize(std::string_view text) {
    std::string tmp(text);
    if (tmp.empty() || tmp.back() != '$') tmp.push_back('$');
    cst_t cst; construct_im(cst, tmp, 1);
    return lzss(cst);
}
} // namespace noLZSS
