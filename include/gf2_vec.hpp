#pragma once

#include <vector>
#include <iostream>
#include <cstdint>
#include <bitset>

namespace zeno {

/// @brief Implementation of the GF(2)^d vector space
namespace GF2 {
    
template<size_t d>
using GF2Vec = std::bitset<d>;

/// @return A list L of GF2 vectors such that the in L[j] the j-th bit is the leftmost one set, 
///         span(L) = span(V), and L forms a basis (if we exclude the 0 vectors of L).
template<size_t d>
std::vector<GF2Vec<d>> get_fullbasis(std::vector<GF2Vec<d>> const &V) {
    std::vector<GF2Vec<d>> B(d, 0);
    for(GF2Vec<d> v: V) {
        for(size_t j = 0; j < d; ++j) {
            if(v[j] == 0) continue;
            if(!B[j].any()) {
                B[j] = v;
                break;
            }
            v ^= B[j];
        }
    }
    return B;
}

/// @return A basis of V.
template<size_t d>
std::vector<GF2Vec<d>> get_basis(std::vector<GF2Vec<d>> const &V) {
    std::vector<GF2Vec<d>> ret;
    for(GF2Vec<d> &b: get_fullbasis(V)) if(b.any()) ret.push_back(b);
    return ret;
}


/// @return true iff v belongs in span(U)
template<size_t d>
bool in_span(std::vector<GF2Vec<d>> const &U, GF2Vec<d> v) {
    std::vector<GF2Vec<d>> B = get_fullbasis(U); 
    for(size_t j = 0; j < d; j++) {
        if(v[j] == 0) continue;
        v ^= B[j];
    }
    return !v.any();
}

/// TODO: Check correctness ???
/// @return A basis for the intersection of span(U) and span(V).
template<size_t d>
std::vector<GF2Vec<d>> intersect(std::vector<GF2Vec<d>> const &U, std::vector<GF2Vec<d>> const &V) {
    std::vector<GF2Vec<d>> I;
    std::vector<GF2Vec<d>> VB = get_fullbasis(V);
    for(GF2Vec<d> uu: U) {
        GF2Vec<d> u = uu;
        for(size_t j = 0; j < d; ++j) if(u[j] == 1) {
            if(VB[j].any()) u ^= VB[j];
        }
        if(u == 0) I.push_back(uu);
    }
    return get_basis(I);
}

/// TODO: Check correctness
/// @return A basis for the nullspace, ker(A), of A.
template<size_t d>
std::vector<GF2Vec<d>> get_nullspace(std::vector<GF2Vec<d>> const &A) {
    std::vector<GF2Vec<d>> B;
    for(int i = 0; i < A.size(); i++) {
        if(A[i].any())
            B.push_back(A[i]);
    }

    std::vector<bool> F(d, true);
    std::vector<int>  p(d, 0);

    for(int i = 0; i < B.size(); i++) {
        for(int j = 0; j < d; j++) if(B[i][j]) {
            F[j] = false, p[i] = j;
        }
    }

    std::vector<GF2Vec<d>> ret;
    for(int f = 0; f < d; f++) if(F[f]) {
        GF2Vec<d> n;
        n[f] = 1;
        for(int i = 0; i < B.size(); i++) {
            n[p[i]] = A[i][f];
        }
        ret.push_back(n);
    }

    return ret;
}


} // namespace GF2

} // namespace zeno