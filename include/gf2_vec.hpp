#ifndef ZENO_GF2_VEC_HPP
#define ZENO_GF2_VEC_HPP

#include <vector>
#include <iostream>
#include <cstdint>
#include <bitset>

namespace zeno {

/// @brief Implementation of the vector space GF(2)^d
template<size_t d>
class GF2Vector {
    static_assert(d > 0, "Template parameter d shoud be positive.");
private:
    // const int logB = 32;
    // const int len  = (d + logB - 1) / logB;
    std::bitset<d> data;
public:


    GF2Vector() {
    }

    template <typename T>
    GF2Vector(std::vector<T> const &vec) {
        /// TODO: ...
    }

    GF2Vector& operator^=(const GF2Vector &rhs) {
        this->data ^= rhs->data;
        return *this;
    }

    GF2Vector& operator+=(const GF2Vector &rhs) {
        this->data ^= rhs->data;
        return *this;
    }

    GF2Vector& operator-=(const GF2Vector &rhs) {
        this->data ^= rhs->data;
        return *this;
    }

    /// @brief Hadamard product.
    GF2Vector& operator*=(const GF2Vector &rhs) {
        this->data = this->data & rhs->data;
        return *this;
    }

    GF2Vector operator+() const { return GF2Vector(*this); }
    GF2Vector operator-() const { return GF2Vector(*this); }

    GF2Vector operator^(const GF2Vector &rhs) const { return GF2Vector(*this) ^= rhs; }
    GF2Vector operator+(const GF2Vector &rhs) const { return GF2Vector(*this) += rhs; }
    GF2Vector operator-(const GF2Vector &rhs) const { return GF2Vector(*this) -= rhs; }
    GF2Vector operator*(const GF2Vector &rhs) const { return GF2Vector(*this) *= rhs; }

    std::bitset::reference& operator[](size_t i) {
        return data[i];
    }


    /// @brief Returns a set of GF2 vectors B of minimum cardinality such that span(B) = span(v).
    ///        For other algorithms in GF2Vector to work, B should be such that if the first 1 in
    ///        B[i] is in the j-th column then all the vectors B[i+1...] have 0s in columns [0, j].
    static std::vector<GF2Vector> get_basis(std::vector<GF2Vector> V) {
        /// TODO: implement it with const in V ???
        std::vector<GF2Vector> B;
        for(size_t j = 0; j < d; j++) {
            if(V[j][j] == 0) continue;
            B.push_back(V[j]);
            for(int i = j + 1; i < V.size(); i++) {
                if(V[i][j] == V[j][j])
                    V[i] -= V[j];
            }
        }
        return B;
    }


    /// @return true iff v belongs in span(U)
    static bool in_span(std::vector<GF2Vector> const &U, GF2Vector v) {
        std::vector<GF2Vector> B = get_basis(U); 
        for(size_t j = 0; j < d; j++) {
            if(v[j] == 0) continue;
            for(int i = 0; i < B.size(); i++) {
                if(B[i][j] == v[j])
                    v -= B[i];
            }
            if(v[j] != 0) 
                return false;
        }
        return true;
    }


    /// @return A basis for the nullspace, ker(A), of A.
    static std::vector<GF2Vector> get_nullspace(std::vector<GF2Vector> const &A) {
        std::vector<GF2Vector> B;
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

        std::vector<GF2Vector> ret;
        for(int f = 0; f < d; f++) if(F[f]) {
            GF2Vector n;
            n[f] = 1;
            for(int i = 0; i < B.size(); i++) {
                n[p[i]] = A[i][f];
            }
            ret.push_back(n);
        }

        return ret;
    }

};



} // namespace zeno


#endif /* ZENO_GF2_VEC_HPP */