#ifndef ZENO_GF2_VEC_HPP
#define ZENO_GF2_VEC_HPP

#include <vector>
#include <iostream>
#include <cstdint>

namespace zeno {

template<size_t d>
class GF2Vector {
    static_assert(d > 0, "Template parameter d shoud be positive.");
private:
    const int logB = 32;
    const int len  = (d + logB - 1) / logB;
    std::vector<uint32_t> data;
public:


    GF2Vector() {
        data.assign(len, 0);
    }

    template <typename T>
    GF2Vector(std::vector<T> const &vec) {
        data.assign(len, 0);
        /// TODO: ...
    }

    GF2Vector& operator^=(const GF2Vector &rhs) {
        for(int i = 0; i < len; i++)
            this->data[i] ^= rhs->data[i];
        return *this;
    }

    GF2Vector& operator+=(const GF2Vector &rhs) {
        for(int i = 0; i < len; i++)
            this->data[i] ^= rhs->data[i];
        return *this;
    }

    GF2Vector& operator-=(const GF2Vector &rhs) {
        for(int i = 0; i < len; i++)
            this->data[i] ^= rhs->data[i];
        return *this;
    }

    /// @brief Hadamard product.
    GF2Vector& operator*=(const GF2Vector &rhs) {
        for(int i = 0; i < len; i++)
            this->data[i] = this->data[i] & rhs->data[i];
        return *this;
    }

    GF2Vector operator+() const { return GF2Vector(*this); }
    GF2Vector operator-() const { return GF2Vector(*this); }

    GF2Vector operator^(const GF2Vector &rhs) const { return GF2Vector(*this) ^= rhs; }
    GF2Vector operator+(const GF2Vector &rhs) const { return GF2Vector(*this) += rhs; }
    GF2Vector operator-(const GF2Vector &rhs) const { return GF2Vector(*this) -= rhs; }
    GF2Vector operator*(const GF2Vector &rhs) const { return GF2Vector(*this) *= rhs; }

    bool& operator[](size_t i) {
        size_t j = i / logB;
        i %= logB;
        return (data[j] >> i & 1);
    }



    /// @brief Returns a set of GF2 vectors B of minimum cardinality such that span(B) = span(v).
    ///        For other algorithms in GF2Vector to work, B should be such that if the first 1 in
    ///        B[i] is in the j-th column then all the vectors B[i+1...] have 0s in columns [0, j].
    static std::vector<GF2Vector> get_basis(std::vector<GF2Vector> V) {
        /// TODO: implement it with const in V
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
    static bool in_span(std::vector<GF2Vector> U, GF2Vector v) {
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

};



} // namespace zeno


#endif /* ZENO_GF2_VEC_HPP */