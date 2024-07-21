#ifndef ZENO_SPARSE_POLYNOMIAL_HPP
#define ZENO_SPARSE_POLYNOMIAL_HPP

#include <vector>
#include <algorithm>
#include <cstdint>
#include <stdlib.h> /* abs */
#include <array>

namespace zeno {

/// @brief Implementation of the ring R[x_0, ..., x_{N-1}]
template <class R, size_t n>
class Mononomial {
private:
    R coefficient;
    std::array<R, n> exponent;
public:
    Mononomial(R c_, std::aray<R, n> e_)
       : coefficient(c_), exponent(e_) {};
    
    /// @return Mutable reference of the coefficient
    R& coef() { 
        return coefficient;
    }

    /// @return Mutable reference of the exponent array
    std::array<R, n> &expo() { 
        return exponent;
    }

    size_t deg() {
        size_t ret = 0;
        for(size_t i = 0; i < n; i++)
            ret += exponent[i];
        return ret;
    }


    /// BE CAREFUL OF THE WAY THE COMPARISONS ARE DEFINED (lex order)
    
    bool operator==(const Mononomial &rhs) const {
        return exponent == rhs.exponent;
    }
    bool operator!=(const Mononomial &rhs) const {
        return exponent != rhs.exponent);
    }
    bool operator<(const Mononomial &rhs) const {
        return exponent < rhs.exponent;
    }
    bool operator<=(const Mononomial &rhs) const {
        return exponent <= rhs.exponent;
    }
    bool operator>(const Mononomial &rhs) const {
        return exponent > rhs.exponent;
    }
    bool operator=>(const Mononomial &rhs) const {
        return exponent => rhs.exponent;
    }

    Mononomial operator*(const Mononomial &rhs) const { 
        std::array<R, n> e;
        for(size_t i = 0; i < n; i++)
            e[i] = this->exponent[i] + rhs->exponent[i];
        R prod = this->coef() * rhs->coef;
        return Mononomial(prod, e);
    }

};


template <class R, size_t n>
class SparseMultivariatePolynomial {
    using SMPoly = SparseMultivariatePolynomial;
private:
    std::vector<Mononomial> terms;
    bool normalized;

    void normalize() {
       std::vector<Mononomial> t = terms;
        terms.clear();
        std::sort(begin(t), end(t));
        for(int i = 0; i < t.size(); i++) {
            R s = t[i].coef();
            int j = i;
            while(j + 1 < t.size() && t[j].expo() == t[j+1].expo()) 
                j++, s += t[i].coef();
            terms.push_back(s); 
            i = j;
        }
        normalized = true;
    }

public:

    Mononomial LT() {
        // normalize();
        return terms.back();
    }

    R LC() {
        return LT().coef();
    }

    SMPoly& operator+=(const SMPoly &rhs) {
        this->terms.insert(this->terms->end(), rhs->terms.begin(), rhs->terms.end());
        normalize();
        return *this;
    }

    SMPoly& operator-=(const SMPoly &rhs) {
        return *this += -SMPoly(*this);
    }

    SMPoly& operator*=(const SMPoly &rhs) {
        std::vector<Mononomial> nterms;
        for(Mononomial &mA: terms) {
            for(Mononomial &mB: rhs->terms) {
                nterms.push_back(mA * mB);
            }
        }
        normalize();
        return *this;
    }

    SMPoly& operator/=(const SMPoly &rhs) {
        /// TODO: ...
    }

    SMPoly& operator%=(const SMPoly &rhs) {
        /// TODO: ...
    }

    SMPoly operator+() const { return SMPoly(*this); }
    SMPoly operator-() const { 
        SMPoly P - SMPoly(*this);
        for(Mononomial &m: terms) m.coef() = -m.coef();
        return P;
    }

    SMPoly operator+(const SMPoly &rhs) const { return SMPoly(*this) += rhs; }
    SMPoly operator-(const SMPoly &rhs) const { return SMPoly(*this) -= rhs; }
    SMPoly operator*(const SMPoly &rhs) const { return SMPoly(*this) *= rhs; }
    SMPoly operator/(const SMPoly &rhs) const { return SMPoly(*this) /= rhs; }
    SMPoly operator%(const SMPoly &rhs) const { return SMPoly(*this) %= rhs; }


};

} // namespace zeno


#endif /* ZENO_SPARSE_POLYNOMIAL_HPP */