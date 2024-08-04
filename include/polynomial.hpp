#pragma once

#include <vector>
#include <algorithm>
#include <iostream>
#include <cstdint>
#include <stdlib.h> /* abs */
#include <array>

namespace zeno {

/// @brief Implementation of mononomials in the ring R[x_0, ..., x_{n-1}]
template <class R, size_t n>
class MononomialTerm {
private:
    R coefficient;
    std::array<int, n> exponent;

public:
    MononomialTerm(R c_)
       : coefficient(c_), exponent{} {};

    MononomialTerm(R c_, std::initializer_list<R> const &l) : coefficient(c_) {
        assert(l.size() == n);
        std::copy(l.begin(), l.end(), exponent.begin());
    };
    MononomialTerm(R c_, std::array<int, n> const &e_)
       : coefficient(c_), exponent(e_) {};
    
    /// @return Mutable reference of the coefficient
    R& coef() { return coefficient; }
    const R& coef() const { return coefficient; }

    /// @return Mutable reference of the exponent array
    std::array<int, n> &expo() { return exponent; }
    const std::array<int, n> &expo() const { return exponent; }

    size_t deg() const {
        size_t ret = 0;
        for(size_t i = 0; i < n; i++)
            ret += exponent[i];
        return ret;
    }


    /// BE CAREFUL OF THE WAY THE COMPARISONS ARE DEFINED (lex order)
    
    bool operator==(const MononomialTerm &rhs) const {
        return exponent == rhs.exponent;
    }
    bool operator!=(const MononomialTerm &rhs) const {
        return exponent != rhs.exponent;
    }
    bool operator<(const MononomialTerm &rhs) const {
        return exponent < rhs.exponent;
    }
    bool operator<=(const MononomialTerm &rhs) const {
        return exponent <= rhs.exponent;
    }
    bool operator>(const MononomialTerm &rhs) const {
        return exponent > rhs.exponent;
    }
    bool operator>=(const MononomialTerm &rhs) const {
        return exponent >= rhs.exponent;
    }

    MononomialTerm& operator+=(const MononomialTerm &rhs) {
        assert(this->exponent == rhs.exponent);
        this->coef() += rhs.coef();
        return *this;
    }

    MononomialTerm& operator-=(const MononomialTerm &rhs) {
        assert(this->exponent == rhs.exponent);
        this->coef() -= rhs.coef();
        return *this;
    }

    MononomialTerm& operator*=(const MononomialTerm &rhs) {
        std::array<int, n> e;
        for(size_t i = 0; i < n; i++)
            e[i] = this->exponent[i] + rhs->exponent[i];
        *this = MononomialTerm(this->coef() * rhs.coef(), e);
        return *this;
    }

    MononomialTerm& operator/=(const MononomialTerm &rhs) {
        std::array<int, n> e;
        for(size_t i = 0; i < n; i++) {
            e[i] = this->exponent[i] - rhs->exponent[i];
            assert(e[i] >= 0);
        }
        *this = MononomialTerm(this->coef() / rhs.coef(), e);
        return *this;
    }

    friend MononomialTerm operator+(const MononomialTerm &lhs, const MononomialTerm &rhs) { return MononomialTerm(lhs) += rhs; }
    friend MononomialTerm operator-(const MononomialTerm &lhs, const MononomialTerm &rhs) { return MononomialTerm(lhs) -= rhs; }
    friend MononomialTerm operator*(const MononomialTerm &lhs, const MononomialTerm &rhs) { return MononomialTerm(lhs) *= rhs; }
    friend MononomialTerm operator/(const MononomialTerm &lhs, const MononomialTerm &rhs) { return MononomialTerm(lhs) /= rhs; }

    MononomialTerm operator+() const { return MononomialTerm(*this); }
    MononomialTerm operator-() const { 
        MononomialTerm m = MononomialTerm(*this);
        m.coef() = -m.coef();
        return m;
    }

    bool divides(MononomialTerm const &m) {
        for(size_t i = 0; i < n; i++)
            if(m.exponent[i] < this->exponent[i]) 
                return false;
        return true;
    }
    

};


/// @brief Sparse implementation of R[x_0, ..., x_{n-1}].
template <class R, size_t n>
class SparseMultivariatePolynomial {
    using SMPoly = SparseMultivariatePolynomial;
    using Mononomial = MononomialTerm<R, n>;
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
                j++, s += t[j].coef();
            terms.emplace_back(s, t[i].expo()); 
            i = j;
        }
        normalized = true;
    }

public:

    SparseMultivariatePolynomial() : terms({}) {}
    SparseMultivariatePolynomial(Mononomial const &m) : terms({m}) {}

    std::vector<Mononomial>& get_terms() {
        return terms;
    }

    const std::vector<Mononomial>& get_terms() const {
        return terms;
    }

    /// @return The leading term.
    Mononomial LT() {
        // normalize();
        return terms.back();
    }

    /// @return The leading mononomial with coefficient R(1).
    Mononomial LM() {
        Mononomial m = LT();
        m.coef() = R(1);
        return m;
    }

    /// @return The leading coefficient.
    R LC() {
        return LT().coef();
    }

    /// @return The multidegree of *this.
    int deg() {
        int ret = -1;
        for(Mononomial &m: terms)
            ret = std::max(ret, m.deg());
        return ret;
    }

    SMPoly& operator+=(const SMPoly &rhs) {
        this->terms.insert(this->terms.end(), rhs.terms.begin(), rhs.terms.end());
        normalize();
        return *this;
    }

    SMPoly& operator-=(const SMPoly &rhs) {
        return *this += -SMPoly(*this);
    }

    SMPoly& operator*=(const SMPoly &rhs) {
        std::vector<Mononomial> nterms;
        for(Mononomial &mA: terms)
            for(Mononomial &mB: rhs.terms)
                nterms.push_back(mA * mB);
        *this = SMPoly(nterms);
        normalize();
        return *this;
    }

    SMPoly& operator/=(const SMPoly &rhs) {
        return division({rhs}).first[0];
    }

    SMPoly& operator%=(const SMPoly &rhs) {
        return division({rhs}).second;
    }

    SMPoly operator+() const { return SMPoly(*this); }
    SMPoly operator-() const { 
        SMPoly P = SMPoly(*this);
        for(Mononomial &m: terms) m.coef() = -m.coef();
        return P;
    }

    SMPoly operator+(const SMPoly &rhs) const { return SMPoly(*this) += rhs; }
    SMPoly operator-(const SMPoly &rhs) const { return SMPoly(*this) -= rhs; }
    SMPoly operator*(const SMPoly &rhs) const { return SMPoly(*this) *= rhs; }
    SMPoly operator/(const SMPoly &rhs) const { return SMPoly(*this) /= rhs; }
    SMPoly operator%(const SMPoly &rhs) const { return SMPoly(*this) %= rhs; }

    std::pair<std::vector<SMPoly>, SMPoly> division(std::vector<SMPoly> f) {
        int s = f.size();
        std::vector<SMPoly> q(s, 0);
        SMPoly r(0), p = SMPoly(*this);
        while(!p.is_zero()) {
            bool div_occured = false;
            for(int i = 0; i < s && div_occured == false; i++) {
                if(f[i].LT().divides(p.LT())) {
                    q[i] += p.LT() / f[i].LT();
                    p -= (p.LT() / f[i].LT()) * f[i];
                    div_occured = true;
                }
            }
            if(!div_occured) {
                r += p.LT();
                p -= p.LT();
            }
        }
        for(SMPoly &e: q) e.normalize();
        r.normalize();
        return {q, r};
    }

};

template <class R, size_t n>
std::ostream& operator<<(std::ostream& os, MononomialTerm<R, n> const &m) {
    os << m.coef() << "x^(";
    for(size_t i = 0; i < n; i++) 
        os << m.expo()[i] << ", "[i + 1 >= n];
    os << ")";
    return os;
}

template <class R, size_t n>
std::ostream& operator<<(std::ostream& os, SparseMultivariatePolynomial<R, n> const &P) {
    os << "{ ";
    for(MononomialTerm<R, n> const &m: P.get_terms()) {
        os << m << " ";
    }
    os << "}";
    return os;
}



} // namespace zeno