#pragma once

#include <vector>
#include <iostream>
#include <cassert>
#include <algorithm> // std::min

#include "convolution.hpp"
#include "internal.hpp"


namespace zeno {

/// @brief Implementation of the R[[x]] ring.
/// @tparam R The coefficient ring
template<typename R>
class FormalPowerSeries {
    using FPS = FormalPowerSeries;

private:
    std::vector<R> a;

public:
    /// @brief Reduce the 0 coefficients of the FPS (i.e remove all trailing zeroes).
    void normalize(){
        while(!a.empty() && a.back() == R(0))
            a.pop_back();
    }

    void resize(size_t n, R v = R(0)) { a.resize(n, v); }
    size_t size()           const { return a.size(); }
    R& at(size_t i)         { return a.at(i); }

    R& operator[](size_t i)             { return a[i]; }
    const R& operator[](size_t i) const { return a[i]; }

    void push_back(R const &v) { a.push_back(v); }
    void pop_back()            { a.pop_back(); }
    R& back()                  { return a.at(a.size()-1); }
    bool empty()         const { return a.size() == 0; }

    /// @brief Constructs a new FPS object with degree -1 (no coefficients).
    FormalPowerSeries(): a({}) {}

    /// @brief Constructs a new FPS object with coefficients from the array coefs.
    FormalPowerSeries(std::vector<R> const &coefs) {
        a = coefs;
        normalize();
    }

    FormalPowerSeries(std::initializer_list<R> const &coefs) {
        a = coefs;
        normalize();
    }

    /// @brief Constructs a new (constant, i.e of 0 degree) FPS object which is just v.
    FormalPowerSeries(const R &v) {
        a.assign(1, v);
        normalize();
    }

    /// @brief Constructss a new FPS object of degree n-1 with initial values v.
    FormalPowerSeries(const int n, const R &v) {
        a.resize(n);
        for(int i = 0; i < n; i++) a.at(i) = v;
        normalize();
    }

    // /// @brief Constructss a new FPS object with coefficients given by the iterators begin, end.
    // FormalPowerSeries(std::vector<R>::iterator begin, std::vector<R>::iterator end) {
    //     this->assign(begin, end);
    //     normalize();
    // }

    /// @return The degree of the FPS (-1 if the FPS is zero)
    int deg() const {
        return this->size() - 1;
    }

    /// @return The order of the FPS (-1 if the FPS is zero)
    int ord() const {
        for(int i = 0; i < this->size(); i++) 
            if(a.at(i) != R(0)) return i;
        return -1;
    }

    int is_zero() const {
        // normalize();
        return this->empty();
    }

    /// @returns A copy of the coefficient of x^i
    R get(const int i) const {
        return (0 <= i && i < a.size() ? a.at(i) : R(0));
    }

    /// @return A copy of the x^i coefficient
    R coef(const int i) const {
        return get(i);
    }

    /// @return Mutable reference of the x^i coefficient. Resizes object if necessary
    R& coef(const int i) {
        assert(i >= 0);
        if(this->size() <= i) this->resize(i + 1, R(0));
        return a.at(i);
    }

    /// @brief Sets the coefficient of x^i
    void set(const int i, const R c) {
        if(c == R(0)) return;
        assert(0 <= i);
        if(this->size() <= i) this->resize(i+1, R(0));
        a.at(i) = c;
    }

    /// @brief Returns a copy of the leading coefficient
    R leadcoef() const {
        if(a.empty()) return 0;
        return a.back();
    }


    /// @brief Addition of two FPSs (*this and P)
    FPS plus(const FPS &P) const {
        return *this + P;
    }

    /// @brief Multiplication of two FPSs (*this and P)
    FPS times(const FPS &P) const {
        return *this * P;
    }

    /// @brief Multiplication of this FPS with the scalar v
    FPS times(const R &v) const {
        return *this * v;
    }

    /// @brief Subtraction of two FPSs (*this - P)
    FPS minus(const FPS &P) const {
        return *this - P;
    }

    FPS& operator+=(const FPS &rhs) {
        int d = (this->deg() > rhs.deg() ? this->deg() : rhs.deg());
        this->resize(d + 1);
        for(int i = 0; i < rhs.size(); i++)
            a[i] += rhs[i];
        this->normalize();
        return *this;
    }

    FPS& operator-=(const FPS &rhs) {
        int d = (this->deg() > rhs.deg() ? this->deg() : rhs.deg());
        this->resize(d + 1);
        for(int i = 0; i < rhs.size(); i++)
            a[i] -= rhs[i];
        this->normalize();
        return *this;
    }

    FPS& operator*=(const FPS &rhs) {
        a = convolution<R>(a, rhs.a);
        this->normalize();
        return *this;
    }

    FPS& operator/=(const FPS &rhs) {
        return *this = div(rhs).first;
    }

    FPS& operator%=(const FPS &rhs) {
        return *this = div(rhs).second;
    }

    FPS& operator+=(const R &rhs) {
        this->coef(0) += rhs;
        return *this;
    }

    FPS& operator-=(const R &rhs) {
        this->coef(0) -= rhs;
        return *this;
    }

    FPS& operator*=(const R &rhs) {
        for(auto &e: this->a) e *= rhs;
        this->normalize();
        return *this;
    }

    FPS& operator/=(const R &rhs) {
        for(auto &e: this) e /= rhs;
        return *this;
    }

    // TODO: >>=, <<= shifts
    // FPS& operator<<=(const int s) {
    //     return *this = (*this).shift(s);
    // }
    // FPS& operator>>=(const int s) {
    //      return operator<<=(-s);
    // }

    FPS operator-() const {
        FPS t = FPS(*this);
        for(auto &e: t) e = -e;
        return t;
    }

    FPS operator+(const FPS &rhs) const { return FPS(*this) += rhs; }
    FPS operator-(const FPS &rhs) const { return FPS(*this) -= rhs; }
    FPS operator*(const FPS &rhs) const { return FPS(*this) *= rhs; }
    FPS operator/(const FPS &rhs) const { return FPS(*this) /= rhs; }
    FPS operator%(const FPS &rhs) const { return FPS(*this) %= rhs; }

    FPS operator+(const R &rhs) const { return FPS(*this) += rhs; }
    FPS operator-(const R &rhs) const { return FPS(*this) -= rhs; }
    FPS operator*(const R &rhs) const { return FPS(*this) *= rhs; }
    FPS operator/(const R &rhs) const { return FPS(*this) /= rhs; }


    static bool _vec_weak_equals(std::vector<R> const &A, std::vector<R> const &B) {
        int i = A.size() - 1, j = B.size() - 1;
        while(i >= 0 && A[i] == R(0)) i--;
        while(j >= 0 && B[j] == R(0)) j--;
        if(i != j) return false;
        while(i >= 0)
            if(A[i] != B[i]) return false;
            else i--;
        return true;
    }

    /// @warning operator perhaps changes lhs, rhs
    friend bool operator==(FPS const &lhs, FPS const &rhs) { 
        return _vec_weak_equals(lhs.a, rhs.a);
    }

    /// @warning operator perhaps changes lhs, rhs
    friend bool operator!=(FPS const &lhs, FPS const &rhs) { 
        return !_vec_weak_equals(lhs.a, rhs.a);
    }

    // bool operator == (const FPS &rhs) const { normalize(); return *this == rhs; }
    // bool operator != (const FPS &rhs) const { normalize(); return *this != rhs; }

    /// @return A copy of *this mod x^k (only the the first k coefficients)
    FPS mod_xk(size_t k) const {
        // std::vector<R> aa(a.begin(), a.begin() + std::min(k, a.size()));
        std::vector<R> aa;
        for(size_t i = 0; i < k && i < a.size(); i++) aa.push_back(a[i]);
        return FPS(aa);
    }
    inline static FPS mod_xk(FPS &f, size_t k) { return f.mod_xk(k); }

    /// @return A copy of *this times x^k (coefficients shifted k steps to the right)
    FPS times_xk(int k) const {
        std::vector<R> coefs = *this;
        coefs.insert(coefs.begin(), k, R(0));
        return FPS(coefs);
    }
    inline static FPS times_xk(FPS &f, int k) { return f.times_xk(k); }

    /// @return A copy of *this with coefficients shifted k steps to the left
    FPS floordiv_xk(int k) const {
        return FPS(std::vector<R>(this->begin() + k, this->end()));
    }
    inline static FPS floordiv_xk(FPS &f, int k) { return f.floordiv_xk(k); }

    /// @return This FPS where the coefficients are shifted by k (k can be both positive/negative).
    FPS shift(int k) const {
        if(k > 0)
            return times_xk(k);
        if(k < 0)
            return floordiv_xk(-k);
        return FPS(*this);
    }
    inline static FPS shift(FPS &f, int k) { return f.shift(k); }


    /// @return P(x^{1/2}) where P = *this.
    FPS sqrt_x() const {
        std::vector<R> ret(this->size(), R(0));
        for(int i = 0; i < this->size(); i++)
            ret[i] = get(2*i);
        return FPS(ret);
    }
    inline static FPS sqrt_x(FPS &f) { return f.sqrt_x(); }

    /// @return P(x^2) where P = *this.
    FPS sqr_x() const {
        std::vector<R> ret(2 * this->size(), R(0));
        for(int i = 0; i < this->size(); i++)
            ret[2*i] = this->a[i];
        return FPS(ret);
    }
    inline static FPS sqr_x(FPS &f) { return f.sqr_x(); }

    /// @return The multiplicative inverse of *this mod x^m.
    FPS inv(int m = -1) const {
        if(m == -1) m = deg() + 1;
        assert(get(0) != R(0));
        if(m == 0) return FPS();
        FPS Q(R(1) / get(0));
        for(int i = 1; i < m; i *= 2)
            Q = (Q * R(2) - Q * Q * this->mod_xk(2*i)).mod_xk(2*i);
        return Q.mod_xk(m);
    }
    inline static FPS inv(FPS &f, int m = deg() + 1) { return f.inv(m); }

    /// @return A copy of *this with the coefficients reversed.
    FPS reverse() const {
        return FPS(std::vector<R>(a.rbegin(), a.rend()));
    }
    inline static FPS reverse(FPS &f) { return f.reverse(); }

    /// @brief Euclidean division of this FPS with the FPS P.
    /// @return pair<FPS, FPS> the quotiend (first) and remainder (second) FPSs.
    std::pair<FPS, FPS> div(const FPS &P) {
        assert(!P.is_zero());
        if(P.deg() > this->deg())
            return {FPS(0), FPS(*this)};
        FPS A_R = this->reverse(), B_R = P.reverse();
        int d = this->deg() - P.deg();
        FPS D = ((A_R).mod_xk(d + 1) * P.inv(d + 1)).mod_xk(d + 1).reverse();
        return {D, *this - D * P};
    }
    inline static std::pair<FPS, FPS> div(FPS &f, const FPS &Q) { return f.div(Q); }

    /// @return A copy of *this with the x value scaled by v [P(x) -> P(vx)].
    FPS scale_x(const R &v) const {
        R c = 1;
        FPS P(*this);
        for(int i = 0; i <= P.deg(); i++, c *= v) { 
            P[i] *= c;
        }
        return  P;
    }
    inline static FPS scale_x(FPS &f, const R &v) { return f.scale_x(v); }



    /// @return The derivative of *this
    FPS deriv(int k = 1) const {
        if(deg() + 1 < k) return FPS(R(0));
        std::vector<R> res(deg() - k +  1, R(0));
        if(k == 1) {
            for(int i = 1; i <= deg(); i++)
                res[i - 1] = a.at(i) * R(i);
        } else {
            for(int i = k; i <= deg(); i++)
                res[i - k] = a.at(i) * internal::fact<R>(i) * internal::inv_fact<R>(i - k);
        }
        return FPS(res);
    }

    inline static FPS deriv(FPS &f, int k = 1) { return f.deriv(k); }


    /// @return The integral of *this (with constant 0)
    FPS integr() const {
        std::vector<R> res(deg() + 2);
        for(int i = 1; i <= deg() + 1; i++) {
            res[i] = a.at(i-1) / i;
        }
        return FPS(res);
    }
    inline static FPS integr(FPS &f) { return f.integr(); }

    /// @return log(*this) (natural logarithm).
    FPS log(int m = -1) const { 
        if(m == -1) m = deg() + 1;
        FPS P(*this);
        return (P.deriv().mod_xk(m) * P.inv(m)).mod_xk(m - 1).integr();
    }
    inline static FPS log(FPS &f, int m) { return f.log(m); }

    /// @return exp(*this).
    FPS exp(int m = -1) const {
        if(m == -1) m = deg() + 1;
        assert(get(0) == 0);
        FPS Q(1), P(*this);
        P[0] = 1;
        for(int i = 1; i < m; i *= 2) { // at the end of each step: F = e^P (mod x^{2^{i+1}})
            Q = (Q * (P - Q.log(2*i))).mod_xk(2*i);
        }
        return Q.mod_xk(m);
    }
    inline static FPS exp(FPS &f, int m) { return f.exp(m); }

    // FPS sqrt(int k = 2) const {

    // }
    // inline static FPS sqrt(FPS &f, int m) { return f.sqrt(m); }

    /// @return The k-th power of the FPS
    FPS pow(int k, int m = -1) const {
        if(m == -1) m = deg() * k;
        if(is_zero()) {
            assert(k > 0);
            return (k == 0 ? FPS(1) : FPS(0));
        }

        if(k < 0)
            return pow(-k, m).inv();

        int ord = this->ord();
        R a = a.at(ord);
        FPS T(this->begin() + ord, this->end());
        T /= a; // P(x) = ax^ord T(x); T(0) = 1
        
        // P(x)^k = a^k x^{k*ord} T(x)^k = a^k x^{k*ord} exp[klnT(x)]
        return (T.log(m) * k).exp(m).shift(k * ord) * internal::pow(a, k);
    }

    /// @return The value of the FPS at point x0 
    R eval(R x0) const {
        R y(0);
        for(int i = this->deg(); i >= 0; i--) {
            y = (x0 * y) + a.at(i);
        }
        return y;
    }

    /// @brief Builds the FPS evaluation tree.
    /// @return The FPS P(x) = (x - x_l) (x - x_{l+1}) ... (x - x_{r-1}).
    static FPS build_poly_tree(std::vector<FPS> &tree, 
        int v, size_t l, size_t r, std::vector<R> const &x) {
            if(r - l == 1) {
                return tree[v] = FPS({-x[l], 1});
            } else {
                size_t m = l + (r - l) / 2;
                return tree[v] = build_poly_tree(tree, 2*v, l, m, x)
                    * build_poly_tree(tree, 2*v+1, m, r, x);
            }
        }

    std::vector<R> _eval(std::vector<FPS> &tree, 
        int v, size_t l, size_t r, std::vector<R> const &x) const {
            if(r - l == 1) {
                return {eval(x[l])};
            } else {
                size_t m = l + (r - l) / 2;
                std::vector<R> A = (*this % tree[2*v])._eval(tree, 2*v, l, m, x);
                std::vector<R> B = (*this % tree[2*v+1])._eval(tree, 2*v+1, m, r, x);
                A.insert(A.end(), B.begin(), B.end());
                return A;
            }
        }

    /// @brief Evaluate the FPS in points x0, x1, ..., x_{n-1} in O(nlog^2 n).
    std::vector<R> eval(std::vector<R> const &x) {
        int n = x.size();
        if(is_zero()) return std::vector<R>(n, R(0));
        std::vector<FPS> tree(4*n);
        build_poly_tree(tree, 1, 0, x.size(), x);
        return _eval(tree, 1, 0, x.size(), x);
    }


    static FPS _interpolate_tree(std::vector<FPS> &Atree,  std::vector<FPS> const &Ptree,
        int v, size_t l, size_t r, std::vector<R> const &x) {
            if(r - l == 1) {
                return Atree[v] = FPS(x[l]);
            } else {
                auto m = l + (r - l) / 2;
                FPS A = _interpolate_tree(Atree, Ptree, 2*v, l, m, x);
                FPS B = _interpolate_tree(Atree, Ptree, 2*v+1, m, r, x);
                return Atree[v] = A * Ptree[2*v+1] + Ptree[2*v] * B;
            }
    }


    /// @brief Polynomial interpolation in O(nlog^2 n)
    /// @pre x All x values should be unique.
    /// @return An interpolated polynomial A(x) s.t. A(x[i]) = y[i] for i = 0, 1, ..., n-1.
    static FPS interpolate(std::vector<R> const &x, std::vector<R> const &y) {
        assert(x.size() == y.size());
        std::vector<FPS> Ptree(4*x.size()), Atree(4*x.size()); 
        // P(x) = (x - x_0) ... (x - x_{n-1}) and A(x) is the interpolation polynomial
        std::vector<R> u = FPS::build_poly_tree(Ptree, 1, 0, x.size(), x).deriv().eval(x);
        for(int i = 0; i < y.size(); i++)
            u[i] = y[i] / u[i];
        return _interpolate_tree(Atree, Ptree, 1, 0, u.size(), u);
    }

    /// @return The FPS P(x) = Prod_{i=0}^n f_i(x).
    static FPS eval_product(std::vector<FPS> &f) {
        std::deque<FPS> q;
        for(FPS &e: f) q.emplace_back(e);
        for(int i = 0; i < f.size() - 1; i++) {
            FPS A = q.front();
            q.pop_front();
            FPS B = q.front();
            q.pop_front();
            q.emplace_back(A * B);
        }
        return q.front();
    }


    /// @return The Borel transformation of the FPS (a_k x^k -> a_k / k! x^k).
    FPS borel() const {
        FPS ret = FPS(*this);
        for(int k = 0; k <= deg(); k++)
            ret[k] *= internal::inv_fact(k);
        return ret;
    }

    /// @return The inverse Borel (Laplace) transformation of the FPS (a_k x^k -> k!a_k x^k).
    FPS inv_borel() const {
        FPS ret = FPS(*this);
        for(int k = 0; k <= deg(); k++)
            ret[k] *= internal::fact(k);
        return ret;
    }


    static std::pair<FPS, FPS> euclid_division(FPS const &A, FPS const &B) {
        return A.div(B);
    }

    static void euclid_division(FPS const &A, FPS const &B, FPS &Qf, FPS &Rf) {
        std::pair<FPS, FPS> p = A.div(B);
        Qf = p.first, Rf = p.second;
    }

    static FPS gcd(FPS const &A, FPS const &B) { // general gcd also works
        return (B.is_zero() ? A : FPS::gcd(B, A % B));
    }

    static FPS extended_gcd(FPS const &A, FPS const &B, FPS &U, FPS &V) {
        FPS D(A), V1 = FPS(0), V3 = FPS(B), T, Q, Rem;
        U = FPS(1);
        while(!V3.is_zero()) {
            FPS::euclid_division(D, V3, Q, Rem);
            T = U - V1 * Q;
            U = V1, D = V3;
            V1 = T, V3 = Rem;
        }
        V = (D - A*U) / B;
    }


    /// @return P(x + c) where P is this FPS.
    FPS taylor_shift(R c) {
        size_t n = this->deg();
        FPS a, b; a.resize(n), b.resize(n);
        R cc = c;
        for(int i = 0; i < n; i++) a[i] = internal::fact(i) * a.at(i);
        for(int i = 0; i < n; i++) b[i] = cc * internal::inv_fact(i), cc *= c;
        return (a * b).borel();
    }

    /// @return The pointwise product of this FPS with f.
    FPS hadamard_product(FPS const &f) {
        FPS ret(0);
        for(int i = 0; i < this->size() && i < f.size(); i++)
            ret.push_back(a.at(i) * f[i]);
        ret.normalize();
        return ret;
    }

};

template<class R>
using FPS = FormalPowerSeries<R>;

template<class R>
using Poly = FormalPowerSeries<R>;

template<class R>
FormalPowerSeries<R>  operator*(const R &lhs, const FormalPowerSeries<R> &rhs) {
    return FPS<R>(rhs) *= lhs;
}

template<class R>
FormalPowerSeries<R> operator/(const R &lhs, const FormalPowerSeries<R> &rhs) {
    return FPS(lhs) /= rhs;
}


template<class R>
std::ostream& operator<<(std::ostream& os, const FormalPowerSeries<R> &P) {
    os << "[ ";
    for(int i = 0; i <= P.deg(); i++)
        os << P.get(i) << " ";
    os << "] ";
    return os;
}

} // namespace zeno