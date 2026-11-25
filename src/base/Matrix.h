// Copyright (c) 2025.
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

//
// Minimal dense row-major matrix and vector, replacing boost::numeric::ublas
// for the subset wetSAXS uses (SVD, Aligner, Ensemble).
//
// Differences from uBLAS worth knowing:
//   * storage is always value-initialised (zeros); uBLAS leaves it undefined
//   * resize(r, c) preserves the overlapping block by default (as uBLAS does)
//     and zero-fills new elements
//   * prod() and trans() return concrete matrices, not expression templates.
//     For each output element prod() sums over k in increasing order, which is
//     the same order uBLAS uses, so results match to the last bit unless the
//     compiler contracts to FMA differently.
//   * operator<< prints the uBLAS format: [r,c]((a,b),(c,d))
//

#ifndef WETSAXS_MATRIX_H
#define WETSAXS_MATRIX_H

#include <algorithm>
#include <cstddef>
#include <ostream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace la {

template <typename T>
class Vector {
public:
    using value_type = T;
    using size_type = std::size_t;

    Vector() = default;
    explicit Vector(size_type n) : v_(n, T(0)) {}
    Vector(size_type n, const T &value) : v_(n, value) {}

    size_type size() const { return v_.size(); }
    void resize(size_type n, bool preserve = true) {
        if (preserve) v_.resize(n, T(0)); else v_.assign(n, T(0));
    }
    void clear() { std::fill(v_.begin(), v_.end(), T(0)); }
    void swap(Vector &o) noexcept { v_.swap(o.v_); }

    T &operator()(size_type i) { return v_[i]; }
    const T &operator()(size_type i) const { return v_[i]; }
    T &operator[](size_type i) { return v_[i]; }
    const T &operator[](size_type i) const { return v_[i]; }

    T *data() { return v_.data(); }
    const T *data() const { return v_.data(); }

private:
    std::vector<T> v_;
};

template <typename T>
class Matrix {
public:
    using value_type = T;
    using size_type = std::size_t;

    Matrix() = default;
    Matrix(size_type rows, size_type cols) : r_(rows), c_(cols), d_(rows * cols, T(0)) {}
    Matrix(size_type rows, size_type cols, const T &value) : r_(rows), c_(cols), d_(rows * cols, value) {}

    size_type size1() const { return r_; }
    size_type size2() const { return c_; }

    T &operator()(size_type i, size_type j) { return d_[i * c_ + j]; }
    const T &operator()(size_type i, size_type j) const { return d_[i * c_ + j]; }

    void resize(size_type rows, size_type cols, bool preserve = true) {
        if (rows == r_ && cols == c_) return;
        std::vector<T> nd(rows * cols, T(0));
        if (preserve) {
            const size_type rr = rows < r_ ? rows : r_;
            const size_type cc = cols < c_ ? cols : c_;
            for (size_type i = 0; i < rr; i++)
                for (size_type j = 0; j < cc; j++)
                    nd[i * cols + j] = d_[i * c_ + j];
        }
        d_.swap(nd);
        r_ = rows;
        c_ = cols;
    }

    void clear() { std::fill(d_.begin(), d_.end(), T(0)); }   // zero, keep size (uBLAS semantics)

    void swap(Matrix &o) noexcept {
        std::swap(r_, o.r_);
        std::swap(c_, o.c_);
        d_.swap(o.d_);
    }

    T *data() { return d_.data(); }               // row-major, size1()*size2()
    const T *data() const { return d_.data(); }

private:
    size_type r_ = 0, c_ = 0;
    std::vector<T> d_;
};

template <typename T>
Matrix<T> trans(const Matrix<T> &a) {
    Matrix<T> t(a.size2(), a.size1());
    for (std::size_t i = 0; i < a.size1(); i++)
        for (std::size_t j = 0; j < a.size2(); j++)
            t(j, i) = a(i, j);
    return t;
}

// C = A * B. i-k-j loop order: cache friendly, and each C(i,j) is still
// accumulated over k = 0..K-1 in order, starting from zero.
template <typename T>
Matrix<T> prod(const Matrix<T> &a, const Matrix<T> &b) {
    if (a.size2() != b.size1())
        throw std::invalid_argument("la::prod: inner dimensions do not match");
    const std::size_t n = a.size1(), kk = a.size2(), m = b.size2();
    Matrix<T> c(n, m);
    const T *pa = a.data();
    const T *pb = b.data();
    T *pc = c.data();
    for (std::size_t i = 0; i < n; i++) {
        T *ci = pc + i * m;
        for (std::size_t k = 0; k < kk; k++) {
            const T aik = pa[i * kk + k];
            const T *bk = pb + k * m;
            for (std::size_t j = 0; j < m; j++)
                ci[j] += aik * bk[j];
        }
    }
    return c;
}

// y = A * x
template <typename T>
Vector<T> prod(const Matrix<T> &a, const Vector<T> &x) {
    if (a.size2() != x.size())
        throw std::invalid_argument("la::prod: matrix-vector dimensions do not match");
    Vector<T> y(a.size1());
    for (std::size_t i = 0; i < a.size1(); i++) {
        T s = T(0);
        for (std::size_t k = 0; k < a.size2(); k++) s += a(i, k) * x(k);
        y(i) = s;
    }
    return y;
}

template <typename T>
std::ostream &operator<<(std::ostream &os, const Matrix<T> &m) {
    os << '[' << m.size1() << ',' << m.size2() << "](";
    for (std::size_t i = 0; i < m.size1(); i++) {
        if (i) os << ',';
        os << '(';
        for (std::size_t j = 0; j < m.size2(); j++) {
            if (j) os << ',';
            os << m(i, j);
        }
        os << ')';
    }
    return os << ')';
}

template <typename T>
std::ostream &operator<<(std::ostream &os, const Vector<T> &v) {
    os << '[' << v.size() << "](";
    for (std::size_t i = 0; i < v.size(); i++) {
        if (i) os << ',';
        os << v(i);
    }
    return os << ')';
}

} // namespace la

#endif // WETSAXS_MATRIX_H
