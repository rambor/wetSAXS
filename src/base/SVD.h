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
// Created by Robert Rambo on 27/02/2025.
//

#ifndef WETSAXS_SVD_H
#define WETSAXS_SVD_H

#include <string>
#include <iostream>
#include <boost/filesystem.hpp>
#include <utility>
#include <fstream>
#include <boost/regex.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <random>
#include <sastools/simde-no-tests-master/simde-math.h>

class SVD{

private:
    const float EPS = 0.00001;
    const int ITER_MAX = 50;

public:
    SVD() {}
    /*
     * adapted from SVD-ublas by Volodymyr Kysenko (vksnk)
     */
    float pythag(float a, float b) {
        float absa = fabs(a);
        float absb = fabs(b);

        if (absa > absb) {
            return absa * sqrtf(1.0f + simde_math_powf(absb / absa, 2));
        } else {
            return absb * sqrtf(1.0f + simde_math_powf(absa / absb, 2));
        }
    }

    /*
     * adapted from SVD-ublas by Volodymyr Kysenko (vksnk)
     */
    void identity(boost::numeric::ublas::matrix < float >&m){
        for (unsigned int i = 0; i < m.size1(); i++)
            for (unsigned int j = 0; j < m.size2(); j++)
                m(i, j) = (i == j) ? 1.0f : 0.0f;
    }

    /*
     * adapted from SVD-ublas by Volodymyr Kysenko (vksnk)
     */
    float norm(boost::numeric::ublas::vector < float >&x){
        float x_norm = 0.0;

        for (unsigned int i = 0; i < x.size(); i++)
            x_norm += simde_math_powf(x(i),2);

        return sqrtf(x_norm);
    }

    /*
     * adapted from SVD-ublas by Volodymyr Kysenko (vksnk)
     */
    void normalize(boost::numeric::ublas::vector < float >&x){
        float invx_norm = 1.0f/norm(x);
        for (unsigned int i = 0; i < x.size(); i++) {
            x(i) *= invx_norm;
        }
    }

    void svd(boost::numeric::ublas::matrix<float> &A,
             boost::numeric::ublas::matrix<float> &QQL,
             boost::numeric::ublas::matrix<float> &sigmas,
             boost::numeric::ublas::matrix<float> &QQR) {

        unsigned int row_num = A.size1();
        unsigned int col_num = A.size2();

        QQL.resize(row_num, row_num);
        sigmas.resize(row_num, col_num);
        QQR.resize(col_num, col_num);

        identity(QQL);
        identity(QQR);

        unsigned int to = std::min(row_num, col_num);

        for (unsigned int i = 0; i < to; i++) {
            householder(A, QQL, i, i, true);
            householder(A, QQR, i, i + 1, false);
        }

        boost::numeric::ublas::vector < float > d(to, 0.0f);
        boost::numeric::ublas::vector < float > s(to + 1, 0.0f);

        for (unsigned int i = 0; i < to; i++) {
            d(i) = A(i, i);
            if (i < (to - 1))
                s(i + 1) = A(i, i + 1);
        }

        svd_qr_shift(QQL, QQR, d, s);

        sigmas.clear();

        for (unsigned int i = 0; i < to; i++)
            sigmas(i, i) = d(i); // singular values
    }


    void householder(boost::numeric::ublas::matrix<float> &A,
                     boost::numeric::ublas::matrix<float> &QQ,
                     unsigned int row_start,
                     unsigned int col_start,
                     bool column) {

        unsigned int size = column ? A.size1() : A.size2();
        unsigned int start = column ? row_start : col_start;

        if (start >= size)
            return;

        boost::numeric::ublas::vector < float > x(size);
        for (unsigned int i = 0; i < size; i++) {
            if (i < start) {
                x(i) = 0;
            } else {
                if (column)
                    x(i) = A(i, col_start);
                else
                    x(i) = A(row_start, i);
            }
        }

        float x_norm = norm(x);

        float alpha = (x(start) >= 0.0f ? 1.0f : -1.0f) * x_norm;

        boost::numeric::ublas::vector<float> v = x;

        v(start) += alpha;
        normalize(v);

        boost::numeric::ublas::matrix < float > Q;

        if (column) {
            for (unsigned int i = 0; i < A.size2(); i++) {
                float sum_Av = 0.0f;

                for (unsigned int j = 0; j < A.size1(); j++)
                    sum_Av = sum_Av + (v(j) * A(j, i));
                for (unsigned int j = 0; j < A.size1(); j++)
                    A(j, i) = A(j, i) - 2 * v(j) * sum_Av;
            }

            for (unsigned int i = 0; i < A.size1(); i++) {
                float sum_Qv = 0.0f;

                for (unsigned int j = 0; j < A.size1(); j++)
                    sum_Qv = sum_Qv + (v(j) * QQ(i, j));
                for (unsigned int j = 0; j < A.size1(); j++)
                    QQ(i, j) = QQ(i, j) - 2 * v(j) * sum_Qv;
            }

        } else {
            for (unsigned int i = 0; i < A.size1(); i++) {
                float sum_Av = 0.0f;

                for (unsigned int j = 0; j < A.size2(); j++)
                    sum_Av = sum_Av + (v(j) * A(i, j));
                for (unsigned int j = 0; j < A.size2(); j++)
                    A(i, j) = A(i, j) - 2 * v(j) * sum_Av;
            }

            for (unsigned int i = 0; i < A.size2(); i++) {
                float sum_Qv = 0.0f;

                for (unsigned int j = 0; j < A.size2(); j++)
                    sum_Qv = sum_Qv + (v(j) * QQ(i, j));
                for (unsigned int j = 0; j < A.size2(); j++)
                    QQ(i, j) = QQ(i, j) - 2 * v(j) * sum_Qv;
            }
        }
    }

    void svd_qr_shift(boost::numeric::ublas::matrix<float> &u,
                      boost::numeric::ublas::matrix<float> &v,
                      boost::numeric::ublas::vector<float> &q,
                      boost::numeric::ublas::vector<float> &e) {

        int n = (int)q.size();
        int m = (int)u.size1();

        bool goto_test_conv = false;

        for (int k = n - 1; k >= 0; k--) {
            //std::cout << "U = " << u << std::endl;

            for (int iter = 0; iter < ITER_MAX; iter++) {
                // test for split
                int l;
                for (l = k; k >= 0; l--) {
                    goto_test_conv = false;
                    if (fabs(e[l]) <= EPS) {
                        // set it
                        goto_test_conv = true;
                        break;
                    }

                    if (fabs(q[l - 1]) <= EPS) {
                        // goto
                        break;
                    }
                }

                if (!goto_test_conv) {
                    float c = 0.0;
                    float s = 1.0;

                    int l1 = l - 1;

                    for (int i = l; i <= k; i++) {
                        float f = s * e[i];
                        e[i] = c * e[i];

                        if (fabs(f) <= EPS) {
                            break;
                        }

                        float g = q[i];
                        float h = pythag(f, g);
                        q[i] = h;
                        c = g / h;
                        s = -f / h;

                        for (int j = 0; j < m; j++) {
                            float y = u(j, l1);
                            float z = u(j, i);
                            u(j, l1) = y * c + z * s;
                            u(j, i) = -y * s + z * c;
                        }
                    }
                }

                float z = q[k];

                if (l == k) {
                    if (z < 0.0f) {
                        q[k] = -z;

                        for(int j = 0; j < n; j++)
                            v(j, k) = -v(j,k);
                    }

                    break;
                }

                if (iter >= ITER_MAX - 1) {
                    break;
                }

                float x = q[l];
                float y = q[k - 1];
                float g = e[k - 1];
                float h = e[k];
                float f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0f * h * y);

                g = pythag(f, 1.0);

                if (f < 0) {
                    f = ((x - z) * (x + z) + h * (y / (f - g) - h)) / x;
                } else {
                    f = ((x - z) * (x + z) + h * (y / (f + g) - h)) / x;
                }

                float c = 1.0;
                float s = 1.0;

                for (int i = l + 1; i <= k; i++) {
                    g = e[i];
                    y = q[i];
                    h = s * g;
                    g = c * g;
                    float zz = pythag(f, h);
                    e[i - 1] = zz;
                    c = f / zz;
                    s = h / zz;
                    f = x * c + g * s;
                    g = -x * s + g * c;
                    h = y * s;
                    y = y * c;

                    for (int j = 0; j < n; j++) {
                        x = v(j, i - 1);
                        zz = v(j, i);
                        v(j, i - 1) = x * c + zz * s;
                        v(j, i) = -x * s + zz * c;
                    }

                    zz = pythag(f, h);
                    q[i - 1] = zz;
                    c = f / zz;
                    s = h / zz;
                    f = c * g + s * y;
                    x = -s * g + c * y;

                    for (unsigned int j = 0; j < m; j++) {
                        float yy = u(j, i - 1);
                        float zzz = u(j, i);
                        u(j, i - 1) = yy * c + zzz * s;
                        u(j, i) = -yy * s + zzz * c;
                    }
                }
                e[l] = 0.0;
                e[k] = f;
                q[k] = x;
            }
        }
    }


    void pseudo_inverse(unsigned int lmax,
                        boost::numeric::ublas::matrix<float> &A_star,
                        boost::numeric::ublas::matrix<float> &um,
                        boost::numeric::ublas::matrix<float> &sm,
                        boost::numeric::ublas::matrix<float> &vm) {

        boost::numeric::ublas::matrix<float> inv_sm(sm.size1(), sm.size2(), 0.0f);
        // form pseudo inverse up to the lmax limit
        // not sure if I should be doing reduced SVD or regular?
        for(unsigned int i=0; i<lmax; i++){
            inv_sm(i,i) = 1.0f/sm(i,i);
        }

        // A_star = vm*inv_sm_T*um_T
        A_star = boost::numeric::ublas::prod(boost::numeric::ublas::trans(inv_sm), boost::numeric::ublas::trans(um));
        A_star = boost::numeric::ublas::prod(vm, A_star);
    }


    /*
     * only 3x3 matrix
     */
    float determinant(boost::numeric::ublas::matrix<float> &a_matrix) {
        return a_matrix(0,0)*(a_matrix(1,1)*a_matrix(2,2) - a_matrix(1,2)*a_matrix(2,1))
        - a_matrix(0,1)*(a_matrix(1,0)*a_matrix(2,2) - a_matrix(1,2)*a_matrix(2,0))
        + a_matrix(0,2)*(a_matrix(1,0)*a_matrix(2,1) - a_matrix(1,1)*a_matrix(2,0));
    }

    template <typename T> int sgn(T val) {
        return (T(0) < val) - (val < T(0));
    }

};

#endif //WETSAXS_SVD_H
