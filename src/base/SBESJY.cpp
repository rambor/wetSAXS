// Copyright (c) 2023.
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
// Created by Robert Rambo on 05/06/2023.
//

#include <string>
#include <cfloat>
#include <stdexcept>
#include <iostream>
#include "SBESJY.h"

// L-T algorithm of continued fractions

SBESJY::SBESJY(float x, int lmax) : qr(x), lmax(lmax) {

    if (lmax < 0)
        throw std::invalid_argument("Invalid argument, LMAX < 0");

    // exception for lmax < 0
    jarray.resize(lmax+1);
    jparray.resize(lmax+1);

    const double ACCUR = 1e-14;

    double invqr = 1.0/qr;
    double j_o = std::sin(qr)*invqr;
    double j_1 = (j_o - std::cos(qr))*invqr;

    if (lmax == 0){
        jarray[0] = j_o;
        return;
    }

    if (lmax == 1){
        jarray[0] = j_o;
        jarray[1] = j_1;
        return;
    }

    // assuming lmax > 2
    jarray[0] = j_o;
    jarray[1] = j_1;


    if (lmax < qr){ // if true, than all l's < lmax are less than qr
        // Miller's algorithm
        double scale, factor, value;
        double invqr2 = invqr*invqr;
        double j_2 = (3.0*invqr2 - 1.0)*j_o - 3.0*std::cos(qr)*invqr2;

        //double j_3 = (15*invqr*invqr*invqr - 6.0*invqr)*j_o - (15*invqr*invqr - 1.0)*std::cos(qr)*invqr;
        jarray[2] = j_2;

        double prev = j_1;
        double current = j_2;

        for (int l=2; l<lmax; l++){

            scale = 1.0;
            factor = (double)(2*l+1)*invqr;
            // do underflow check
            if ((factor > 1) && ( (std::numeric_limits<double>::max() - fabs(prev))/factor < fabs(current) )){
                scale /= current;
                prev /= current;
                current = 1;
            }

            value = factor * current - prev;
            prev = current;
            current = value;

            jarray[l+1] = value/scale;
        }

    } else {

        double twoinvx = 2.0 * invqr;
        double sl = (double)lmax * invqr;
        double tk =  2.0 * sl  + invqr * 3.0;

        double cf1 = sl;  // initialize value
        double den = 1.0; // unnormalised j(Lmax,x)

        double c = cf1;                   // inverse ratio of A convergents
        double d = 0.0;                  // direct  ratio of B convergents

        double dcf1;
        const double small = 1.e-150;

        //
        limit = 20001;
        for(int l=1; l<=limit; l++){
            c = tk - 1.0/c;
            d = tk - d;

            if (std::fabs(c) < small)
                c = small;

            if (std::fabs(d) < small)
                d = small;

            d = 1.0/d;

            dcf1 = c*d;
            cf1 *= dcf1;

            if (d < 0.0)
                den = -den;

            if (fabs(dcf1 - 1.0) <= ACCUR ){
                break;
            }

            tk += twoinvx;
            nfp = l; // number in loop
        }

        if (nfp == limit) { // no convergence
            // throw exception
            throw std::invalid_argument("Invalid argument, LIMIT EXCEEDED");
        }

        double error = ACCUR * sqrt(double(nfp)); // error estimate

        jarray[lmax] = den;
        jparray[lmax] = cf1 * den;

        double * pValue;
        // !------ DOWNWARD RECURSION TO L=0  AS SPHERICAL BESSEL FUNCTIONS
        for (int l=lmax; l > 0; l--){
            pValue = & jarray[l];

            jarray[l-1] = (sl + invqr) * *pValue + jparray[l];

            sl = sl - invqr;
            jparray[l-1] = sl * jarray[l-1] - *pValue;
        }


        den = jarray[0];
        // CALCULATE THE L=0 SPHERICAL BESSEL FUNCTIONS DIRECTLY

        double * pJArray = jarray.data();
        *pJArray =  invqr * sin(qr);

        double omega = *pJArray/den;

        double * pJPArray = jparray.data();
        double y_o = -invqr * cos(qr);
        *pJPArray = -y_o - invqr * *pJArray;

        for (int l=1; l<=lmax; l++){
            ++pJArray;
            ++pJPArray;
            *pJArray *= omega;
            *pJPArray *= omega;
        }






    } // end of the continued fraction


}