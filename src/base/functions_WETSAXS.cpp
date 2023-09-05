// Copyright (c) 2022.
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
// Created by Robert Rambo on 05/07/2022.
//



#include "functions_WETSAXS.h"

// assumes q is measured as 4PI*sin(theta)/lambda
float functions_WETSAXS::asf ( int atomicNumber, float q) {
    /*
	 * Taken from Waasmaier, D., and Kirfel, A. "New Analytical Scattering-Factor Functions for Free Atoms and Ions" Acta Cryst. (1995) A51, 416-431
	   The last element of the array contains the four gaussian terms for water taken from: Hajdu, F. Acta Cryst (1972). A28, 250. 	 */
    q = (float)(q/(4.0*M_PI));  //need to convert scattering vector to 4PI*sin(theta)/lambda
    float * atom = functions_WETSAXS::coeffs[atomicNumber];
    float q_squared = q*q;

    return	  *(atom  )*expf(-*(atom+6 )*q_squared)
               + *(atom+1)*expf(-*(atom+7 )*q_squared)
               + *(atom+2)*expf(-*(atom+8 )*q_squared)
               + *(atom+3)*expf(-*(atom+9 )*q_squared)
               + *(atom+4)*expf(-*(atom+10)*q_squared)
               + *(atom+5);
}

// assumes q is measured as 4PI*sin(theta)/lambda
float functions_WETSAXS::asf_at_q_zero ( int atomicNumber) {
    /*
	 * Taken from Waasmaier, D., and Kirfel, A. "New Analytical Scattering-Factor Functions for Free Atoms and Ions" Acta Cryst. (1995) A51, 416-431
	   The last element of the array contains the four gaussian terms for water taken from: Hajdu, F. Acta Cryst (1972). A28, 250. 	 */
    float * atom = functions_WETSAXS::coeffs[atomicNumber];
    return *(atom  ) + *(atom+1) + *(atom+2) + *(atom+3) + *(atom+4) + *(atom+5);
}
