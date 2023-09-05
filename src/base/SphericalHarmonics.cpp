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
// Created by Robert Rambo on 04/05/2022.
//

#include "SphericalHarmonics.h"

SphericalHarmonics::SphericalHarmonics(unsigned int lmax, unsigned int numAtoms) : lmax(lmax), numAtoms(numAtoms) {
    ylm_size = (lmax+1)*(lmax+1)*numAtoms;
    y_lm_real.resize((unsigned long)ylm_size);
    y_lm_imag.resize((unsigned long)ylm_size);

//    y_lm_real_per.resize(numAtoms);
//    y_lm_imag_per.resize(numAtoms);

//    plmSize = ((lmax+1)*(lmax+2))/2 + 1; // (lmax+1)*(lmax+2)/2, add one to compensate fo
//    plm_values.resize((numAtoms*(plmSize-1)));

}

/*
 * given a set of atoms in the AtomisticModel container, calculate SHE table
 *
 */
void SphericalHarmonics::populateSHETable( const float * pThetas,  const float * pPhis) {

    SASTOOLS_UTILS_H::logger("CALCULATING", "SPHERICAL HARMONICS");
    SASTOOLS_UTILS_H::logger("", std::to_string(ylm_size));

//    double * const ptrPLM = plm_values.data();
//    float phiValue, sinp, cosp;

//    auto totalsincos = (unsigned int) ((lmax+1)*numAtoms);
//    std::vector<float> costerm(totalsincos);
//    std::vector<float> sinterm(totalsincos);
//    float * const pCos = costerm.data();
//    float * const pSin = sinterm.data();
//
//    int count = 0;
//    for(unsigned int i = 0; i<numAtoms; i++){
//        float phi = model.getPhi(i);
//
//        costerm[count] = 1; // 0*phi
//        sinterm[count] = 0; // 0*phi
//        count++;
//
//        for(int m=1; m <= lmax; m++) {
//            phiValue = m*phi;
//            __sincosf(phiValue, &sinp, &cosp);
//            pCos[count] = cosp;
//            pSin[count] = sinp;
//            count++;
//        }
//    }

    float * const ptrR = y_lm_real.data();
    float * const ptrI = y_lm_imag.data();

    /*
     * Assemble the spherical harmonics for all the atoms
     * for each l, m calculate P_lm*[cos(m*theta) + i * sin(m*theta)]
     */
    unsigned int index_lm = 0;//, lmaxIndex;

    for(int l=0; l <= lmax; l++) { // for each l value, calculate spherical harmonic for each atom

        for(int m=-l; m <= l; m++) { //assemble spherical harmonics at constant l,m for all atoms
            //int lm_index = (l*(l+1))/2 + abs(m);
            for (int i=0 ; i < numAtoms; i++) {
                // includes negative angles for m < 0
                auto complex = boost::math::spherical_harmonic(l, m, pThetas[i], pPhis[i]);
                ptrR[ index_lm ] = (float) (complex.real());//(float) (cosp * result);
                ptrI[ index_lm ] = (float) (complex.imag());//(float) (sinp * result);
                index_lm++;
            }
        }
    }
    // index_lm should equal (lmax+1)*(lmax+1)*numAtoms
}


/*
 * given a set of atoms in the AtomisticModel container, calculate SHE table
 *
 */
//void SphericalHarmonics::calculateSHE(int lvalue, int mvalue, const float * pThetas,  const float * pPhis) {
//
//    float * const ptrR = y_lm_real_per.data();
//    float * const ptrI = y_lm_imag_per.data();
//
//    std::complex<float> complex;
//    /*
//     * Assemble the spherical harmonics for all the atoms
//     * for each l, m calculate P_lm*[cos(m*theta) + i * sin(m*theta)]
//     */
//    for (int i=0 ; i < numAtoms; i++) {
//        // includes negative angles for m < 0
//        complex = boost::math::spherical_harmonic(lvalue, mvalue, pThetas[i], pPhis[i]);
//
//        ptrR[ i ] = (float) (complex.real()*factor);//(float) (cosp * result);
//        ptrI[ i ] = (float) (complex.imag()*factor);//(float) (sinp * result);
//    }
//}