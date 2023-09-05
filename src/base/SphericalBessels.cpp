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
// Created by Robert Rambo on 20/05/2022.
//



#include "SphericalBessels.h"
#include "SBESJY.h"

SphericalBessels::SphericalBessels(){
    //besselalq = nullptr;
}

SphericalBessels::SphericalBessels(int lmax, int qvaluesSize, int numatoms, std::vector < float > & qvalues, std::vector < float > * rvalues) :
    lmax(lmax),
    qvaluesSize(qvaluesSize),
    numAtoms(numatoms) {

    SASTOOLS_UTILS_H::logger("CALCULATING", "SPHERICAL BESSELS");
    bessel_size = qvaluesSize*(lmax+1)*numAtoms;

    std::vector<float> qr(numAtoms*qvaluesSize);
    float * const pQR = qr.data();

    int index = 0;

    for(int q = 0; q < qvaluesSize; q++) {
        float * pq = &qvalues[q];
        for (int i=0 ; i < numAtoms; i++) {
            pQR[ index ] = (float)(*pq * (*rvalues)[i]);
            index++;
        }
    }

    calculate(qr);

    //besselalq = new alignas(16) float[bessel_size];
//    besselalq = (float *)_aligned_malloc(sizeof(float)*bessel_size + 16, 16); // 16 byte aligned
//
//    std::clock_t startTime = std::clock();

//    index = 0;
//
//    for(int q = 0; q < qvaluesSize; q++) {
//        int start_at = q*numAtoms;
//        for(int l=0; l <= lmax; l++) { // for a water, calculate bessel's, y_lm's and w_q (scattering factor)
//            // index = (lmax+1)*numAtoms*q + numAtoms*l;
//            // unrolled ?
//            for (int i=0 ; i < numAtoms; i++) {
//                besselalq[ index ] =  boost::math::sph_bessel(l, pQR[start_at+i]);
//                index++;
//            }
//        }
//    }

//    for(int q = 0; q < qvaluesSize; q++) {
//
//        int start_at = q*numAtoms;
//
//        for (int i=0 ; i < numAtoms; i++) {
//
//            float * pVal = &pQR[start_at + i];
//
//            SBESJY temp = SBESJY(*pVal, lmax);
//
//            double * values = temp.getBesselValues();
//
//            int base = (lmax+1)*numAtoms*q + i;
//
//            for(int l=0; l <= lmax; l++) { // for a water, calculate bessel's, y_lm's and w_q (scattering factor)
//                besselalq[ base + numAtoms*l ] = values[l];
//            }
//        }
//    }
//
//    char buffer [10];
//    int n = sprintf(buffer, "%.3f", ((std::clock() - startTime)/(double) CLOCKS_PER_SEC ));
//    SASTOOLS_UTILS_H::logger("TIME", buffer);
//    SASTOOLS_UTILS_H::logger("TOTAL BESSELS", std::to_string(bessel_size));

}


/**
 * Used by Waters since no PDB model is read in
 * @param lmax
 * @param qvaluesSize
 * @param numAtoms
 * @param qvalues
 * @param coords
 */
SphericalBessels::SphericalBessels(int lmax, int qvaluesSize, int numAtoms, std::vector < float > & qvalues, std::vector<Coords> & coords) : lmax(lmax), qvaluesSize(qvaluesSize), numAtoms(numAtoms) {

    bessel_size = qvaluesSize*(lmax+1)*numAtoms;

    std::vector<float> qr(numAtoms*qvaluesSize);
    float * const pQR = qr.data();
    Coords * pCoord = coords.data();

    unsigned int index = 0;
    // pre-calculate q*r
    for(auto & qval : qvalues){
        for (int i=0 ; i < numAtoms; i++) {
            pQR[ index ] = (qval * pCoord[i].r);
            index++;
        }
    } // n-atoms per q value

    calculate(qr);

}

void SphericalBessels::calculate(std::vector<float> & qr){

    float * const pQR = qr.data();

    //besselalq = new alignas(16) float[bessel_size];
    besselalq.resize(bessel_size);
    //besselalq = (float *)_aligned_malloc(sizeof(float)*bessel_size + 16, 16); // 16 byte aligned

    std::clock_t startTime = std::clock();

    for(int q = 0; q < qvaluesSize; q++) {

        int start_at = q*numAtoms;

        for (int i=0 ; i < numAtoms; i++) {

            float * pVal = &pQR[start_at + i];

            SBESJY temp = SBESJY(*pVal, lmax);

            double * values = temp.getBesselValues();

            int base = (lmax+1)*numAtoms*q + i;

            for(int l=0; l <= lmax; l++) { // for a water, calculate bessel's, y_lm's and w_q (scattering factor)
                besselalq[ base + numAtoms*l ] = values[l];
            }
        }
    }

    char buffer [10];
    int n = sprintf(buffer, "%.3f", ((std::clock() - startTime)/(double) CLOCKS_PER_SEC ));
    SASTOOLS_UTILS_H::logger("TIME", buffer);
    SASTOOLS_UTILS_H::logger("TOTAL BESSELS", std::to_string(bessel_size));
}


void SphericalBessels::recalculate(int totalQvalues, std::vector<float> & qvalues, std::vector < float > * rvalues){

    qvaluesSize = totalQvalues;
    bessel_size = qvaluesSize*(lmax+1)*numAtoms;

    std::vector<float> qr(numAtoms*qvaluesSize);
    float * const pQR = qr.data();

    int index = 0;

    for(int q = 0; q < qvaluesSize; q++) {
        float * pq = &qvalues[q];
        for (int i=0 ; i < numAtoms; i++) {
            pQR[ index ] = (float)(*pq * (*rvalues)[i]);
            index++;
        }
    }


    //besselalq = new alignas(16) float[bessel_size];
    besselalq.resize(bessel_size);
//    _aligned_free(besselalq);
//    besselalq = (float *)_aligned_malloc(sizeof(float)*bessel_size + 16, 16); // 16 byte aligned

    std::clock_t startTime = std::clock();

    for(int q = 0; q < qvaluesSize; q++) {

        int start_at = q*numAtoms;

        for (int i=0 ; i < numAtoms; i++) {

            float * pVal = &pQR[start_at + i];

            SBESJY temp = SBESJY(*pVal, lmax);

            double * values = temp.getBesselValues();

            int base = (lmax+1)*numAtoms*q + i;

            for(int l=0; l <= lmax; l++) { // for a water, calculate bessel's, y_lm's and w_q (scattering factor)
                besselalq[ base + numAtoms*l ] = values[l];
            }
        }
    }

    char buffer [10];
    int n = sprintf(buffer, "%.3f", ((std::clock() - startTime)/(double) CLOCKS_PER_SEC ));
    SASTOOLS_UTILS_H::logger("TIME", buffer);
    SASTOOLS_UTILS_H::logger("TOTAL BESSELS", std::to_string(bessel_size));
}