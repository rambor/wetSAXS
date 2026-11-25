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
// Created by Robert Rambo on 19/05/2022.
//

#include <sastools/PDBModel.h>
#include <random>
#include <cmath>
#include <complex>
#include "gtest/gtest.h"
#include "support.hpp"
#include "../src/base/SphericalHarmonics.h"
#include "../src/base/AtomisticModel.h"

/*
 * Independent reference Y_l^m(theta, phi) in double precision, so the tests do
 * not depend on Boost (macOS libc++ has no std::sph_legendre).
 * Condon-Shortley phase included; Y_l^{-m} = (-1)^m conj(Y_l^m).
 */
static std::complex<double> referenceYlm(int l, int m, double theta, double phi){
    const int am = std::abs(m);
    if (am > l) return {0.0, 0.0};
    const double x = std::cos(theta);
    const double sx = std::sqrt(std::max(0.0, 1.0 - x*x));

    // P_m^m = (-1)^m (2m-1)!! (1-x^2)^{m/2}
    double pmm = 1.0;
    for (int i = 1; i <= am; i++) pmm *= -(2.0*i - 1.0) * sx;
    double plm = pmm;
    if (l > am) {
        double pm1 = x * (2.0*am + 1.0) * pmm;   // P_{m+1}^m
        plm = pm1;
        for (int ll = am + 2; ll <= l; ll++) {
            const double pll = ((2.0*ll - 1.0) * x * pm1 - (ll + am - 1.0) * pmm) / (ll - am);
            pmm = pm1;
            pm1 = pll;
            plm = pll;
        }
    }
    // (l-|m|)!/(l+|m|)! computed as a product to avoid overflow
    double ratio = 1.0;
    for (int k = l - am + 1; k <= l + am; k++) ratio /= k;
    const double norm = std::sqrt((2.0*l + 1.0) / (4.0*M_PI) * ratio);
    std::complex<double> y = norm * plm * std::complex<double>(std::cos(am*phi), std::sin(am*phi));
    if (m < 0) y = ((am % 2) ? -1.0 : 1.0) * std::conj(y);
    return y;
}


class SphericalHarmonicsTests : public ::testing::Test {

public:
    PDBModel bsaModel;
    // setup
    SphericalHarmonicsTests() : ::testing::Test(),
                            bsaModel( PDBModel(fixture(bsa.pdb), false, false) ){
    }
};

TEST_F(SphericalHarmonicsTests, checkYlmSize){
    unsigned int lmax = 5;
    unsigned int numAtoms = 1231;
    unsigned int answer = ((lmax+1)*(lmax+1)*numAtoms);
    SphericalHarmonics she = SphericalHarmonics(lmax, numAtoms);
    EXPECT_EQ(she.getYLMSize(), answer) << "(lmax+1)^2 * numAtoms" << answer;
}


TEST_F(SphericalHarmonicsTests, calculateSHETable){
    unsigned int lmax = 5;
    unsigned int numAtoms = bsaModel.getTotalCoordinates();

    SphericalHarmonics she = SphericalHarmonics(lmax, numAtoms);
    AtomisticModel pdb = AtomisticModel(fixture(bsa.pdb), false, false);

    she.populateSHETable(pdb.getThetas(), pdb.getPhis());

    // given l and m, calcualte for a given coordinate in bsaModel

    // random number based on numAtoms
    std::random_device rd;
    std::mt19937 gen(rd());
    std::vector<unsigned int> indices(numAtoms);
    for(unsigned int i=0; i<numAtoms; i++){
        indices[i] = i;
    }
    std::shuffle(indices.begin(), indices.end(), gen);

    unsigned int lvalue = 0, mvalue = 0;
    for(unsigned int i=0; i<100; i++){
        auto complex = referenceYlm(lvalue, mvalue, 0,0);
        /*
         * m and l at zero should return the normalization constant for real
         * constant = sqrt(1/(4PI))
         * imaginary part is zero due to e^(i*m*phi) where m = 0 implies sin(0) is 0
         */
        EXPECT_NEAR(0.28209479, complex.real(), 0.000001);
    }

    const float * pThetas = pdb.getThetas();
    const float * pPhis = pdb.getPhis();
    float theta, phi;

    /*
     * l = 1, m = 0
     */
    lvalue = 1;
    for(unsigned int i=0; i<100; i++){
        theta = pThetas[indices[i]];
        phi = pPhis[indices[i]];
        auto complex = referenceYlm(lvalue, mvalue, theta, phi);
        /*
         * m = 0 implies sin(0) is 0
         * constant = sqrt(13/(4PI)) = 0.48860251
         *
         */
        float costheta = cosf(theta)*0.48860251;
        EXPECT_NEAR(costheta, complex.real(), 0.000001);
    }

    /*
     * l = 1, m = 1
     */
    mvalue = 1;
    for(unsigned int i=0; i<100; i++){
        theta = pThetas[indices[i]];
        phi = pPhis[indices[i]];
        auto complex = referenceYlm(lvalue, mvalue, theta, phi);
        /*
         * m = 1 implies imaginary term is sin(m*phi)
         * constant = sqrt(3/(4PI)) = 0.48860251
         */
        float costheta = cosf(theta);
        float term = -sqrtf(1.0f - costheta*costheta); // legendre
        float constant = sqrtf((2.0*lvalue+1.0)/(4.0*M_PI)*0.5);

        float real_term = term*constant*cosf(mvalue*phi);
        float imag_term = term*constant*sinf(mvalue*phi);

        EXPECT_NEAR(real_term, complex.real(), 0.000001) << " Real term " ;
        EXPECT_NEAR(imag_term, complex.imag(), 0.000001) << " Imag term " ;
    }
}

/*
 * The library's SHE table against the independent reference, all l <= 12 and
 * -l <= m <= l, on real BSA atom angles. Layout: ((l^2 + l + m) * numAtoms + atom).
 */
TEST_F(SphericalHarmonicsTests, tableMatchesReference){
    const int lmax = 12;
    AtomisticModel pdb = AtomisticModel(fixture(bsa.pdb), false, false);
    const int numAtoms = (int) bsaModel.getTotalCoordinates();

    SphericalHarmonics she = SphericalHarmonics(lmax, numAtoms);
    she.populateSHETable(pdb.getThetas(), pdb.getPhis());
    const float * pR = she.getpDataYLMReal();
    const float * pI = she.getpDataYLMImag();
    const float * pThetas = pdb.getThetas();
    const float * pPhis = pdb.getPhis();

    double maxErr = 0.0;
    for (int l = 0; l <= lmax; l++){
        for (int m = -l; m <= l; m++){
            const size_t base = (size_t)(l*l + l + m) * numAtoms;
            for (int i = 0; i < numAtoms; i += 7){
                auto ref = referenceYlm(l, m, pThetas[i], pPhis[i]);
                maxErr = std::max(maxErr, std::fabs(ref.real() - pR[base + i]));
                maxErr = std::max(maxErr, std::fabs(ref.imag() - pI[base + i]));
            }
        }
    }
    EXPECT_LT(maxErr, 2.0e-6) << "max |Y_lm(table) - Y_lm(reference)| = " << maxErr;
}
