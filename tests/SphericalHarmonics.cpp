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
#include "gtest/gtest.h"
#include "support.hpp"
#include "../src/base/SphericalHarmonics.h"
#include "../src/base/AtomisticModel.h"


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
        auto complex = boost::math::spherical_harmonic(lvalue, mvalue, 0,0);
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
        auto complex = boost::math::spherical_harmonic(lvalue, mvalue, theta, phi);
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
        auto complex = boost::math::spherical_harmonic(lvalue, mvalue, theta, phi);
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
