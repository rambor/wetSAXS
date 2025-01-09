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
// Created by Robert Rambo on 06/06/2023.
//


#include <sastools/IofQData.h>
#include "gtest/gtest.h"
#include "support.hpp"
#include "../src/base/SBESJY.h"
#include "AtomisticModel.h"

class SBESJYTests : public ::testing::Test {

};


TEST_F(SBESJYTests, initializationTest){
    //SphericalBessels(int lmax, int qvaluesSize, int numAtoms, std::vector < float > & qvalues, AtomisticModel &model)
    int lmax = 37;

    IofQData iofqdata = IofQData(fixture(BSA_0p7.dat), false);
    iofqdata.extractData();
    auto qvalues = iofqdata.getQvalues();

    AtomisticModel md = AtomisticModel(fixture(bsa.pdb), false, false);

    std::vector<float> * pRvalues = md.getRValuesVector();

    double third = -1.0/3.0;



    for(auto & qval : qvalues){

        for(auto rval : *pRvalues){

            double qr = qval*rval;

            SBESJY sbj = SBESJY(qr, lmax);

            auto sbjvalues = sbj.getBesselValues();

            for(int l=0; l<=lmax; l++){
                double value = std::sph_bessel(l, qr);
                //std::cout << l << " " << qr <<  " value " << sbjvalues[l] << " " << value <<  std::endl;
                ASSERT_NEAR( sbjvalues[l], value, FLT_EPSILON) << " qr : " << qr << " l => " << l << " "  << std::endl;
            }
        }
    }
}