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
#include <random>
#include "gtest/gtest.h"
#include "support.hpp"
#include "../src/base/SphericalBessels.h"
#include "../src/base/AtomisticModel.h"

class SphericalBesselsTests : public ::testing::Test {

public:
    PDBModel bsaModel;
    // setup
    SphericalBesselsTests() : ::testing::Test(),
                                bsaModel( PDBModel(fixture(bsa.pdb), false, false) ){
    }
};

TEST_F(SphericalBesselsTests, initializationTest){
    //SphericalBessels(int lmax, int qvaluesSize, int numAtoms, std::vector < float > & qvalues, AtomisticModel &model)
    int lmax = 17;
    std::vector <float> qvalues;

    qvalues.push_back(0.001);
    qvalues.push_back(0.007);
    qvalues.push_back(0.01);
    qvalues.push_back(0.331);
    qvalues.push_back(0.71);

    AtomisticModel pdb = AtomisticModel(fixture(bsa.pdb), false, false);

    SphericalBessels sbj = SphericalBessels(lmax, qvalues.size(), pdb.getTotalAtoms(), qvalues, pdb.getRValuesVector());

    EXPECT_EQ(421380, sbj.getSize());
}