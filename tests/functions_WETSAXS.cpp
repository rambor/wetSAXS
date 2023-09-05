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
// Created by Robert Rambo on 08/07/2022.
//

#include <sastools/PDBModel.h>
#include "gtest/gtest.h"
#include "support.hpp"
#include "../src/base/functions_WETSAXS.h"


class FunctionsWETSAXSTests : public ::testing::Test {

public:
    PDBModel bsaModel;
    // setup
    FunctionsWETSAXSTests() : ::testing::Test(),
                            bsaModel( PDBModel(fixture(bsa.pdb), false, false) ){
    }
};

TEST_F(FunctionsWETSAXSTests, asfTest){

    std::map<int, float> asf;

    asf.insert(std::pair<int, float>(1, 0.0f));
    asf.insert(std::pair<int, float>(4, 0.0f));
    asf.insert(std::pair<int, float>(6, 0.0f));
    asf.insert(std::pair<int, float>(7, 0.0f));


    for(auto & ff : asf){ // calculate asf at q of 0.4321f
        ff.second = functions_WETSAXS::asf(ff.first, 0.4321f);
    }

    for(auto & ff : asf){
        EXPECT_LE(ff.second, ff.first);
    }

}