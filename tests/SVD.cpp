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
#include <sastools/PDBModel.h>
#include "gtest/gtest.h"
#include "support.hpp"
#include "SVD.h"

class SVDTests : public ::testing::Test {

public:

};

TEST_F(SVDTests, determinantTest) {
    boost::numeric::ublas::matrix<float> a_matrix(3,3);
    a_matrix(0,0) = 2;
    a_matrix(0,1) = -3;
    a_matrix(0,2) = 1;
    a_matrix(1,0) = 2;
    a_matrix(1,1) = 0;
    a_matrix(1,2) = -1;
    a_matrix(2,0) = 1;
    a_matrix(2,1) = 4;
    a_matrix(2,2) = 5;
    SVD svd;
    float determinant = svd.determinant(a_matrix);
    ASSERT_EQ(49, determinant);
}

TEST_F(SVDTests, signTest) {

    SVD svd;
    ASSERT_EQ(-1, svd.sgn(-1));
    ASSERT_EQ(1, svd.sgn(1));
}