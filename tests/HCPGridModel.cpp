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
// Created by Robert Rambo on 26/02/2025.
//
#include <boost/numeric/ublas/matrix.hpp>
#include "gtest/gtest.h"
#include "support.hpp"
#include "HCPGridModel.h"

class HCPTests : public ::testing::Test {

public:

};

TEST_F(HCPTests, getRadiusTest) {
    HCPGridModel hcp(2);
    ASSERT_EQ(hcp.getRadius(),2);
}

TEST_F(HCPTests, getInvRadiusTest) {
    HCPGridModel hcp(2);
    ASSERT_NEAR(hcp.getInvRadius(), 0.5, 0.000001);
}

TEST_F(HCPTests, getInv3RadiusTest) {
    HCPGridModel hcp(2);
    ASSERT_NEAR(hcp.getInvRSqrt3(), 0.2886751294, 0.000001);
}