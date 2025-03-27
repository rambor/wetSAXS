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
#include <sastools/PDBModel.h>
#include <boost/numeric/ublas/matrix.hpp>
#include "gtest/gtest.h"
#include "support.hpp"
#include "Aligner.h"

class AlignerTests : public ::testing::Test {

public:

};

TEST_F(AlignerTests, readFile) {

    boost::filesystem::path full_path(boost::filesystem::current_path());
    std::cout << "Current path is : " << full_path << std::endl;
    Aligner align(fixture(aligner.txt), full_path.string());

    ASSERT_EQ(align.getStartResID(),10);
    ASSERT_EQ(align.getEndResID(),101);
}

TEST_F(AlignerTests, parsePDBLines) {

    std::string full_path = test_dir();
    std::cout << "Current path is : " << full_path << std::endl;
    Aligner align(fixture(aligner.txt), full_path);

    ASSERT_EQ(align.getTotalPDBs(),6);
}

TEST_F(AlignerTests, alignSVDFiles) {

    std::string full_path = test_dir();
    std::cout << "Current path is : " << full_path << std::endl;
    Aligner align(fixture(aligner.txt), full_path);

    align.align();
    //align.median_SVD_alignment();
}

TEST_F(AlignerTests, svdAlignmentTest) {

    std::string full_path = test_dir();
    std::cout << "Current path is : " << full_path << std::endl;
    Aligner align(fixture(aligner.txt), full_path);

    boost::numeric::ublas::matrix<float> rotation;
    boost::numeric::ublas::matrix<float> i_m(3,3);

    i_m(0,0) = 1.0;
    i_m(1,1) = 1.0;
    i_m(2,2) = 1.0;

    boost::numeric::ublas::matrix<float> p_ref(5, 3);
    boost::numeric::ublas::matrix<float> q_tar(3, 5);
//ATOM      5  CA  ASP     1      20.990  -3.760  25.490  1.00  0.00           C
    p_ref(0,0) = 20.990;
    p_ref(0,1) = -3.760;
    p_ref(0,2) = 25.490;
// ATOM     17  CA  THR     2      22.600  -0.340  25.240  1.00  0.00           C
    p_ref(1,0) = 22.60;
    p_ref(1,1) = -0.340;
    p_ref(1,2) = 25.24;
// ATOM     31  CA  HIE     3      23.480   1.400  21.840  1.00  0.00           C
    p_ref(2,0) = 23.480;
    p_ref(2,1) = 1.40;
    p_ref(2,2) = 21.840;
//  ATOM     48  CA  LYS     4      26.760   1.490  19.820  1.00  0.00           C
    p_ref(3,0) = 26.760;
    p_ref(3,1) = 1.490;
    p_ref(3,2) = 19.820;
// ATOM     81  CA  GLU     6      25.100   7.750  20.080  1.00  0.00           C
    p_ref(4,0) = 25.10;
    p_ref(4,1) = 7.750;
    p_ref(4,2) = 20.080;

// ATOM      5  CA  ASP X   1      28.935  -3.377  53.288  1.00  0.00           C
    q_tar(0,0) = 28.9350;
    q_tar(1,0) = -3.337;
    q_tar(2,0) = 53.2880;
// ATOM     17  CA  THR X   2      30.390   0.102  52.922  1.00  0.00           C
    q_tar(0,1) = 30.390;
    q_tar(1,1) = 0.102;
    q_tar(2,1) = 52.922;
// ATOM     31  CA  HIE X   3      33.187   1.022  50.335  1.00  0.00           C
    q_tar(0,2) = 33.187;
    q_tar(1,2) = 1.022;
    q_tar(2,2) = 50.335;
// ATOM     48  CA  LYS X   4      37.017   1.211  50.712  1.00  0.00           C
    q_tar(0,3) = 37.017;
    q_tar(1,3) = 1.211;
    q_tar(2,3) = 50.712;
// ATOM     81  CA  GLU X   6      35.603   6.873  47.892  1.00  0.00           C
    q_tar(0,4) = 35.603;
    q_tar(1,4) = 6.873;
    q_tar(2,4) = 47.892;

    vector3 ref_center(0,0,0);
    vector3 tar_center(0,0,0);

    for(unsigned int j=0; j<5; j++){
        ref_center += vector3(p_ref(j,0), p_ref(j,1), p_ref(j,2));
        tar_center += vector3(q_tar(0,j), q_tar(1,j), q_tar(2,j));
    }

    ref_center /= (float)5;
    tar_center /= (float)5;

    // fill matrix after centering coordiantes
    for(unsigned int j=0; j<5; j++){
        p_ref(j,0) -= ref_center.x;
        p_ref(j,1) -= ref_center.y;
        p_ref(j,2) -= ref_center.z;

        q_tar(0,j) -= tar_center.x;
        q_tar(1,j) -= tar_center.y;
        q_tar(2,j) -= tar_center.z;
    }

    boost::numeric::ublas::matrix<float> a_matrix = boost::numeric::ublas::prod(q_tar, p_ref);

    align.getRotationMatrix(a_matrix, i_m, rotation);

    std::cout << rotation << std::endl;

    for(unsigned int j=0; j<5; j++){

        float * px = &q_tar(0,j);
        float * py = &q_tar(1,j);
        float * pz = &q_tar(2,j);

        float new_x =*px * rotation(0,0) + *py * rotation(0,1) + *pz * rotation(0,2);
        float new_y =*px * rotation(1,0) + *py * rotation(1,1) + *pz * rotation(1,2);
        float new_z =*px * rotation(2,0) + *py * rotation(2,1) + *pz * rotation(2,2);

        ASSERT_NEAR(p_ref(j,0), new_x, 0.025);
        ASSERT_NEAR(p_ref(j,1), new_y, 0.025);
        ASSERT_NEAR(p_ref(j,2), new_z, 0.025);
    }

}