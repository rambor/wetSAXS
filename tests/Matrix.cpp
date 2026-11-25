//
// Tests for the internal dense matrix (src/base/Matrix.h) that replaced
// boost::numeric::ublas in SVD, Aligner and Ensemble.
//

#include <cmath>
#include <random>
#include <sstream>
#include "gtest/gtest.h"
#include "Matrix.h"
#include "SVD.h"

static la::Matrix<float> randomMatrix(size_t r, size_t c, unsigned seed){
    std::mt19937 gen(seed);
    std::uniform_real_distribution<float> dist(-1.0f, 1.0f);
    la::Matrix<float> m(r, c);
    for (size_t i = 0; i < r; i++)
        for (size_t j = 0; j < c; j++)
            m(i, j) = dist(gen);
    return m;
}

TEST(MatrixTests, constructZeroAndFill){
    la::Matrix<float> a(3, 4);
    EXPECT_EQ(a.size1(), 3u);
    EXPECT_EQ(a.size2(), 4u);
    for (size_t i = 0; i < 3; i++)
        for (size_t j = 0; j < 4; j++)
            EXPECT_EQ(a(i, j), 0.0f);

    la::Matrix<float> b(2, 2, 1.5f);
    EXPECT_EQ(b(1, 1), 1.5f);
    b.clear();
    EXPECT_EQ(b(1, 1), 0.0f);
    EXPECT_EQ(b.size1(), 2u);
}

TEST(MatrixTests, resizePreservesBlock){
    la::Matrix<float> a = randomMatrix(3, 3, 1);
    la::Matrix<float> keep(a);
    a.resize(5, 2, true);
    EXPECT_EQ(a.size1(), 5u);
    EXPECT_EQ(a.size2(), 2u);
    for (size_t i = 0; i < 3; i++)
        for (size_t j = 0; j < 2; j++)
            EXPECT_EQ(a(i, j), keep(i, j));
    EXPECT_EQ(a(4, 1), 0.0f);

    a.resize(2, 2, false);
    EXPECT_EQ(a(0, 0), 0.0f);
}

TEST(MatrixTests, prodAndTransMatchNaive){
    la::Matrix<float> a = randomMatrix(7, 5, 2);
    la::Matrix<float> b = randomMatrix(5, 4, 3);
    la::Matrix<float> c = la::prod(a, b);
    ASSERT_EQ(c.size1(), 7u);
    ASSERT_EQ(c.size2(), 4u);
    for (size_t i = 0; i < 7; i++)
        for (size_t j = 0; j < 4; j++){
            double s = 0.0;
            for (size_t k = 0; k < 5; k++) s += (double)a(i, k) * b(k, j);
            EXPECT_NEAR(c(i, j), s, 1e-5);
        }

    la::Matrix<float> at = la::trans(a);
    ASSERT_EQ(at.size1(), 5u);
    ASSERT_EQ(at.size2(), 7u);
    EXPECT_EQ(at(4, 6), a(6, 4));

    // (AB)^T = B^T A^T
    la::Matrix<float> lhs = la::trans(c);
    la::Matrix<float> rhs = la::prod(la::trans(b), la::trans(a));
    for (size_t i = 0; i < 4; i++)
        for (size_t j = 0; j < 7; j++)
            EXPECT_NEAR(lhs(i, j), rhs(i, j), 1e-5);

    EXPECT_THROW(la::prod(a, a), std::invalid_argument);
}

TEST(MatrixTests, swapAndAliasedAssign){
    la::Matrix<float> a = randomMatrix(3, 3, 4);
    la::Matrix<float> b(2, 5, 2.0f);
    la::Matrix<float> a0(a);
    a.swap(b);
    EXPECT_EQ(a.size2(), 5u);
    EXPECT_EQ(b(2, 2), a0(2, 2));

    // A = prod(B, A) must use the old A (uBLAS semantics without noalias)
    la::Matrix<float> m = randomMatrix(3, 3, 5);
    la::Matrix<float> expected = la::prod(b, m);
    m = la::prod(b, m);
    for (size_t i = 0; i < 3; i++)
        for (size_t j = 0; j < 3; j++)
            EXPECT_EQ(m(i, j), expected(i, j));
}

TEST(MatrixTests, streamFormatMatchesUblas){
    la::Matrix<float> a(2, 2);
    a(0, 0) = 1; a(0, 1) = 2; a(1, 0) = 3; a(1, 1) = 4;
    std::ostringstream os;
    os << a;
    EXPECT_EQ(os.str(), "[2,2]((1,2),(3,4))");
}

// SVD from SVD.h on the new matrix type: A = U S V^T and U, V orthogonal
TEST(MatrixTests, svdReconstructs){
    const size_t r = 8, c = 5;
    la::Matrix<float> a = randomMatrix(r, c, 6);
    la::Matrix<float> work(a), um, sm, vm;
    SVD svd;
    svd.svd(work, um, sm, vm);

    la::Matrix<float> back = la::prod(um, la::prod(sm, la::trans(vm)));
    for (size_t i = 0; i < r; i++)
        for (size_t j = 0; j < c; j++)
            EXPECT_NEAR(back(i, j), a(i, j), 1e-4);

    la::Matrix<float> utu = la::prod(la::trans(um), um);
    for (size_t i = 0; i < r; i++)
        for (size_t j = 0; j < r; j++)
            EXPECT_NEAR(utu(i, j), i == j ? 1.0f : 0.0f, 1e-4);

}
