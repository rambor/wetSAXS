//
// Created by xos81802 on 30/07/2026.
//
#include <sastools/PDBModel.h>
#include "gtest/gtest.h"
#include "support.hpp"
#include "Mask.h"
#include "Waters.h"
#include "omp.h"

class MaskTests : public ::testing::Test {

public:
    PDBModel bsaModel;
    // setup
    MaskTests() : ::testing::Test(),
                      bsaModel( PDBModel(fixture(bsa.pdb), false, false) ){
    }
};


TEST_F(MaskTests, maskBSATest) {

    Waters waters = Waters();
    AtomisticModel md = AtomisticModel(fixture(bsa.pdb), false, false);
    waters.hydrateAtomisticModel(md);
    waters.createSphericalCoordinateOfHydration();

    Mask mask = Mask();
    mask.readInWaterModel(waters);
    mask.createParticleMask(md);

    mask.writeMaskToPDB();

}

TEST_F(MaskTests, maskBSAQValueTest) {

    Waters waters = Waters();
    AtomisticModel md = AtomisticModel(fixture(bsa.pdb), false, false);
    waters.hydrateAtomisticModel(md);
    waters.createSphericalCoordinateOfHydration();

    Mask mask = Mask();
    mask.readInWaterModel(waters);
    mask.createParticleMask(md);

    float cross_Ap_Am, cross_Ax_Am;

    mask.calculateMaskCrossTerm(0.1, md, &cross_Ap_Am, &cross_Ax_Am);

    mask.writeMaskToPDB();

}

TEST_F(MaskTests, maskBSAQValuesTest) {

    Waters waters = Waters();
    AtomisticModel md = AtomisticModel(fixture(bsa.pdb), false, false);
    waters.hydrateAtomisticModel(md);
    waters.createSphericalCoordinateOfHydration();

    Mask mask = Mask();
    mask.readInWaterModel(waters);
    mask.createParticleMask(md);

    float diff = (0.5-0.001)/700;
    float q_at = 0.001;
    std::vector<float> qvalues;
    while (q_at < 0.5) {
        qvalues.push_back(q_at);
        q_at += diff;
    }

    unsigned int totalq = qvalues.size();
    // Initialize the std::vector v by arr

    std::clock_t startTime = std::clock();
    mask.calculateMaskCrossTerm(totalq, qvalues, md);

    double runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;
    std::cout<< " RUN TIME: " << runtime << std::endl;
    mask.writeMaskToPDB();

//     const int N = 100'000'000;
//     std::vector<double> vec(N);
//
//     startTime = std::clock();
//     int i;
// #pragma omp parallel for default(none) num_threads(10) private(i) shared(vec, N)
//     for( i=0; i< N; ++i) {
//         vec[i] = sin(M_PI*float(i)/N);
//     }
//
//     runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;
//     std::cout<< " RUN TIME looped: " << runtime << " " << vec.size() << std::endl;

}