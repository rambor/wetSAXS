//
// Created by xos81802 on 30/07/2026.
//
#include <sastools/PDBModel.h>
#include "gtest/gtest.h"
#include "support.hpp"
#include "FCCModel.h"
#include "Waters.h"

class FCCModelTests : public ::testing::Test {

public:
    PDBModel bsaModel;
    // setup
    FCCModelTests() : ::testing::Test(),
                      bsaModel( PDBModel(fixture(bsa.pdb), false, false) ){
    }
};

TEST_F(FCCModelTests, fccBSATest) {

    Waters waters = Waters();

    AtomisticModel md = AtomisticModel(fixture(bsa.pdb), false, false);
    waters.hydrateAtomisticModel(md);


    FCCModel fcc = FCCModel();
    fcc.addPDBModel(md);
    //fcc.addWaterModel(waters);
    fcc.writeFCCModelToPDB();

}


TEST_F(FCCModelTests, fccNumberDensityTest) {

    FCCModel fcc = FCCModel();
    float n = fcc.number_density_per_Angstrom3(1.0, 18);
    ASSERT_NEAR(n, 0.03345, 0.0001);

}