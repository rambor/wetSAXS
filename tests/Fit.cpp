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
// Created by Robert Rambo on 27/07/2022.
//

#include "gtest/gtest.h"
#include "support.hpp"
#include "../src/base/Fit.h"

class FitTests : public ::testing::Test {
// BSA_main_peak.dat
// BSA_high.dat
public:
    IofQData iofqdata, rnaIofQdata;
    // setup
    // BSA_0p7.dat
    // BSA_high.dat
    FitTests() : ::testing::Test(),
                 iofqdata(IofQData(fixture(BSA_0p7.dat), false)),
                 rnaIofQdata(IofQData(fixture(P4P6_1_merged_refined_sx.dat), false)){
    }
};



TEST_F(FitTests, calculatePartials){

    AtomisticModel md = AtomisticModel(fixture(bsa_md_1_7.pdb), false, false);

    iofqdata.extractData(); // assume dmax is not set as data from other facilities will not be compatible with Scatter format
    iofqdata.setDmax(md.getDmax());
    iofqdata.makeWorkingSet(2);

    // need q-values and lmax for the fit
    std::vector<float> qvalues = iofqdata.getWorkingSetQvalues();

    unsigned int lmax = (unsigned int)((iofqdata.getQmax()*md.getDmax())/SIMDE_MATH_PI) + 1;
    SASTOOLS_UTILS_H::logger("LMAX", std::to_string(lmax));
    SASTOOLS_UTILS_H::logger("", "Calculating Amplitudes");
    md.calculatePartialAmplitudes(lmax, iofqdata.getTotalInWorkingSet(), qvalues);

    // create water model
    Waters waters = Waters();
    waters.hydrateAtomisticModel(md); // add waters to PDB
    waters.createSphericalCoordinateOfHydration();

    waters.calculatePartialAmplitudes(lmax, iofqdata.getTotalInWorkingSet(), qvalues);
    waters.writeWatersToFile("water_model");

    Fit fitme = Fit(lmax, 500,0.1,1.8f, 120);
    fitme.search(iofqdata, md, waters);
    fitme.printBestScores();

}


TEST_F(FitTests, calculatePartialsRNA){

    AtomisticModel md = AtomisticModel(fixture(p4p6.pdb), true, false);
    AtomisticModel bent = AtomisticModel(fixture(bent_6.pdb), true, false);

    rnaIofQdata.extractData(); // assume dmax is not set as data from other facilities will not be compatible with Scatter format
    rnaIofQdata.setDmax(md.getDmax());
    rnaIofQdata.makeWorkingSet(4);

    // need q-values and lmax for the fit
    std::vector<float> qvalues = rnaIofQdata.getWorkingSetQvalues();

    unsigned int lmax = (unsigned int)((rnaIofQdata.getQmax()*md.getDmax())/SIMDE_MATH_PI) + 1;
    SASTOOLS_UTILS_H::logger("LMAX", std::to_string(lmax));
    SASTOOLS_UTILS_H::logger("", "Calculating Amplitudes for WT P4P6");
    md.calculatePartialAmplitudes(lmax, rnaIofQdata.getTotalInWorkingSet(), qvalues);

    // create water model
    Waters waters = Waters();
    waters.hydrateAtomisticModel(md); // add waters to PDB
    waters.createSphericalCoordinateOfHydration();
    waters.calculatePartialAmplitudes(lmax, rnaIofQdata.getTotalInWorkingSet(), qvalues);

    SASTOOLS_UTILS_H::logger("", "Calculating Amplitudes for BENT P4P6");
    bent.calculatePartialAmplitudes(lmax, rnaIofQdata.getTotalInWorkingSet(), qvalues);

    Waters waters_bent = Waters();
    waters_bent.hydrateAtomisticModel(bent); // add waters to PDB
    waters_bent.createSphericalCoordinateOfHydration();
    waters_bent.calculatePartialAmplitudes(lmax, rnaIofQdata.getTotalInWorkingSet(), qvalues);

    Fit fitme = Fit(lmax, 333, 0.1f,1.7, 85);
    fitme.search(rnaIofQdata, md, waters);

    SASTOOLS_UTILS_H::logger("", "BENT P4P6");
    fitme.search(rnaIofQdata, bent, waters);
}

TEST_F(FitTests, calculatePartialsAB){

    AtomisticModel mdA = AtomisticModel(fixture(confA.pdb), false, false);
    AtomisticModel mdB = AtomisticModel(fixture(confB.pdb), false, false);

    IofQData data = IofQData( fixture(confAorB.dat), false);

    data.extractData(); // assume dmax is not set as data from other facilities will not be compatible with Scatter format
    data.setDmax(mdA.getDmax());
    data.makeWorkingSet();

    // need q-values and lmax for the fit
    std::vector<float> qvalues = data.getWorkingSetQvalues();

    unsigned int lmax = (unsigned int)((data.getQmax()*mdA.getDmax())/SIMDE_MATH_PI) + 1;
    SASTOOLS_UTILS_H::logger("LMAX", std::to_string(lmax));
    SASTOOLS_UTILS_H::logger("", "Calculating Amplitudes");
    mdA.calculatePartialAmplitudes(lmax, data.getTotalInWorkingSet(), qvalues);

    SASTOOLS_UTILS_H::logger("", "Calculating Amplitudes");
    mdB.calculatePartialAmplitudes(lmax, data.getTotalInWorkingSet(), qvalues);

    // create water model
    Waters waters = Waters();
    waters.hydrateAtomisticModel(mdA); // add waters to PDB
    waters.createSphericalCoordinateOfHydration();
    waters.calculatePartialAmplitudes(lmax, data.getTotalInWorkingSet(), qvalues);

    Waters watersB = Waters();
    watersB.hydrateAtomisticModel(mdB); // add waters to PDB
    watersB.createSphericalCoordinateOfHydration();
    watersB.calculatePartialAmplitudes(lmax, data.getTotalInWorkingSet(), qvalues);

    Fit fitme = Fit(lmax, 333,0.1,1.7f, 85);
//    fitme.search(data, mdA, waters);
    fitme.chiFreeSearch(7, data, mdA, waters);

    SASTOOLS_UTILS_H::logger("", "Fitting confB");
    //fitme.search(data, mdB, watersB);
    fitme.chiFreeSearch(7, data, mdB, watersB);
}

TEST_F(FitTests, calculatePartialsDNA){

    AtomisticModel md = AtomisticModel(fixture(25dsDNA.pdb), true, false);
    AtomisticModel rna = AtomisticModel(fixture(6yyt.pdb), true, false);

    IofQData dnaIofQdata = IofQData(fixture(dsDNA.dat), false);
    dnaIofQdata.extractData(); // assume dmax is not set as data from other facilities will not be compatible with Scatter format
    dnaIofQdata.setDmax(md.getDmax());
    dnaIofQdata.makeWorkingSet(4);

    // need q-values and lmax for the fit
    std::vector<float> qvalues = dnaIofQdata.getWorkingSetQvalues();

    unsigned int lmax = (unsigned int)((dnaIofQdata.getQmax()*md.getDmax())/SIMDE_MATH_PI) + 1;
    SASTOOLS_UTILS_H::logger("LMAX", std::to_string(lmax));
    SASTOOLS_UTILS_H::logger("", "Calculating Amplitudes");
    md.calculatePartialAmplitudes(lmax, dnaIofQdata.getTotalInWorkingSet(), qvalues);

    // create water model
    Waters waters = Waters();
    waters.hydrateAtomisticModel(md); // add waters to PDB
    waters.createSphericalCoordinateOfHydration();
    waters.calculatePartialAmplitudes(lmax, dnaIofQdata.getTotalInWorkingSet(), qvalues);

    SASTOOLS_UTILS_H::logger("", "Calculating Amplitudes for BENT P4P6");
    rna.calculatePartialAmplitudes(lmax, dnaIofQdata.getTotalInWorkingSet(), qvalues);

    // create water model
    Waters watersB = Waters();
    watersB.hydrateAtomisticModel(rna); // add waters to PDB
    watersB.createSphericalCoordinateOfHydration();
    watersB.calculatePartialAmplitudes(lmax, dnaIofQdata.getTotalInWorkingSet(), qvalues);


    Fit fitme = Fit(lmax, 300, 0.1,1.7, 120);
    SASTOOLS_UTILS_H::logger("", "Fitting DNA helix");
    fitme.search(dnaIofQdata, md, waters);

    SASTOOLS_UTILS_H::logger("", "Fitting RNA helix");
    fitme.search(dnaIofQdata, rna, watersB);

}

// FITTING SCHEMES
TEST_F(FitTests, chiFreeMisFitTest){

    AtomisticModel md = AtomisticModel(fixture(6yyt.pdb), true, false);

    IofQData dnaIofQdata = IofQData(fixture(dsDNA.dat), false);
    dnaIofQdata.extractData(); // assume dmax is not set as data from other facilities will not be compatible with Scatter format
    dnaIofQdata.setDmax(md.getDmax());
    dnaIofQdata.makeWorkingSet(4);

    // need q-values and lmax for the fit
    std::vector<float> qvalues = dnaIofQdata.getWorkingSetQvalues();

    unsigned int lmax = (unsigned int)((dnaIofQdata.getQmax()*md.getDmax())/SIMDE_MATH_PI) + 1;
    SASTOOLS_UTILS_H::logger("LMAX", std::to_string(lmax));
    SASTOOLS_UTILS_H::logger("", "Calculating Amplitudes");
    md.calculatePartialAmplitudes(lmax, dnaIofQdata.getTotalInWorkingSet(), qvalues);

    // create water model
    Waters waters = Waters();
    waters.hydrateAtomisticModel(md); // add waters to PDB
    waters.createSphericalCoordinateOfHydration();
    waters.calculatePartialAmplitudes(lmax, dnaIofQdata.getTotalInWorkingSet(), qvalues);

    Fit fitme = Fit(lmax, 300, 0.1,1.7, 120);
    SASTOOLS_UTILS_H::logger("", "Fitting RNA helix");
    fitme.chiFreeSearch(17, dnaIofQdata, md, waters);

}

// FITTING SCHEMES
TEST_F(FitTests, chiFreeTest){

    AtomisticModel md = AtomisticModel(fixture(25dsDNA.pdb), true, false);

    IofQData dnaIofQdata = IofQData(fixture(dsDNA.dat), false);
    dnaIofQdata.extractData(); // assume dmax is not set as data from other facilities will not be compatible with Scatter format
    dnaIofQdata.setDmax(md.getDmax());
    dnaIofQdata.makeWorkingSet(1);

    // need q-values and lmax for the fit
    std::vector<float> qvalues = dnaIofQdata.getWorkingSetQvalues();

    unsigned int lmax = (unsigned int)((dnaIofQdata.getQmax()*md.getDmax())/SIMDE_MATH_PI) + 1;
    SASTOOLS_UTILS_H::logger("LMAX", std::to_string(lmax));
    SASTOOLS_UTILS_H::logger("", "Calculating Amplitudes");
    md.calculatePartialAmplitudes(lmax, dnaIofQdata.getTotalInWorkingSet(), qvalues);

    // create water model
    Waters waters = Waters();             //
    waters.hydrateAtomisticModel(md);     // add waters to PDB
    waters.createSphericalCoordinateOfHydration();
    waters.calculatePartialAmplitudes(lmax, dnaIofQdata.getTotalInWorkingSet(), qvalues);

    Fit fitme = Fit(lmax, 300, 0.1, 1.7, 120);
    SASTOOLS_UTILS_H::logger("", "Fitting DNA helix");
    fitme.chiFreeSearch(7, dnaIofQdata, md, waters);

    // produce final fit using bestCX and BestBfactor
    dnaIofQdata.setAllDataToWorkingSet();

    // recalculate Bessels
    qvalues = dnaIofQdata.getWorkingSetQvalues();

    // make instance of new models
    AtomisticModel mdT = AtomisticModel(fixture(25dsDNA.pdb), true, false);
    mdT.calculatePartialAmplitudes(lmax, dnaIofQdata.getTotalInWorkingSet(), qvalues);

    Waters watersT = Waters();
    watersT.hydrateAtomisticModel(mdT); // add waters to PDB
    watersT.createSphericalCoordinateOfHydration();
    watersT.calculatePartialAmplitudes(lmax, dnaIofQdata.getTotalInWorkingSet(), qvalues);

    // recalculate Partials for the initial models
    waters.calculatePartialAmplitudes(lmax, qvalues.size(), qvalues, true);
    md.calculatePartialAmplitudes(lmax, qvalues.size(), qvalues, true);

    auto mdPartialR = md.getModelPartialsReal();        // updated model
    auto mdTPartialR = mdT.getModelPartialsReal();      // new model
    auto mdPartialI = md.getModelPartialsImag();        // updated model
    auto mdTPartialI = mdT.getModelPartialsImag();      // new model

    auto mdWPartialR = waters.getModelPartialsReal();   // updated model
    auto mdWTPartialR = watersT.getModelPartialsReal(); // new model
    auto mdWPartialI = waters.getModelPartialsImag();   //
    auto mdWTPartialI = watersT.getModelPartialsImag();
//
    unsigned int totalPartials = (lmax+1)*(lmax+1)*qvalues.size();
//
    for(unsigned int i=0; i<totalPartials; i++){
        ASSERT_EQ(mdWPartialR[i], mdWTPartialR[i]);
        ASSERT_EQ(mdWPartialI[i], mdWTPartialI[i]);
    }

    for(unsigned int i=0; i<totalPartials; i++){
        ASSERT_EQ(mdPartialR[i], mdTPartialR[i]);
        ASSERT_EQ(mdPartialI[i], mdTPartialI[i]);
    }

    int totalAtomsWaters = waters.getTotalWaters();
    int bessel_index;
    int q_index = 0;
    const auto pSBJWaters = waters.getSBJ().getpSphericalBessels();
    const auto pSBJTWaters = watersT.getSBJ().getpSphericalBessels();

    int totalAtoms = md.getTotalAtoms();
    const auto pSBJ = md.getSBJ().getpSphericalBessels();
    const auto pSBJT = mdT.getSBJ().getpSphericalBessels();

    for(auto & qvalue : qvalues){
        for(int l=0; l<=lmax; l++){
            bessel_index = (lmax+1)*totalAtomsWaters*q_index + totalAtomsWaters*l;
                for(unsigned int n=0; n<totalAtomsWaters; n++){
                    ASSERT_EQ(pSBJWaters[bessel_index + n], pSBJTWaters[bessel_index + n]);
                }
        }

        for(int l=0; l<=lmax; l++){
            bessel_index = (lmax+1)*totalAtoms*q_index + totalAtoms*l;
            for(unsigned int n=0; n<totalAtoms; n++){
                ASSERT_EQ(pSBJ[bessel_index + n], pSBJT[bessel_index + n]);
            }
        }
        q_index++;
    }

    // assemble intensities for each model
}

/*
 * u52
 */
TEST_F(FitTests, cxTestsBSA){

    AtomisticModel md = AtomisticModel(fixture(bsa.pdb), false, false);

    iofqdata.extractData(); // assume dmax is not set as data from other facilities will not be compatible with Scatter format
    iofqdata.setDmax(md.getDmax());
    iofqdata.makeWorkingSet(3);

    // need q-values and lmax for the fit
    std::vector<float> qvalues = iofqdata.getWorkingSetQvalues();

    int totalQvalues = qvalues.size();

    unsigned int lmax = (unsigned int)((iofqdata.getQmax()*md.getDmax())/SIMDE_MATH_PI) + 1;
    SASTOOLS_UTILS_H::logger("LMAX", std::to_string(lmax));
    SASTOOLS_UTILS_H::logger("", "Calculating Amplitudes");

    md.calculatePartialAmplitudes(lmax, totalQvalues, qvalues);

    // create water model
    Waters waters = Waters();
    waters.hydrateAtomisticModel(md); // add waters to PDB
    waters.createSphericalCoordinateOfHydration();
    waters.calculatePartialAmplitudes(lmax, totalQvalues, qvalues);
    waters.writeWatersToFile("water_model");

    Fit fitme = Fit(lmax, 133,0.4,1.8f, 87);

    SASTOOLS_UTILS_H::logger("", "Fitting BSA");

    fitme.search(iofqdata, md, waters);

    fitme.printBestScores();
}

TEST_F(FitTests, simulateMicelle){

    float qmax = 0.67;
    float qmin = 0.002;

    std::vector<float> qvalues;
    float delta = (qmax-qmin)/(float)601;
    float val = qmin;
    while (val < qmax){
        qvalues.push_back(val);
        val += delta;
    }

    // load first model, need it for dmax
    // if dmax is not specified, must retrieve from PDB
    // dmax should really be based on the empirical value not model
    AtomisticModel md = AtomisticModel(fixture(chen_DDM.pdb), false, false);
    float dmax = md.getDmax();
    SASTOOLS_UTILS_H::logger("DMAX SET", std::to_string(dmax));

    unsigned lmax = (unsigned int)((qmax*dmax)/SIMDE_MATH_PI) + 1;
    SASTOOLS_UTILS_H::logger("LMAX", std::to_string(lmax));
    SASTOOLS_UTILS_H::logger("", "Calculating Amplitudes");
    md.calculatePartialAmplitudes(lmax, qvalues.size(), qvalues);

    Waters waters = Waters();
    waters.hydrateAtomisticModel(md); // add waters to PDB

    waters.calculatePartialAmplitudes(lmax, qvalues.size(), qvalues, false);
    waters.translateAndWriteWatersToFile("chen_DDM", md.getPDBModel().getCenteringVector());

    float cx = 0.5;
    float bfactor = 17;

    Fit simulated = Fit(lmax);
    simulated.simulate(cx, bfactor, qvalues, md, waters);

}


TEST_F(FitTests, simulateBSA){

    AtomisticModel md = AtomisticModel(fixture(bsa.pdb), false, false);
    float qmin = 0.001;
    float qmax = 0.9;
    float ns = 501;

    std::vector<float> qvalues;
    float delta = (qmax-qmin)/(float)ns;

    float val = qmin;
    while (val < qmax){
        qvalues.push_back(val);
        val += delta;
    }

    unsigned int lmax = (unsigned int)((qmax*md.getDmax())/SIMDE_MATH_PI) + 1;
    md.calculatePartialAmplitudes(lmax, qvalues.size(), qvalues);

    // hydrate the model
    Waters waters = Waters();
    waters.hydrateAtomisticModel(md); // add waters to PDB
    waters.createSphericalCoordinateOfHydration();
    waters.calculatePartialAmplitudes(lmax, qvalues.size(), qvalues);

    float cx = 0.5;
    float bfactor = 1;

    Fit simulated = Fit(lmax);
    simulated.simulate(cx, bfactor, qvalues, md, waters);

    bfactor = 300;
    simulated.simulate(cx, bfactor, qvalues, md, waters);

}



TEST_F(FitTests, cxTestsDNA){

    AtomisticModel md = AtomisticModel(fixture(25dsDNA.pdb), true, false);

    IofQData dnaIofQdata = IofQData(fixture(dsDNA.dat), false);
    dnaIofQdata.extractData(); // assume dmax is not set as data from other facilities will not be compatible with Scatter format
    dnaIofQdata.setDmax(md.getDmax());
    dnaIofQdata.makeWorkingSet(4);

    // need q-values and lmax for the fit
    //std::vector<float> qvalues = dnaIofQdata.getWorkingSetQvalues();
    std::vector<float> qvalues = dnaIofQdata.getQvalues();

    unsigned int lmax = (unsigned int)((dnaIofQdata.getQmax()*md.getDmax())/SIMDE_MATH_PI) + 1;
    SASTOOLS_UTILS_H::logger("LMAX", std::to_string(lmax));
    SASTOOLS_UTILS_H::logger("", "Calculating Amplitudes");
    md.calculatePartialAmplitudes(lmax, dnaIofQdata.getTotal(), qvalues);

    // create water model
    Waters waters = Waters();
    waters.hydrateAtomisticModel(md); // add waters to PDB
    waters.createSphericalCoordinateOfHydration();
    waters.calculatePartialAmplitudes(lmax, dnaIofQdata.getTotal(), qvalues);

    Fit fitme = Fit(lmax, 300, 0.1,1.7, 120);
    SASTOOLS_UTILS_H::logger("", "Fitting DNA helix");

}

// for deciding between competing models, should also use common cx and B-factor

// ensemble fitting, need to fit a small subset to estimate parameters, cx and B-factor, then apply this to the entire set
// calculate SAXS curves for each model (using subset) and then