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
// Created by Robert Rambo on 04/05/2022.
//

#include <sastools/PDBModel.h>
#include "gtest/gtest.h"
#include "support.hpp"
#include "AtomisticModel.h"
#include "Waters.h"

class AtomisticModelTests : public ::testing::Test {

public:
    PDBModel bsaModel;
    // setup
    AtomisticModelTests() : ::testing::Test(),
                      bsaModel( PDBModel(fixture(bsa.pdb), false, false) ){
    }
};

TEST_F(AtomisticModelTests, getFilename){
    EXPECT_EQ(bsaModel.getFilename(), "bsa.pdb") << "File name should be bsa.pdb " << bsaModel.getFilename();
}

TEST_F(AtomisticModelTests, getTotalAtoms){
    AtomisticModel md = AtomisticModel(fixture(bsa.pdb), false, false);
    EXPECT_EQ(4682, md.getTotalAtoms()) << "Total coordinates should be " << bsaModel.getTotalCoordinates();
}

TEST_F(AtomisticModelTests, getTotalUniqueAtoms){
    AtomisticModel md = AtomisticModel(fixture(bsa.pdb), false, false);

    EXPECT_EQ(5, md.getTotalUniqueAtoms()) << "Total Unique atoms should be " << bsaModel.getTotalUniqueAtoms();
}

TEST_F(AtomisticModelTests, calculatePartials){
    AtomisticModel md = AtomisticModel(fixture(bsa.pdb), false, false);

    unsigned int lmax = 2;
    unsigned int totalqvalues = 5;
    std::vector < float > qvalues(totalqvalues);
    qvalues[0] = 0.0001;
    qvalues[1] = 0.001;
    qvalues[2] = 0.01;
    qvalues[3] = 0.1;
    qvalues[4] = 0.41;

    md.calculatePartialAmplitudes(lmax, totalqvalues, qvalues);

    //EXPECT_EQ(5, md.getTotalUniqueAtoms()) << "Total Unique atoms should be " << bsaModel.getTotalUniqueAtoms();
}


TEST_F(AtomisticModelTests, calculateMicelleInVacuo){
    float qmax = 0.67;
    float qmin = 0.002;

    std::vector<float> qvalues;
    float delta = (qmax-qmin)/(float)601;

    float val = qmin;
    for(int k=0; k<601; k++){
        qvalues.push_back(val);
        val += delta;
    }

    AtomisticModel md = AtomisticModel(fixture(chen_DDM.pdb), false, false);
    // load first model, need it for dmax
    // if dmax is not specified, must retrieve from PDB
    // dmax should really be based on the empirical value not model
    float dmax = md.getDmax();
    SASTOOLS_UTILS_H::logger("DMAX SET", std::to_string(dmax));

    unsigned lmax = (unsigned int)((qmax*dmax)/SIMDE_MATH_PI) + 1;
    SASTOOLS_UTILS_H::logger("LMAX", std::to_string(lmax));
    SASTOOLS_UTILS_H::logger("", "Calculating Amplitudes");

    md.writeInVacuoModel( lmax, qvalues.size(), qvalues);

}

TEST_F(AtomisticModelTests, calculateDNAInVacuo){
    float qmax = 0.67;
    float qmin = 0.002;

    std::vector<float> qvalues;
    float delta = (qmax-qmin)/(float)601;
    for (float val = qmin; val<qmax; val+=delta){
        qvalues.push_back(val);
    }


    AtomisticModel md = AtomisticModel(fixture(25dsDNA.pdb), false, false);
    // load first model, need it for dmax

    // if dmax is not specified, must retrieve from PDB
    // dmax should really be based on the empirical value not model
    float dmax = md.getDmax();


    SASTOOLS_UTILS_H::logger("DMAX SET", std::to_string(dmax));

    unsigned lmax = (unsigned int)((qmax*dmax)/SIMDE_MATH_PI) + 1;
    SASTOOLS_UTILS_H::logger("LMAX", std::to_string(lmax));
    SASTOOLS_UTILS_H::logger("", "Calculating Amplitudes");

    md.writeInVacuoModel( lmax, qvalues.size(), qvalues);
}

// Calculate partials only calculates for positive m values, we use symmetry to flip signs for negative m values
TEST_F(AtomisticModelTests, validatePartials){

    AtomisticModel model = AtomisticModel(fixture(bsa.pdb), false, false);
    int totalAtoms = model.getTotalAtoms();

    auto atomList = model.getAtomList();
    std::map<int, float> asf;
    for(auto & atom : atomList){
        asf.insert(std::pair<int, float>(atom, 0.0f));
    }

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

    unsigned int lmax = (unsigned int)((qmax*model.getDmax())/SIMDE_MATH_PI) + 1;
    model.calculatePartialAmplitudes(lmax, qvalues.size(), qvalues);

    // get the partials to check they are correct
    const float * ptrRAlm = model.getModelPartialsReal();
    const float * ptrIAlm = model.getModelPartialsImag();
    const float * ptrREx = model.getExcludedVolumePartialsReal();
    const float * ptrIEx = model.getExcludedVolumePartialsImag();

    float *pYLM;

    int lmax_plus_1 = lmax + 1;
    unsigned int lmax_plus_1_totalAtoms = lmax_plus_1*totalAtoms;

    SphericalHarmonics & she = model.getSHE();
    auto pReal = she.getpDataYLMReal();
    auto pImag = she.getpDataYLMImag();

    SphericalBessels & sbj = model.getSBJ();
    float * pSBJ = sbj.getpSphericalBessels(), * pBessel;

    float bessel_product, bessel_product_Ex;

    auto pAtomicNumbers = model.getPDBModel().getAtomicNumberVec();
    const float * pfractionOccAtoms = model.getFractionalVolume();

    // validate partials
    int q_index=0;
    int count=0;
    for(auto & qvalue : qvalues){

        for(auto & ff : asf){ // calculate atomic scattering form factors for given qvalue
            ff.second = functions_WETSAXS::asf(ff.first, qvalue);
        }

        float waterASF = functions_WETSAXS::asf(99, qvalue);

        // populate qr values

        for(int l=0; l<=lmax; l++){

            int bessel_index = lmax_plus_1_totalAtoms*q_index + totalAtoms*l;
            int base_index = (l*l)*totalAtoms; // always positive

            // biggest time saving will be here, need to exploit the symmetry of spherical hamormonics (for m = -m)
            for(int m=-l; m<=l; m++) {
                // for given (l,m) calculate Spherical Harmonic
                float partialSumRAp = 0.0f;
                float partialSumIAp = 0.0f;

                // excluded volume term based on atomic coordinates
                float partialSumREx = 0.0f;
                float partialSumIEx = 0.0f;

                int ylm_index = base_index + (l + m) * totalAtoms; //

                // for each (l,m) calculate partial sum at given q value
                for (unsigned int n = 0; n < totalAtoms; n++) { // over all atoms

                    pBessel = &pSBJ[bessel_index + n];

                    bessel_product = asf[pAtomicNumbers[n]] * *pBessel; // for each m value, using same bessel values at l
                    bessel_product_Ex = pfractionOccAtoms[n] * *pBessel; // for estimating excluded volume using water ASF

                    // real term
                    pYLM = &pReal[ylm_index];
                    partialSumRAp += bessel_product * *pYLM;     // particle, real
                    partialSumREx += bessel_product_Ex * *pYLM;  // excluded volume, real

                    //imaginary term
                    pYLM = &pImag[ylm_index];
                    partialSumIAp += bessel_product * *pYLM;    // particle, imaginary
                    partialSumIEx += bessel_product_Ex * *pYLM; // excluded volume, imaginary

                    ++ylm_index;
                }

                partialSumREx *= waterASF;
                partialSumIEx *= waterASF;

                float ulp = abs((ptrRAlm[count] - partialSumRAp)/partialSumRAp);
                EXPECT_LE(ulp, 0.0001) << "error RAp " << count << " l " << l << " m " << m << " " << ptrRAlm[count]  << " " << partialSumRAp << " " << ulp;
                ulp = abs((ptrIAlm[count] - partialSumIAp)/partialSumIAp);
                EXPECT_LE(ulp, 0.0001) << "error IAp " << count << " l " << l << " m " << m << " " << ptrIAlm[count]  << " " << partialSumIAp << " " << ulp;

//                EXPECT_NEAR(ptrREx[count], partialSumREx, 0.0001) << "error REx " << count << " l " << l << " m " << m << " " << ptrREx[count]  << " " << partialSumREx;
                ulp = abs((ptrREx[count] - partialSumREx)/partialSumREx);
                EXPECT_LE(ulp, 0.0001) << "error REx " << count << " l " << l << " m " << m << " " << ptrREx[count]  << " " << partialSumREx << " " << ulp;

                ulp = abs((ptrIEx[count] - partialSumIEx)/partialSumIEx);
                EXPECT_LE(ulp, 0.0001) << "error IEx " << count << " l " << l << " m " << m << " " << ptrIEx[count]  << " " << partialSumIEx << " " << ulp;

                ++count;
            }
        }
        ++q_index;
    }

}



