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
// Created by Robert Rambo on 14/06/2022.
//

#include "gtest/gtest.h"
#include "support.hpp"
#include <Waters.h>

class WatersTests : public ::testing::Test {

public:
    PDBModel bsaModel;
    // setup
    WatersTests() : ::testing::Test(),
                              bsaModel( PDBModel(fixture(bsa.pdb), false, false) ){
    }
};



TEST_F(WatersTests, ConstructorTest){

    Waters waters = Waters();
    auto pWaters = waters.getWaters();

    for (std::pair<std::string, std::vector<Coords>> element : pWaters) {
        std::string residue = element.first;
        auto coords = element.second;
        if (residue == "PHE"){
            ASSERT_EQ(coords.size(), 9) << residue << " PHE is 9 " << coords.size();
        } else if (residue.compare("rC") == 0){
            ASSERT_EQ(coords.size(), 17);
        }
    }
}

TEST_F(WatersTests, GetSideChainByName){

    Waters waters = Waters();
    auto pSideChain = waters.getSideChainByName("PHE");

    /*
    sideChains["PHE"].insert({"N",Coords(-0.0683,1.9033,-1.4177,"N",1.0)});
    sideChains["PHE"].insert({"CA",Coords(-0.7693,1.6893,-0.1377,"CA",1.0)});
    sideChains["PHE"].insert({"C",Coords(-1.1843,3.0113,0.4973,"C",1.0)});
    sideChains["PHE"].insert({"O",Coords(-0.4523,4.0003,0.3973,"O",1.0)});
    sideChains["PHE"].insert({"CB",Coords(0.0857,0.8453,0.8333,"CB",1.0)});
    sideChains["PHE"].insert({"CG",Coords(0.2457,-0.5797,0.3863,"CG",1.0)});
    sideChains["PHE"].insert({"CD1",Coords(-0.8303,-1.4677,0.4413,"CD1",1.0)});
    sideChains["PHE"].insert({"CD2",Coords(1.4687,-1.0317,-0.0977,"CD2",1.0)});
    sideChains["PHE"].insert({"CE1",Coords(-0.6763,-2.7917,0.0323,"CE1",1.0)});
    sideChains["PHE"].insert({"CE2",Coords(1.6217,-2.3467,-0.5127,"CE2",1.0)});
    sideChains["PHE"].insert({"CZ",Coords(0.5587,-3.2317,-0.4217,"CZ",1.0)});
    */

    auto pIt = pSideChain.find("N");
    auto pCZIt = pSideChain.find("CZ");
    ASSERT_NEAR(pIt->second.x, -0.0683, 0.0001) << " N x pos failed";
    ASSERT_NEAR(pIt->second.y, 1.9033, 0.0001)  << " N y pos failed";
    ASSERT_NEAR(pIt->second.z, -1.4177, 0.0001) << " N z pos failed";
    ASSERT_NEAR(pCZIt->second.x, 0.5587, 0.0002) << "CZ x pos Failed";
    ASSERT_NEAR(pCZIt->second.y, -3.2317, 0.0002) << "CZ y pos Failed";
    ASSERT_NEAR(pCZIt->second.z, -0.4217, 0.0002) << "CZ z pos Failed";

}


TEST_F(WatersTests, HydrateTest){

    Waters waters = Waters();
    AtomisticModel md = AtomisticModel(fixture(bsa.pdb), false, false);
    waters.hydrateAtomisticModel(md);
    waters.createSphericalCoordinateOfHydration();
}

TEST_F(WatersTests, HydrateTestProtein2){

    Waters waters = Waters();
    AtomisticModel md = AtomisticModel(fixture(confA.pdb), false, false);
    waters.hydrateAtomisticModel(md);
    waters.writeWatersToFile("confA_old");
    waters.createSphericalCoordinateOfHydration();

    md.getPDBModel().writeCenteredCoordinatesToFile("confA_centered");
}


TEST_F(WatersTests, HydrateTestRNA){

    Waters waters = Waters();
    AtomisticModel md = AtomisticModel(fixture(p4p6.pdb), true, false);
    waters.hydrateAtomisticModel(md);

    waters.writeWatersToFile("RNA");
    waters.createSphericalCoordinateOfHydration();

    md.getPDBModel().writeCenteredCoordinatesToFile("p4p6_centered");
}

TEST_F(WatersTests, CheckDistances){

    Waters waters = Waters();
    AtomisticModel md = AtomisticModel(fixture(p4p6.pdb), true, false);

    auto sidechains = waters.getSideChains();
    auto waterlist = waters.getWaters();

    for(auto & chain : sidechains){
        auto waterCoords = waterlist[chain.first];
        auto atomcoords = chain.second;
        unsigned int total = waterCoords.size();
        float max = FLT_MAX;
        for(auto & atom : atomcoords){
            auto & c1 = atom.second;
            for(unsigned int j=0; j<total; j++){
                auto & c2 = waterCoords[j];
                float value = (vector3(c1.x, c1.y, c1.z) - vector3(c2.x, c2.y, c2.z)).length();
                if (value < max){
                    max = value;
                }
            }
        }
    }
}

TEST_F(WatersTests, calculateHydratedMicelle){

    float qmax = 0.67;
    float qmin = 0.002;

    std::vector<float> qvalues;
    float delta = (qmax-qmin)/(float)601;
    float val = qmin;
    while (val < qmax){
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

    Waters waters = Waters();
    waters.hydrateAtomisticModel(md); // add waters to PDB

    waters.translateAndWriteWatersToFile("micelle_hydration_model", md.getPDBModel().getCenteringVector());

}


TEST_F(WatersTests, validatePartials){

    AtomisticModel model = AtomisticModel(fixture(bsa.pdb), false, false);
    Waters waters = Waters();
    waters.hydrateAtomisticModel(model);
    waters.createSphericalCoordinateOfHydration();

    int totalwaters = waters.getTotalWaters();

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
    waters.calculatePartialAmplitudes(lmax, qvalues.size(), qvalues);

    // get the partials to check they are correct
    const float * ptrRAlm = waters.getModelPartialsReal();
    const float * ptrIAlm = waters.getModelPartialsImag();

    int lmax_plus_1 = lmax + 1;
    unsigned int lmax_plus_1_totalAtoms = lmax_plus_1*totalwaters;

    SphericalHarmonics & she = waters.getSHE();
    auto pRealYlm = she.getpDataYLMReal();
    auto pImagYlm = she.getpDataYLMImag();

    SphericalBessels & sbj = waters.getSBJ();
    float * pSBJ = sbj.getpSphericalBessels();

    float bessel_product;

    // validate partials
    int q_index=0;
    int count=0;

    std::vector<Coords> hydration = waters.getHydration();
    float occupancies [totalwaters];
    for (unsigned int i=0; i< totalwaters; i++){
        occupancies[i] = hydration[i].occ;
    }

    for(auto & qvalue : qvalues){

        float asf_at_q = functions_WETSAXS::asf(99, qvalue);

        for(int l=0; l<=lmax; l++){

            int bessel_index = lmax_plus_1_totalAtoms*q_index + totalwaters*l;
            int base_index = (l*l)*totalwaters;

            for(int m=-l; m<=l; m++){
                // for given (l,m) calculate Spherical Harmonic
                // she.calculateSHE(l,m,thetas.data(), phis.data());

                // over all atoms
                float partialSumR = 0.0f;
                float partialSumI = 0.0f;

                int ylm_index = base_index + (l+m)*totalwaters;

                // for each (l,m) calculate partial sum at given q value
                for(unsigned int n=0; n<totalwaters; n++){

                    bessel_product = pSBJ[bessel_index + n] * occupancies[n];

                    // real term
                    partialSumR += bessel_product * pRealYlm[ylm_index];
                    //imaginary term
                    partialSumI += bessel_product * pImagYlm[ylm_index];

                    ++ylm_index;
                }

                partialSumR *= asf_at_q;
                partialSumI *= asf_at_q;

                // particle scattering
                EXPECT_EQ(ptrRAlm[count], partialSumR) << "error RAp " << count << " l " << l << " m " << m << " " << ptrRAlm[count]  << " " << partialSumR;
                EXPECT_EQ(ptrIAlm[count], partialSumI) << "error IAp " << count << " l " << l << " m " << m;

                ++count;
            }
        }

        q_index++;
    }

}