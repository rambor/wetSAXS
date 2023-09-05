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

#ifndef WETSAXS_ATOMISTICMODEL_H
#define WETSAXS_ATOMISTICMODEL_H

#include <vector>
#include <string>
#include <sastools/vector3.h>
#include <sastools/PDBModel.h>
#include <set>
#include "SphericalHarmonics.h"
#include "SphericalBessels.h"
#include "functions_WETSAXS.h"


class AtomisticModel {
private:

    std::string filename;
    bool forceRNA = false;
    bool keepWaters = false;

    PDBModel pdbModel;

    unsigned int totalAtoms;  // maximum number of waters per residue, should match declaration in waters
    unsigned int totalDistances;
    unsigned int lmax;

    float matchPoint = 1.0;
    float forwardScatteringParticle, forwardScatteringExVolume;

    // setup for dynamically allocated arrays
    std::vector<float> rvalues;
    std::vector<float> thetas; // could be aligned memory
    std::vector<float> phis;
    std::vector<float> binned_distances;
    std::vector<float> almRealSum; // partial sums
    std::vector<float> almImagSum;
    std::vector<float> almRealExcludedVolSum; // partial sums
    std::vector<float> almImagExcludedVolSum;
    std::vector<float> binnedPartials;

    std::set <int> atomList;

    std::map <std::string, int> hydrogens;

    std::vector<float> fractionOccAtoms;

    SphericalHarmonics she;
    SphericalBessels sbj;

public:
    AtomisticModel()=default;

    AtomisticModel(std::string filename, bool forceRNA, bool keepWaters);

    // copy constructor - prevent copying
    AtomisticModel(const AtomisticModel & model)= delete;
    // copy assignment - prevent copying
    AtomisticModel & operator=(const AtomisticModel & model)= delete;

    // move assignment operator
    AtomisticModel & operator=(AtomisticModel && model) noexcept {

        if (&model == this){
            return *this;
        }

        filename = std::move(model.filename);
        forceRNA = model.forceRNA;
        keepWaters = model.keepWaters;

        pdbModel = std::move(model.pdbModel);
        totalAtoms = model.totalAtoms;
        totalDistances = model.totalDistances;

        lmax = model.lmax;

        matchPoint = model.matchPoint;

        rvalues = std::move(model.rvalues);
        thetas = std::move(model.thetas);
        phis = std::move(model.phis);
        binned_distances = std::move(model.binned_distances);
        almRealSum = std::move(model.almRealSum);
        almImagSum = std::move(model.almImagSum);

        almRealExcludedVolSum = std::move(model.almRealExcludedVolSum); // partial sums
        almImagExcludedVolSum = std::move(model.almImagExcludedVolSum);
        binnedPartials = std::move(model.binnedPartials);

        atomList = std::move(model.atomList);
        hydrogens = std::move(model.hydrogens);

        fractionOccAtoms = std::move(model.fractionOccAtoms);

        she = std::move(model.she);
        sbj = std::move(model.sbj);

        return *this;
    }

    AtomisticModel (AtomisticModel && model) noexcept {
        *this = std::move(model);
    }

    ~AtomisticModel()=default;

    /*
     * total atoms in calculation should equal size of rvalues, thetas and phis
     */
    unsigned int getTotalAtoms(){ return totalAtoms;}

    /**
     * total number of unique atoms by atom type or atomic number
     * @return
     */
    unsigned int getTotalUniqueAtoms(){ return atomList.size();}
    std::set <int> & getAtomList(){ return atomList; }

    const float * getThetas(){ return thetas.data();}
    const float * getPhis(){ return phis.data();}
    const float * getRValues() { return rvalues.data();}

    const float * getModelPartialsReal() { return almRealSum.data(); } // partial sums
    const float * getModelPartialsImag() { return almImagSum.data();}

    const float * getExcludedVolumePartialsReal() { return almRealExcludedVolSum.data();} // partial sums
    const float * getExcludedVolumePartialsImag(){ return almImagExcludedVolSum.data();}

    const float * getBinnedPartials() { return binnedPartials.data();}

    std::vector<float> * getRValuesVector() { return &rvalues;}

    void createSphericalHarmonics(unsigned int lmax);

    void createSphericalBessels(unsigned int totalqvalues, std::vector<float> &qvalues);

    const SphericalBessels &getSbj() const;

    PDBModel & getPDBModel() { return pdbModel; }
    //const PDBModel * getPdbModel() const { return &pdbModel; }

    float getDmax(){ return pdbModel.getDmax(); }

    void calculatePartialAmplitudes(unsigned int lmax,
                                    unsigned int totalqvalues,
                                    std::vector<float> &qvalues,
                                    bool recalculate = false);

    void writeInVacuoModel(unsigned int lmax, unsigned int totalqvalues, std::vector<float> &qvalues);

    void calculatePartialAmplitudesNoExcludedVolume(unsigned int lmax,
                                               unsigned int totalqvalues,
                                               std::vector<float> &qvalues);

    void setMatchPointScale(int * pAtomicNumbers, float * fractionalOccupancies);
    void setMatchPointScale(unsigned int totalWatersInHydration);

    float getMatchPoint(){ return matchPoint;}

    void populateBinnedDistances(unsigned int lmax);

    inline void babinet_at_q(float qvalue, std::vector<float> & babinets, float * pAtomicVol, std::vector<float> & cbrft_vol2 ){

        float q2 = qvalue*qvalue*M_PI;
        for(unsigned int i=0; i<totalAtoms; i++){ // this is a function of q
            babinets[i] = pAtomicVol[i]*expf(-q2*cbrft_vol2[i]);
        }
    }

    void populateSphBesselsArray(int l, std::vector<float> &qr, float *bessel_alq);

    void setFractionalVolume(float value){
        fractionOccAtoms.resize(totalAtoms);
        auto pAtomicVolumes = pdbModel.getAtomicVolumeVec();
        float inv = 1.0f/value;

        for(unsigned int i=0; i<totalAtoms; i++){
            fractionOccAtoms[i] =  pAtomicVolumes[i]*inv; // normalized to volume of water, looking at fractions of water volume
        }
    }


    const float * getFractionalVolume() {

        if (fractionOccAtoms.empty()){
            setFractionalVolume(29.9); // 29.9 volume for H2O
        }

        return fractionOccAtoms.data();
    }


    SphericalHarmonics & getSHE(){ return she;}

    SphericalBessels & getSBJ() { return sbj;}
};


#endif //WETSAXS_ATOMISTICMODEL_H
