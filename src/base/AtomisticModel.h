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
#include <sastools/Bead.h>
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
    std::vector<float> bead_rvalues;
    std::vector<float> bead_thetas; // could be aligned memory
    std::vector<float> bead_phis;

    std::vector<float> binned_distances;
    std::vector<float> almRealSum; // partial sums
    std::vector<float> almImagSum;
    std::vector<float> almRealExcludedVolSum; // partial sums
    std::vector<float> almImagExcludedVolSum;

    std::vector<float> clmNonNeighbors, clmNeighbors;

    std::vector<float> binnedPartials;

    std::set <int> atomList;

    //std::map <std::string, int> hydrogens;
    std::map<unsigned int, std::vector<unsigned int>> neighbors;
    std::map<unsigned int, std::vector<unsigned int>> non_neighbors;

    std::vector<float> fractionOccAtoms;
    std::vector<Bead> beads;

    SphericalHarmonics she, she_x;
    SphericalBessels sbj, sbj_x;

    float rg = 0;

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
        bead_rvalues = std::move(model.bead_rvalues);
        bead_thetas = std::move(model.bead_thetas);
        bead_phis = std::move(model.bead_phis);

        binned_distances = std::move(model.binned_distances);
        almRealSum = std::move(model.almRealSum);
        almImagSum = std::move(model.almImagSum);

        almRealExcludedVolSum = std::move(model.almRealExcludedVolSum); // partial sums
        almImagExcludedVolSum = std::move(model.almImagExcludedVolSum);
        binnedPartials = std::move(model.binnedPartials);

        clmNonNeighbors = std::move(model.clmNonNeighbors);
        clmNeighbors = std::move(model.clmNeighbors);

        atomList = std::move(model.atomList);
//        hydrogens = std::move(model.hydrogens);
        neighbors = std::move(model.neighbors);
        non_neighbors = std::move(model.non_neighbors);

        fractionOccAtoms = std::move(model.fractionOccAtoms);

        beads = std::move(model.beads);


        she = std::move(model.she);
        sbj = std::move(model.sbj);

        rg = std::move(model.rg);

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

    const float * getExcludedVolumeIntensitiesNeighbors(){ return clmNeighbors.data(); }
    const float * getExcludedVolumeIntensitiesNonNeighbors(){ return clmNonNeighbors.data(); }

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

    void setMatchPointScale(int * pAtomicNumbers);
    void setMatchPointScale(unsigned int totalWatersInHydration);

    float getMatchPoint(){ return matchPoint;}

    void createHCPGridModel(float bin_width);

    void populateBinnedDistances(unsigned int lmax);

    void populateSphBesselsArray(int l, std::vector<float> &qr, float *bessel_alq);

    void setFractionalVolume(float volume){

        fractionOccAtoms.resize(totalAtoms);
        auto pAtomicVolumes = pdbModel.getAtomicVolumeVec();

        float inv = 1.0f/volume;

        for(unsigned int i=0; i<totalAtoms; i++){
            fractionOccAtoms[i] =  pAtomicVolumes[i]*inv; // normalized to volume of water, looking at fractions of water volume
        }
    }

    /*
     * Convert the Gerstein-Voronoi atomic radii into Gaussian radii
     *
     */
//    void setGaussianRadii(){
//        auto pAtomicVolumes = pdbModel.getAtomicVolumeVec();
//
//        gaussian_radii.resize(totalAtoms);
//        float * pvol;
//        float invpi3halves = 1.0f/std::sqrtf(M_PI*M_PI*M_PI);
//
//        for(unsigned int i=0; i<totalAtoms; i++){
//            pvol = &pAtomicVolumes[i];
//            gaussian_radii[i] =  std::cbrt(*pvol * invpi3halves); // normalized to volume of water, looking at fractions of water volume
//        }
//    }

    /*
     * Gaussian volumes where v = r^3 * pi^3/2
     */
    void setFractionalVolumeAtQ(float q_val){

        fractionOccAtoms.resize(totalAtoms);

        float * gaussian_radii = pdbModel.getAtomicGaussianRadii();

        float factor, diff, inv = 1.0f/(4*M_PI);
        float pw_vol = 29.9;
        float qval2 = q_val*q_val;

        for(unsigned int i=0; i<totalAtoms; i++){

            diff = gaussian_radii[i]/1.75f;
            factor = pw_vol*(1 - diff*diff*diff);

            fractionOccAtoms[i] = 0.334*factor * std::expf(-std::cbrt(factor*factor)*inv*qval2); // normalized to volume of water, looking at fractions of water volume
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

    void createNeighborhoods(float cutoff);

    void createSphericalHarmonicsEx();

    void createSphericalBesselsEx(unsigned int totalqvalues, std::vector<float> &qvalues);

    void calculateExVolPartialAmplitudes(unsigned int lmax, unsigned int totalqvalues, std::vector<float> &qvalues,
                                         bool recalculate);


    float calculateRg(){
        const int * pAtomicNumbers = pdbModel.getAtomicNumberVec();
        rg = 0;
        for(unsigned int n=0; n < totalAtoms; n++){

            float rval = rvalues[n];
            rg += rval*rval; // unique atomic numbers
        }
        rg *= 1/(float)totalAtoms;
        rg = std::sqrtf(rg);
        return rg;
    }

    float getRg(){ return (rg ==0) ? this->calculateRg(): rg ;}

};


#endif //WETSAXS_ATOMISTICMODEL_H
