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

#ifndef WETSAXS_WATERS_H
#define WETSAXS_WATERS_H

#include <map>
#include "SphericalHarmonics.h"
#include "SphericalBessels.h"
#include "AtomisticModel.h"
#include <sastools/svd3.h>
#include <sastools/vector3.h>
#include <stdexcept>
#include <random>

class AtomisticModel;

class Waters {
private:

    const int totalTempWaters = 30;
    unsigned int totalwaters, totalwatersInExcludedVolume;
    float rotation[9];
    bool writeWaters = false;

    std::map < std::string, float> minima;

    std::map < std::string, std::unordered_map<std::string, Coords> > sideChains;
    std::map <std::string, std::string> centeringAtoms;

private:
    std::map < std::string, std::vector<Coords> >waters;

    std::vector<Coords> hydration;

    std::vector<float> rvalues;

private:
    std::vector<float> thetas; // could be aligned memory
    std::vector<float> phis;

    std::vector<float> almRealSum; // partial sums, bound water A_w
    std::vector<float> almImagSum;
    std::vector<float> almRealDisplacedSum; // partial sums, displaced bulk water at the sites A_w0
    std::vector<float> almImagDisplacedSum;

    SphericalHarmonics she;
    SphericalBessels sbj;
    unsigned int lmax;

    float convert = (M_PI/180.0f);
    std::mt19937 gen;
    std::uniform_int_distribution<> randomIndex;//(0,360); // guaranteed unbiased
    std::uniform_int_distribution<> randomBeta;//(0,180); // guaranteed unbiased

    Coords * tempWaters;

public:

    Waters();

    void hydrateAtomisticModel(AtomisticModel & model);

    void extractWatersFromAtomisticModel(AtomisticModel & model);

    void createBackbonesAndSideChains();

    void hydrateResidueDirect(std::string residue, int resid, int atomsInResidue, unsigned int startAt,
                              const std::string *atomType,
                              const float *xpos,
                              const float *ypos,
                              const float *zpos);

    void resetTempWaters();

    void setRotationMatrix(int totalAtomsInResidue,
                           float *ref_x,
                           float *ref_y,
                           float *ref_z,
                           float *tar_x,
                           float *tar_y,
                           float *tar_z);

    int placeWaters(std::string residue, int resid, float delX, float delY, float delZ);

    [[nodiscard]] const std::map<std::string, std::vector<Coords>> & getWaters() const {
        return waters;
    }

    [[nodiscard]] const std::unordered_map<std::string, Coords> & getSideChainByName(const std::string& name) const {
        auto findIt = sideChains.find(name);
        return findIt->second;
//        return (findIt != sideChains.end()) ? (findIt->second) : std::make_pair("", Coords(0,0,0,"",0));
    }


    void writeWatersToFile(std::string name);

    void createSphericalCoordinateOfHydration();

    void populateAmplitudesExcludedVolume(unsigned int lmax, unsigned int totalqvalues, std::vector<float> &qvalues);

    void createSphericalHarmonics(unsigned int lmax);

    void createSphericalBessels(unsigned int totalqvalues, std::vector<float> &qvalues);

    std::string assignBackbone(std::string residue);

    void calculatePartialAmplitudes(unsigned int lmax,
                                    unsigned int totalqvalues,
                                    std::vector<float> &qvalues,
                                    bool recalculate = false);

    const float * getModelPartialsReal() { return almRealSum.data(); } // partial sums
    const float * getModelPartialsImag() { return almImagSum.data();}
    const float * getDisplacedPartialsReal() { return almRealDisplacedSum.data(); } // A_w0, buffer side
    const float * getDisplacedPartialsImag() { return almImagDisplacedSum.data(); }

    [[nodiscard]] const std::map<std::string, std::unordered_map<std::string, Coords>> & getSideChains() const {
        return sideChains;
    }

    float getForwardScatteringOfHydrationModel();

    [[nodiscard]] unsigned int getTotalWaters() const { return totalwaters; }

    [[nodiscard]] const std::vector<float> & getRvalues() const {
        return rvalues;
    }
    [[nodiscard]] const std::vector<float> & getPhiValues() const {
        return phis;
    }
    [[nodiscard]] const std::vector<float> & getThetaValues() const {
        return thetas;
    }


    void translateAndWriteWatersToFile(std::string name, const vector3 *);

    SphericalHarmonics & getSHE(){ return she;}
    SphericalBessels & getSBJ() { return sbj;}

    const std::vector<Coords> & getHydration(){ return hydration; }

    void
    rotate(float alpha, float beta, float gamma, std::vector<vector3> & tetrahedron,
           std::vector<vector3> & transformedCoordinates,
           vector3 &translateToVector);
};


#endif //WETSAXS_WATERS_H
