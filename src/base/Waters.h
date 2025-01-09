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

class Waters {
private:

    const int totalTempWaters = 30;
    unsigned int totalwaters, totalwatersInExcludedVolume;
    float rotation[9];
    bool writeWaters = false;

    std::map < std::string, float> minima;

    std::map < std::string, std::unordered_map<std::string, Coords> > sideChains;


private:
    std::map < std::string, std::vector<Coords> >waters;

    std::vector<Coords> hydration;

    std::vector<float> rvalues;
    std::vector<float> thetas; // could be aligned memory
    std::vector<float> phis;
    std::vector<float> almRealSum; // partial sums
    std::vector<float> almImagSum;

    SphericalHarmonics she;
    SphericalBessels sbj;
    unsigned int lmax;

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

    const std::map<std::string, std::vector<Coords>> & getWaters() const {
        return waters;
    }

    const std::unordered_map<std::string, Coords> & getSideChainByName(const std::string& name) const {
        auto findIt = sideChains.find(name);
        return findIt->second;
//        return (findIt != sideChains.end()) ? (findIt->second) : std::make_pair("", Coords(0,0,0,"",0));
    }


    void writeWatersToFile(std::string name);

    void createSphericalCoordinateOfHydration();

    void populateAmplitudesExcludedVolume(unsigned int lmax, unsigned int totalqvalues, std::vector<float> &qvalues);

    void createSphericalHarmonics(unsigned int lmax);

    void createSphericalBessels(unsigned int totalqvalues, std::vector<float> &qvalues);

    void calculatePartialAmplitudes(unsigned int lmax,
                                    unsigned int totalqvalues,
                                    std::vector<float> &qvalues,
                                    bool recalculate = false);

    const float * getModelPartialsReal() { return almRealSum.data(); } // partial sums
    const float * getModelPartialsImag() { return almImagSum.data();}

    const std::map<std::string, std::unordered_map<std::string, Coords>> & getSideChains() const {
        return sideChains;
    }

    float getForwardScatteringOfHydrationModel();

    unsigned int getTotalWaters(){ return totalwaters; }

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
