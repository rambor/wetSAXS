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
// Created by Robert Rambo on 13/07/2022.
//

#ifndef WETSAXS_FIT_H
#define WETSAXS_FIT_H

#include <sastools/IofQData.h>
#include "AtomisticModel.h"
#include "Waters.h"
#include "Result.h"
#include "Mask.h"

class Fit {

    unsigned int lmax, total_cx;
    float lowcx, highcx;
    float deltaCX;

    float bestScore;
    float bestCx, bestOcc, bestDW, bestChi, bestBackground;
    float bestScaleCoef = -1;

    float high_c_occ = 1.21, low_c_occ = 0.1;
    const float inv16Pi2 = 1.0/(16.0*M_PI_2);

    std::vector<unsigned int> partialRows; // working-set row -> partials row (empty: identical)
    std::vector<float> norm_aPs;
    std::vector<float> norm_aWs;
    std::vector<float> norm_aCs;
    std::vector<float> aXW_cross_term;
    std::vector<float> aPX_cross_term;
    std::vector<float> aPW_cross_term;

    std::vector<float> i_best;
    std::vector<float> i_obs;
    std::vector<float> i_obs_over_variance;
    std::vector<float> inv_variance;

    std::string header;

public:
    Fit(unsigned int lmax);
    Fit(unsigned int lmax, unsigned int division, float lowcx, float highcx, float low_c_occ, float high_c_occ);

    void search(IofQData & iofqdata, AtomisticModel & model, Waters & waterModel);
    void searchWithMask(IofQData & iofqdata, AtomisticModel & model, Waters & waterModel, Mask & mask);

    void searchCOFactor(IofQData & iofqdata,
                       AtomisticModel & model,
                       Waters & waterModel,
                       float fixedCx,
                       float colo,
                       float cohi,
                       int totalCO);

    float getBestCx(){ return bestCx;} // excluded volume
    float getBestDW(){ return bestDW; }
    void printBestScores();

    float calculateDurbinWatson(float * residuals, unsigned int total);

    void chiFreeSearch(unsigned int totalRounds, IofQData &iofqdata, AtomisticModel &model, Waters &waterModel);

    void writeIcalcAtCXToFile(std::string name, unsigned int total, const std::vector<float> & qvalues, std::vector<float> & icalc,  std::vector<float> & diffA);

    void assemblePartials(AtomisticModel &model, Waters &waterModel,
        std::vector<float> & norm_Pamps,
        std::vector<float> & norm_Wamps,
        std::vector<float> & norm_Camps,
        std::vector<float> & aPW_cross_terms,
        std::vector<float> & aPX_cross_terms,
        std::vector<float> & aXW_cross_terms,
        unsigned int totalInWorkingSet);

    void assemblePartials(AtomisticModel &model, Waters &waterModel,
                              unsigned int totalInWorkingSet);

    void assemblePartialsWithMask(AtomisticModel &model,
                               Mask &mask,
                               unsigned int totalInWorkingSet);

    float calculateIofQfromModelWithVariances(unsigned int totalInWorkingSet,
                                          const float *pQ,
                                          const float *pAp,
                                          const float *pAc,
                                          const float *pAw,
                                          const float *pAXW_cross_term,
                                          const float *pAPX_cross_term,
                                          const float *pAPW_cross_term,
                                          float cx,
                                          float co,
                                          float *pIcalc);

    float calculateIofQfromModelWithVariancesMask(unsigned int totalInWorkingSet,
                                      const float *pQ,
                                      const float *pAp,
                                      const float *pAc,
                                      const float *pAXW_cross_term,
                                      const float *pAPX_cross_term,
                                      float cx,
                                      float co,
                                      float *pIcalc);

    float calculateIofQfromModelWithVariances(unsigned int totalInWorkingSet,
                                              const float *pQ,
                                              const float *pApw,
                                              const float *pAc,
                                              const float *pAXW_cross_term,
                                              float cx,
                                              float b_factor,
                                              float *pIcalc);

    float calculateIofQfromModelWithBackground(unsigned int totalInWorkingSet,
                                              const float *pQ,
                                              const float *pApw,
                                              const float *pAc,
                                              const float *pAXW_cross_term,
                                              float cx,
                                              float b_factor,
                                              float *pIcalc);


    void calculateIofQfromModel(unsigned int totalInWorkingSet,
                                           const float * pQ,
                                           const float * pAp,
                                           const float * pAw,
                                           const float * pAc,
                                           const float * pAPW_cross_term,
                                           const float * pAPX_cross_term,
                                           const float * pAXW_cross_term,
                                           float cx,
                                           float co,
                                           float * pIcalc
                                           );


    void writeModelToFile(float cx,
                          float bfactor,
                          std::vector<float> &qvalues,
                          AtomisticModel & model,
                          std::vector<float> &icalc,
                          const float * pAp,  // particle amplitude
                          const float * pAw,  // hydrated particle amplitude
                          const float * pAc,  // excluded volume
                          const float * pPW_cross_term, // cross-term
                          const float * pPX_cross_term, // cross-term
                          const float * pXW_cross_term // cross-term
                          );

    void simulate(float cx, float bfactor, std::vector<float> &qvalues, AtomisticModel &model, Waters &waterModel);

    void setUpForFitting(IofQData &iofqdata, AtomisticModel &model, Waters &waterModel);
    void setUpForFittingWithMask(IofQData &iofqdata, AtomisticModel &model, Mask &maskModel);

    void
    searchCX(IofQData &iofqdata, AtomisticModel &model, Waters &waterModel, float fixedBF, float lowerCX, float upperCX,
             int divisions);

    void fitFixedCoAndCx(IofQData &iofqdata, AtomisticModel &model, Waters &waterModel, float fixedCO, float fixedCX);


    void assemblePartialsFastHalf(AtomisticModel &model, Waters &waterModel, std::vector<float> &norm_aPs,
                                  std::vector<float> &norm_aPWs, std::vector<float> &norm_aWs,
                                  std::vector<float> &norm_aCs,
                                  std::vector<float> &aXW_cross_term, unsigned int totalInWorkingSet);

    void writeBestModelToFile(std::vector<float> &qvalues, std::string filename, AtomisticModel &model);

    void writeBestChiFreeModel(AtomisticModel &model,
                               const std::vector<unsigned int> &selectedIndices,
                               IofQData & iofqdata,
                               unsigned int totalRNDS);



    std::string getScoreText(std::string filename){

        char buffer [50];

        sprintf (buffer, "%8.3f %6.3f %5.3f %5.3f %5.1f %s\n", bestScore, bestChi,
                 bestDW, bestCx, bestOcc, filename.c_str()); // change the number to a string for name

        return buffer;
    }

    float
    calculateIofQfromModelWithVariancesMooreBfactor(unsigned int totalInWorkingSet, const float *pQ, const float *pApw,
                                                    const float *pNeighbors, const float *pNonNeighbors,
                                                    const float *pAXW_cross_term, float cx, float b_factor,
                                                    float *pIcalc);
};


#endif //WETSAXS_FIT_H
