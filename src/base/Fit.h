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

class Fit {

    unsigned int lmax, total_cx;
    float lowcx, highcx;
    float deltaCX, deltaCW;
    float upper_Bfactor;

    float bestScore;
    float bestBfactor, bestCx, bestScaleCoef, bestDW, bestChi;
    const float inv16Pi2 = 1.0/(16.0*M_PI_2);

    std::vector<float> norm_aPs;
    std::vector<float> norm_aPWs;
    std::vector<float> norm_aWs;
    std::vector<float> norm_aCs;
    std::vector<float> aXW_cross_term;

    std::vector<float> i_best;
    std::vector<float> i_obs;
    std::vector<float> i_obs_over_variance;
    std::vector<float> inv_variance;

    std::string header;

public:
    Fit(unsigned int lmax);
    Fit(unsigned int lmax, unsigned int division, float lowcx, float highcx, float upper_Bfactor);

    void search(IofQData & iofqdata, AtomisticModel & model, Waters & waterModel);
    void searchBFactor(IofQData & iofqdata,
                       AtomisticModel & model,
                       Waters & waterModel,
                       float fixedCx,
                       float upperB,
                       int totalB);

    float getBestCx(){ return bestCx;} // excluded volume
    float getBestBfactor(){ return bestBfactor; }
    float getBestDW(){ return bestDW; }
    void printBestScores();

    float calculateDurbinWatson(float * residuals, unsigned int total);

    void chiFreeSearch(unsigned int totalRounds, IofQData &iofqdata, AtomisticModel &model, Waters &waterModel);

    void writeIcalcAtCXToFile(std::string name, unsigned int total, const std::vector<float> & qvalues, std::vector<float> & icalc,  std::vector<float> & diffA);

    void assemblePartials(AtomisticModel &model, Waters &waterModel,
                          std::vector<float> &norm_aPs,
                          std::vector<float> &norm_aPWs,
                          std::vector<float> &norm_aWs,
                          std::vector<float> &norm_aCs,
                          std::vector<float> &aXW_cross_term,
                          unsigned int totalInWorkingSet);

    void assemblePartials(AtomisticModel &model, Waters &waterModel,
                              unsigned int totalInWorkingSet);

    float calculateIofQfromModelWithVariances(unsigned int totalInWorkingSet,
                                              const float *pQ,
                                              const float *pApw,
                                              const float *pAc,
                                              const float *pAXW_cross_term,
                                              float cx,
                                              float b_factor,
                                              float *pIcalc);

    void calculateIofQfromModel(unsigned int totalInWorkingSet,
                                const float *pQ, // q values (vector)
                                const float *pApw, // hydrated particle squared amplitude (vector)
                                const float *pAc,  // squared amplitude excluded volume (vector)
                                const float *pAXW_cross_term, // cross-term (vector)
                                 float cx, // contrast
                                 float b_factor, // b-factor
                                 float *pIcalc // output vector
                                 );

    void writeModelToFile(float cx,
                          float bfactor,
                          std::vector<float> &qvalues,
                          AtomisticModel & model,
                          std::vector<float> &icalc,
                          const float * aP,  // particle amplitude
                          const float * aC,  // excluded volume
                          const float * aCXP // cross-term
                          );

    void simulate(float cx, float bfactor, std::vector<float> &qvalues, AtomisticModel &model, Waters &waterModel);

    void setUpForFitting(IofQData &iofqdata, AtomisticModel &model, Waters &waterModel);

    void
    searchCX(IofQData &iofqdata, AtomisticModel &model, Waters &waterModel, float fixedBF, float lowerCX, float upperCX,
             int divisions);

    void
    fitFixedCxandBFactor(IofQData &iofqdata, AtomisticModel &model, Waters &waterModel, float fixedBF, float fixedCX);

    void assemblePartialsFastHalf(AtomisticModel &model, Waters &waterModel, std::vector<float> &norm_aPs,
                                  std::vector<float> &norm_aPWs, std::vector<float> &norm_aWs,
                                  std::vector<float> &norm_aCs,
                                  std::vector<float> &aXW_cross_term, unsigned int totalInWorkingSet);

    void
    writeBestModelToFile(std::vector<float> &qvalues, AtomisticModel &model, const float *aP,
                         const float *aC, const float *aCXP);

    void writeBestChiFreeModel(AtomisticModel &model,
                               const std::vector<unsigned int> &selectedIndices,
                               IofQData & iofqdata,
                               unsigned int totalRNDS);

    float * getAP() { return norm_aPWs.data(); }
};


#endif //WETSAXS_FIT_H
