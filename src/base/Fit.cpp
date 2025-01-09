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

#include "Fit.h"

Fit::Fit(unsigned int lmax) : lmax(lmax) {

}


Fit::Fit(unsigned int lmax, unsigned int division, float lowcx, float highcx, float upperB) :
        lmax(lmax),
        total_cx(division),
        lowcx(lowcx),
        highcx(highcx),
        upper_Bfactor(upperB) {

    deltaCX = (highcx-lowcx)/(float)division;
    bestScore = FLT_MAX;
    bestCx = lowcx;

}


void Fit::assemblePartials(AtomisticModel &model,
                           Waters &waterModel,
                           std::vector<float> & norm_Pamps,
                           std::vector<float> & norm_PWamps,
                           std::vector<float> & norm_Wamps,
                           std::vector<float> & norm_Camps,
                           std::vector<float> & XW_cross_term_amps,
                           unsigned int totalInWorkingSet){

    const float * almPR = model.getModelPartialsReal();
    const float * almPI = model.getModelPartialsImag();
    const float * xlmR = model.getExcludedVolumePartialsReal();
    const float * xlmI = model.getExcludedVolumePartialsImag();

    const float * wlmR = waterModel.getModelPartialsReal();
    const float * wlmI = waterModel.getModelPartialsImag();

    float * pAp = norm_Pamps.data();
    float * pApw = norm_PWamps.data();

    float * pAw = norm_Wamps.data();
    float * pAc = norm_Camps.data();

    // cross-terms
    float * pAXW_cross_term = XW_cross_term_amps.data();

    float a_p_norm, a_x_norm, a_w_norm, aPWR_term, aPWI_term, aXWR_term, aXWI_term;
    float x_w_r, x_w_i;
    float aPXR_term, aPXI_term;
    const float * pAlmPR, * pAlmPI, * pXlmR, * pXlmI, * pWlmR, *pWlmI;

    unsigned int totalpartialsperq = (lmax+1)*(lmax+1);
    // calculate squared norm terms
    for(unsigned int q =0; q < totalInWorkingSet; q++) {

        a_p_norm = 0.0;
        a_x_norm = 0.0;
        a_w_norm = 0.0;

        aPWR_term = 0.0;
        aPWI_term = 0.0;
        aXWR_term = 0.0;
        aXWI_term = 0.0;

        aPXR_term = 0.0;
        aPXI_term = 0.0;

        unsigned int startAt = q*totalpartialsperq; // always positive

        for(unsigned int lm=0; lm< totalpartialsperq; lm++){ // sums the L and M's per q value

            unsigned int index = startAt + lm;

            pAlmPR = &almPR[index];
            pAlmPI = &almPI[index];
            pXlmR = &xlmR[index];
            pXlmI = &xlmI[index];
            pWlmR = &wlmR[index];
            pWlmI = &wlmI[index];

            a_p_norm +=  *pAlmPR * *pAlmPR + *pAlmPI * *pAlmPI;     // sum of the A*Atranspose
            a_w_norm +=  *pWlmR * *pWlmR + *pWlmI * *pWlmI;         // water layer

            x_w_r = (*pXlmR + *pWlmR);
            x_w_i = (*pXlmI + *pWlmI);

            a_x_norm +=  x_w_r*x_w_r + x_w_i*x_w_i; // excluded volume, includes the water that form the scattered object

            // particle - water cross term
            aPWR_term += *pAlmPR * *pWlmR;
            aPWI_term += *pAlmPI * *pWlmI;

            // excluded volume and water cross-term
            aXWR_term += x_w_r * *pWlmR;
            aXWI_term += x_w_i * *pWlmI; // opaque body ?

            // excluded volume and particle cross-term
            aPXR_term += x_w_r * *pAlmPR; //
            aPXI_term += x_w_i * *pAlmPI; //
        }

        // PI constant gets absorbed into the scaling constant below during the chi2 fit
        pAp[q] = a_p_norm; //*pi_constant;
        pAw[q] = a_w_norm; //*pi_constant;
        pAc[q] = a_x_norm; //*pi_constant; // includes hydration waters

        pApw[q] = a_p_norm + a_w_norm + 2.0f*(aPWR_term + aPWI_term); // hydratedParticle

        pAXW_cross_term[q] = -2.0f*(aXWR_term + aXWI_term + aPXR_term + aPXI_term); // invert phase

    } // end loop over q values
}


void Fit::assemblePartials(AtomisticModel &model,
                           Waters &waterModel,
                           unsigned int totalInWorkingSet){

    const float * almPR = model.getModelPartialsReal();
    const float * almPI = model.getModelPartialsImag();
    const float * xlmR = model.getExcludedVolumePartialsReal();
    const float * xlmI = model.getExcludedVolumePartialsImag();

    const float * wlmR = waterModel.getModelPartialsReal();
    const float * wlmI = waterModel.getModelPartialsImag();

    float * pAp = norm_aPs.data();
    float * pApw = norm_aPWs.data();
    float * pAw = norm_aWs.data();
    float * pAc = norm_aCs.data();

    // cross-terms
    float * pAXW_cross_term = aXW_cross_term.data();

    float a_p_norm, a_x_norm, a_w_norm, aPWR_term, aPWI_term, aXWR_term, aXWI_term;
    float x_w_r, x_w_i;
    float aPXR_term, aPXI_term;
    const float * pAlmPR, * pAlmPI, * pXlmR, * pXlmI, * pWlmR, *pWlmI;

    unsigned int totalpartialsperq = (lmax+1)*(lmax+1);
    // calculate squared norm terms
    for(unsigned int q =0; q < totalInWorkingSet; q++) {

        a_p_norm = 0.0;
        a_x_norm = 0.0;
        a_w_norm = 0.0;

        aPWR_term = 0.0;
        aPWI_term = 0.0;
        aXWR_term = 0.0;
        aXWI_term = 0.0;

        aPXR_term = 0.0;
        aPXI_term = 0.0;

        unsigned int startAt = q*totalpartialsperq; // always positive

        for(unsigned int lm=0; lm< totalpartialsperq; lm++){ // sums the L and M's per q value

            unsigned int index = startAt + lm;

            pAlmPR = &almPR[index];
            pAlmPI = &almPI[index];
            pXlmR = &xlmR[index];
            pXlmI = &xlmI[index];
            pWlmR = &wlmR[index];
            pWlmI = &wlmI[index];

            a_p_norm +=  *pAlmPR * *pAlmPR + *pAlmPI * *pAlmPI;     // sum of the A*Atranspose
            a_w_norm +=  *pWlmR * *pWlmR + *pWlmI * *pWlmI;         // water layer

            x_w_r = (*pXlmR + *pWlmR);
            x_w_i = (*pXlmI + *pWlmI);

            a_x_norm +=  x_w_r*x_w_r + x_w_i*x_w_i; // excluded volume, includes the water that form the scattered object

            // particle - water cross term
            aPWR_term += *pAlmPR * *pWlmR;
            aPWI_term += *pAlmPI * *pWlmI;

            // excluded volume and water cross-term
            aXWR_term += x_w_r * *pWlmR;
            aXWI_term += x_w_i * *pWlmI; // opaque body ?

            // excluded volume and particle cross-term
            aPXR_term += x_w_r * *pAlmPR; //
            aPXI_term += x_w_i * *pAlmPI; //
        }

        // PI constant gets absorbed into the scaling constant below during the chi2 fit
        pAp[q] = a_p_norm; //*pi_constant;
        pAw[q] = a_w_norm; //*pi_constant;
        pAc[q] = a_x_norm; //*pi_constant;

        pApw[q] = a_p_norm + a_w_norm + 2.0f*(aPWR_term + aPWI_term); // hydratedParticle
        pAXW_cross_term[q] = -2.0f*(aXWR_term + aXWI_term + aPXR_term + aPXI_term); // invert phase
    } // end loop over q values
}





void Fit::assemblePartialsFastHalf(AtomisticModel &model,
                           Waters &waterModel,
                           std::vector<float> & norm_aPs,
                           std::vector<float> & norm_aPWs,
                           std::vector<float> & norm_aWs,
                           std::vector<float> & norm_aCs,
                           std::vector<float> & aXW_cross_term,
                           unsigned int totalInWorkingSet){

    const float * almPR = model.getModelPartialsReal();
    const float * almPI = model.getModelPartialsImag();
    const float * xlmR = model.getExcludedVolumePartialsReal();
    const float * xlmI = model.getExcludedVolumePartialsImag();

    const float * wlmR = waterModel.getModelPartialsReal();
    const float * wlmI = waterModel.getModelPartialsImag();

    float * pAp = norm_aPs.data();
    float * pAw = norm_aWs.data();
    float * pAc = norm_aCs.data();
    float * pApw = norm_aPWs.data();

    // cross-terms
    float * pAXW_cross_term = aXW_cross_term.data();

    float a_p_norm, a_x_norm, a_w_norm, aPWR_term, aPWI_term, aXWR_term, aXWI_term;
    float x_w_r, x_w_i;
    float aPXR_term, aPXI_term;
    const float * pAlmPR, * pAlmPI, * pXlmR, * pXlmI, * pWlmR, *pWlmI;

    unsigned int totalpartialsperq = (lmax+1)*(lmax+1);
    // calculate squared norm terms
    for(unsigned int q =0; q < totalInWorkingSet; q++) {

        a_p_norm = 0.0;
        a_x_norm = 0.0;
        a_w_norm = 0.0;

        aPWR_term = 0.0;
        aPWI_term = 0.0;
        aXWR_term = 0.0;
        aXWI_term = 0.0;

        aPXR_term = 0.0;
        aPXI_term = 0.0;

        unsigned int startAt = q*totalpartialsperq;

        for(unsigned int lm=0; lm< totalpartialsperq; lm++){ // sums the L and M's per q value

            unsigned int index = startAt + lm;

            pAlmPR = &almPR[index];
            pAlmPI = &almPI[index];
            pXlmR = &xlmR[index];
            pXlmI = &xlmI[index];
            pWlmR = &wlmR[index];
            pWlmI = &wlmI[index];

            a_p_norm +=  *pAlmPR * *pAlmPR + *pAlmPI * *pAlmPI;     // sum of the A*Atranspose
            a_w_norm +=  *pWlmR * *pWlmR + *pWlmI * *pWlmI;         // water layer

            x_w_r = (*pXlmR + *pWlmR);
            x_w_i = (*pXlmI + *pWlmI);

            a_x_norm +=  x_w_r*x_w_r + x_w_i*x_w_i; // excluded volume, includes the water that form the scattered object

            // particle - water cross term
            aPWR_term += *pAlmPR * *pWlmR;
            aPWI_term += *pAlmPI * *pWlmI;

            // excluded volume and water cross-term
            aXWR_term += x_w_r * *pWlmR;
            aXWI_term += x_w_i * *pWlmI; // opaque body ?

            // excluded volume and particle cross-term
            aPXR_term += x_w_r * *pAlmPR; //
            aPXI_term += x_w_i * *pAlmPI; //
        }

        // PI constant gets absorbed into the scaling constant below during the chi2 fit
        pAp[q] = a_p_norm; //*pi_constant;
        pAw[q] = a_w_norm; //*pi_constant;
        pAc[q] = a_x_norm; //*pi_constant;

        // can I get same values from above exploiting symmetry of spherical harmonics
        // total SH partials (lmax+1)*(lmax+1)*totalqvalues
        float sumItR = 0;
        float sumITaXWR_term = 0;
        for(int l=0; l<=lmax; l++){

            int base_index = (lmax+1)*(lmax+1)*q + l*l + l;
            // add up the m=0 term
            pAlmPR = &almPR[base_index];
            pAlmPI = &almPI[base_index];
            pWlmR = &wlmR[base_index];

            pXlmR = &xlmR[base_index];

            x_w_r = (*pXlmR + *pWlmR);

            sumItR += *pAlmPR * *pAlmPR + *pAlmPI * *pAlmPI;
            sumITaXWR_term += x_w_r * *pWlmR;


            for(int m=1; m<=l; m++){
                pAlmPR = &almPR[base_index + m];
                pAlmPI = &almPI[base_index + m];
                sumItR += 2 * (*pAlmPR * *pAlmPR + *pAlmPI * *pAlmPI);
            }
        }

        pApw[q] = a_p_norm + a_w_norm + 2.0f*(aPWR_term + aPWI_term); // hydratedParticle
        pAXW_cross_term[q] = -2.0f*(aXWR_term + aXWI_term + aPXR_term + aPXI_term); // invert phase
    } // end loop over q values
}

// set up the containers for holding complex norms for calculating intensities
void Fit::setUpForFitting(IofQData &iofqdata, AtomisticModel &model, Waters &waterModel) {

    unsigned int totalInWorkingSet = iofqdata.getTotalInWorkingSet();

    norm_aPs.resize(totalInWorkingSet);
    norm_aPWs.resize(totalInWorkingSet);
    norm_aWs.resize(totalInWorkingSet);
    norm_aCs.resize(totalInWorkingSet);
    aXW_cross_term.resize(totalInWorkingSet);

    assemblePartials(model, waterModel, totalInWorkingSet);

    // populate experimental observation information
    //std::vector<Datum> & workingSet = const_cast<std::vector<Datum> &>(iofqdata.getWorkingSet());
    const std::vector<Datum> & workingSet = iofqdata.getWorkingSet();

    i_obs_over_variance.resize(totalInWorkingSet);
    inv_variance.resize(totalInWorkingSet);
    i_obs.resize(totalInWorkingSet);
    i_best.resize(totalInWorkingSet);

    for(unsigned int q =0; q < totalInWorkingSet; q++) {
        const Datum *pW = &workingSet[q];
        i_obs[q] = pW->getI();
        i_obs_over_variance[q] = pW->getI() * pW->getInvVar();
        inv_variance[q] = pW->getInvVar();
    }
}


void Fit::search(IofQData &iofqdata, AtomisticModel &model, Waters &waterModel) {

    setUpForFitting(iofqdata, model, waterModel); // takes < 0.001 seconds

    float * pAc = norm_aCs.data();
    float * pApw = norm_aPWs.data();
    // cross-terms
    float * pAXW_cross_term = aXW_cross_term.data();

//    const float * excludedVolumeNeighbors = model.getExcludedVolumeIntensitiesNeighbors();
//    const float * excludedVolumeNonNeighbors = model.getExcludedVolumeIntensitiesNonNeighbors();

    // perform grid search
    float sum, diff, scalar;
    bestCx = 0.0f;
    bestBfactor = 0.0f;
    bestScore = FLT_MAX;

    auto & qvalues = iofqdata.getWorkingSetQvalues();
    auto pQ = qvalues.data();

    unsigned int totalInWorkingSet = qvalues.size();
    std::vector<float> residuals(totalInWorkingSet);
    std::vector<float> i_calc(totalInWorkingSet);
    auto pIcalc = i_calc.data();

    float delta_b = 0.331f; // actual equation is exp(-q^2*d_w/3), d_w goes to 100, divided by 3 is 33

    //if (excludedVolumeNonNeighbors[totalInWorkingSet-1] < 0){

//        float neighbor, nonneighbor;
//        float maxdiff = FLT_MAX;
//        float top = 0;
//        float bottom = 1;
//        float q2=1;
//
//        for (unsigned int i=0; i<totalInWorkingSet; i++){
//
//             neighbor = excludedVolumeNeighbors[i];
//             nonneighbor = excludedVolumeNonNeighbors[i];
//
//             if (nonneighbor < 0){
//              float mdiff = -nonneighbor/neighbor;
//                 if (mdiff < maxdiff ){
//                     maxdiff = mdiff;
//                     q2 = qvalues[i];
//                 }
//             }
//        }
//
//        q2 *= q2;
//        std::cout << " q2 :::::  " << q2 << " " << maxdiff << std::endl;
//        upper_Bfactor = 0.99*sqrtf(-std::logf(abs(maxdiff))*2/q2);
//        delta_b = upper_Bfactor/33;

    //}

    std::cout << " upper B :::::  " << upper_Bfactor << std::endl;

    float chi, b_factor, dw;
    auto total_b = (int)std::ceil(upper_Bfactor/delta_b); // divided by 3

    // setting match point requires total waters in hydration model
    model.setMatchPointScale(waterModel.getTotalWaters());
    float contrastMatchScale = model.getMatchPoint();
    SASTOOLS_UTILS_H::logger("Match Point", formatNumber(model.getMatchPoint(),3));

    float cx =  lowcx*contrastMatchScale;
    float inv_n_s = 1.0f/(float)(totalInWorkingSet-2);

    for(int incr_cx = 0; incr_cx < total_cx; incr_cx++){

        b_factor = 0.0f;

        for (int incr_b = 0; incr_b < total_b; incr_b++){

//            scalar = calculateIofQfromModelWithBackground(totalInWorkingSet,
//                                                         pQ,
//                                                         pApw,
//                                                         pAc,
//                                                         pAXW_cross_term,
//                                                         cx,
//                                                         b_factor,
//                                                         pIcalc
//            );

            scalar = calculateIofQfromModelWithVariances(totalInWorkingSet,
                                                         pQ,
                                                         pApw,
                                                         pAc,
                                                         pAXW_cross_term,
                                                         cx,
                                                         b_factor,
                                                         pIcalc
            );


//            scalar = calculateIofQfromModelWithVariancesNewMethod(totalInWorkingSet,
//                                                         pQ,
//                                                         pApw,
//                                                         pAc,
//                                                         pAXW_cross_term,
//                                                         cx,
//                                                         b_factor,
//                                                         pIcalc
//            );

//            scalar = calculateIofQfromModelWithVariancesMooreBfactor(totalInWorkingSet,
//                    pQ,
//                    pApw,
//                    excludedVolumeNeighbors,
//                    excludedVolumeNonNeighbors,
//                    pAXW_cross_term,
//                    cx,
//                    b_factor,
//                    pIcalc
//                    );

            // score the model
            sum = 0.0f;
            for(unsigned int q =0; q < totalInWorkingSet; q++) {
                diff = i_obs[q] - scalar*pIcalc[q];
                sum += diff*diff*inv_variance[q];
                residuals[q] = diff;
            }

            sum *= inv_n_s;
            chi = sum;
            dw = fabs(2 - calculateDurbinWatson(residuals.data(), totalInWorkingSet)); // should be near zero
            sum += dw;

            if (sum < bestScore){
                bestScore = sum;
                bestChi = chi;
                bestBfactor = b_factor;
                bestCx = cx;
                bestScaleCoef = scalar;
                bestDW = dw;

                // do point transfer and do the scalar after finishing
                for(unsigned int q =0; q < totalInWorkingSet; q++) {
                    i_best[q] = scalar*pIcalc[q];
                }
            }

            b_factor += delta_b;
        }

        cx += deltaCX*contrastMatchScale;
    }

    SASTOOLS_UTILS_H::logger("CX RANGE", formatNumber(lowcx*contrastMatchScale,3) + " to " + formatNumber(cx,3));

    std::cout << "_______________________________________________________________________" << std::endl;

    printBestScores();

//    header = iofqdata.getHeader();
//    writeBestModelToFile(const_cast<std::vector<float> &>(qvalues), model, pApw, pAc, pAXW_cross_term);

}

void Fit::printBestScores(){
    SASTOOLS_UTILS_H::logger("Best Score", formatNumber(bestScore,3));
    SASTOOLS_UTILS_H::logger("chi-squared", formatNumber(bestChi,3));
    SASTOOLS_UTILS_H::logger("Durbin-Watson", formatNumber(bestDW,3));
    SASTOOLS_UTILS_H::logger("B-factor", formatNumber(bestBfactor,1));
    SASTOOLS_UTILS_H::logger("contrast(c_x)", formatNumber(bestCx,2));
}


void Fit::searchBFactor(IofQData &iofqdata, AtomisticModel &model, Waters &waterModel, float fixedCx, float upperB, int totalB) {

    setUpForFitting(iofqdata, model, waterModel);

    upper_Bfactor = totalB;
    const float delta_b = upper_Bfactor/(float)totalB; // actual equation is exp(-q^2*d_w/3), d_w goes to 100, divided by 3 is 33
    unsigned int totalInWorkingSet = iofqdata.getTotalInWorkingSet();

    float * pAc = norm_aCs.data();
    float * pApw = norm_aPWs.data();
    // cross-terms
    float * pAXW_cross_term = aXW_cross_term.data();

    // perform grid search
    float sum, diff, scalar;
    float chi, dw;

    bestCx = 0.0f;
    bestBfactor = 0.0f;
    bestScore = FLT_MAX;

    auto & qvalues = iofqdata.getWorkingSetQvalues();
    auto pQ = qvalues.data();

    std::vector<float> residuals(totalInWorkingSet);
    std::vector<float> i_calc(totalInWorkingSet);
    auto pIcalc = i_calc.data();

    // setting match point requires total waters in hydration model
    model.setMatchPointScale(waterModel.getTotalWaters());
    float contrastMatchScale = model.getMatchPoint();
    SASTOOLS_UTILS_H::logger("Match Point", formatNumber(model.getMatchPoint(),3));

    float cx =  fixedCx*contrastMatchScale;
    bestCx = cx;

    float inv_n_s = 1.0f/(float)(totalInWorkingSet-2);

    float b_factor = 0.0f;

    for (int incr_b = 0; incr_b < totalB; incr_b++){

        scalar = calculateIofQfromModelWithVariances(totalInWorkingSet,
                                                     pQ,
                                                     pApw,
                                                     pAc,
                                                     pAXW_cross_term,
                                                     cx,
                                                     b_factor,
                                                     pIcalc
        );

        // score the model
        sum = 0.0f;
        for(unsigned int q =0; q < totalInWorkingSet; q++) {
            diff = i_obs[q] - scalar * pIcalc[q];
            sum += diff*diff*inv_variance[q];
            residuals[q] = diff;
        }

        sum *= inv_n_s;
        chi = sum;
        dw = fabs(2 - calculateDurbinWatson(residuals.data(), totalInWorkingSet));
        sum += dw;

        if (sum < bestScore){
            bestScore = sum;
            bestChi = chi;
            bestBfactor = b_factor;
            bestScaleCoef = scalar;
            bestDW = dw;

            for(unsigned int q =0; q < totalInWorkingSet; q++) {
                i_best[q] = scalar * pIcalc[q];
            }
        }

        b_factor += delta_b;
    }
}


void Fit::searchCX(IofQData &iofqdata, AtomisticModel &model, Waters &waterModel, float fixedBF, float lowerCX, float upperCX, int divisions) {

    lowcx = lowerCX;
    highcx = upperCX;

    setUpForFitting(iofqdata, model, waterModel);

    deltaCX = (highcx-lowcx)/(float)divisions;

    unsigned int totalInWorkingSet = iofqdata.getTotalInWorkingSet();

    float * pAc = norm_aCs.data();
    float * pApw = norm_aPWs.data();
    // cross-terms
    float * pAXW_cross_term = aXW_cross_term.data();

    // perform grid search
    float sum, diff, scalar;
    float chi, dw;

    bestCx = 0.0f;
    bestBfactor = fixedBF;
    bestScore = FLT_MAX;

    auto & qvalues = iofqdata.getWorkingSetQvalues();
    auto pQ = qvalues.data();

    std::vector<float> residuals(totalInWorkingSet);
    std::vector<float> i_calc(totalInWorkingSet);
    auto pIcalc = i_calc.data();

    // setting match point requires total waters in hydration model
    model.setMatchPointScale(waterModel.getTotalWaters());
    float contrastMatchScale = model.getMatchPoint();
    SASTOOLS_UTILS_H::logger("Match Point", formatNumber(model.getMatchPoint(),3));

    float inv_n_s = 1.0f/(float)(totalInWorkingSet-2);

    float cx =  lowcx*contrastMatchScale; // fractions of the contrastMatchScale
    for(int incr_cx = 0; incr_cx < total_cx; incr_cx++){

        scalar = calculateIofQfromModelWithVariances(totalInWorkingSet,
                                                     pQ,
                                                     pApw,
                                                     pAc,
                                                     pAXW_cross_term,
                                                     cx,
                                                     fixedBF,
                                                     pIcalc
        );

        // score the model
        sum = 0.0f;
        for(unsigned int q =0; q < totalInWorkingSet; q++) {
            diff = i_obs[q] - scalar * pIcalc[q];
            sum += diff*diff*inv_variance[q];
            residuals[q] = diff;
        }

        sum *= inv_n_s;
        chi = sum;
        dw = fabs(2 - calculateDurbinWatson(residuals.data(), totalInWorkingSet));
        sum += dw;

        if (sum < bestScore){
            bestScore = sum;
            bestChi = chi;
            bestCx = cx;
            bestScaleCoef = scalar;
            bestDW = dw;

            for(unsigned int q =0; q < totalInWorkingSet; q++) {
                i_best[q] = scalar * pIcalc[q];
            }
        }

        cx += deltaCX*contrastMatchScale;
    }
}


// Cx returned by fitting is already scaled to matchpoint
void Fit::fitFixedCxandBFactor(IofQData &iofqdata,
                               AtomisticModel &model,
                               Waters &waterModel,
                               float fixedBF,
                               float fixedCX) {

    setUpForFitting(iofqdata, model, waterModel);

    float * pAc = norm_aCs.data();
    float * pApw = norm_aPWs.data();

//    const float * excludedVolumeNeighbors = model.getExcludedVolumeIntensitiesNeighbors();
//    const float * excludedVolumeNonNeighbors = model.getExcludedVolumeIntensitiesNonNeighbors();

    // cross-terms
    float * pAXW_cross_term = aXW_cross_term.data();

    auto & qvalues = iofqdata.getWorkingSetQvalues();
    auto pQ = qvalues.data();

    unsigned int totalInWorkingSet = iofqdata.getTotalInWorkingSet();
    std::vector<float> residuals(totalInWorkingSet);
    std::vector<float> i_calc(totalInWorkingSet);
    auto pIcalc = i_calc.data();

    // setting match point requires total waters in hydration model
    // model.setMatchPointScale(waterModel.getTotalWaters());
    // SASTOOLS_UTILS_H::logger("Match Point", formatNumber(model.getMatchPoint(),3));

    float inv_n_s = 1.0f/(float)(totalInWorkingSet-2);

    float cx =  fixedCX;//*model.getMatchPoint(); // fractions of the contrastMatchScale

    float scalar = calculateIofQfromModelWithVariances(totalInWorkingSet,
                                                 pQ,
                                                 pApw,
                                                 pAc,
                                                 pAXW_cross_term,
                                                 cx,
                                                 fixedBF,
                                                 pIcalc
    );

//    float scalar = calculateIofQfromModelWithVariancesMooreBfactor(totalInWorkingSet,
//                                                             pQ,
//                                                             pApw,
//                                                             excludedVolumeNeighbors,
//                                                             excludedVolumeNonNeighbors,
//                                                             pAXW_cross_term,
//                                                             cx,
//                                                             fixedBF,
//                                                             pIcalc
//    );

    // score the model
    float sum = 0.0f, diff;
    for(unsigned int q =0; q < totalInWorkingSet; q++) {
        diff = i_obs[q] - scalar * pIcalc[q];
        sum += diff*diff*inv_variance[q];
        residuals[q] = diff;
    }

    sum *= inv_n_s;
    float chi = sum;
    float dw = fabs(2 - calculateDurbinWatson(residuals.data(), totalInWorkingSet));
    sum += dw;

    bestScore = sum;
    bestChi = chi;
    bestCx = fixedCX;
    bestBfactor = fixedBF;
    bestScaleCoef = scalar;
    bestDW = dw;

    SASTOOLS_UTILS_H::logger("Overall BestScore", formatNumber(bestScore, 2));
    SASTOOLS_UTILS_H::logger("Overall CHI2", formatNumber(bestChi, 3));
    SASTOOLS_UTILS_H::logger("Overall DW", formatNumber(bestDW, 2));

    for(unsigned int q =0; q < totalInWorkingSet; q++) {
        i_best[q] = scalar * pIcalc[q];
    }
}

/*
 * use this for calculating models with no input q-range
 * simulated waters, must specify contrast and bfactor
 *
 */
void Fit::writeModelToFile(float cx,
                           float bfactor,
                           std::vector<float> & qvalues,
                           AtomisticModel & model,
                           std::vector<float> & icalc,
                           const float * aP,  // particle amplitude
                           const float * aC,  // excluded volume
                           const float * aCXP // cross-term
                           ){

    std::string name = model.getPDBModel().getFileStemName();

    char buffer [50];

    sprintf (buffer, "B_%.2f", bfactor); // change the number to a string for name
    std::string sbf = buffer;
    std::replace(sbf.begin(), sbf.end(), '.', 'p');

    sprintf (buffer, "CX_%.2f", cx); // change the number to a string for name
    std::string cxs = buffer;
    std::replace(cxs.begin(), cxs.end(), '.', 'p');

    std::string nameOf = name + "_" + cxs + "_" + sbf + "_calc.dat";

    FILE * pFile = fopen(nameOf.c_str(), "w");
    unsigned int total = qvalues.size();

    fprintf(pFile, "# REMARK           :: %s\n", model.getPDBModel().getFilename().c_str());
    fprintf(pFile, "# REMARK        CX :: %.4f \n", cx);
    fprintf(pFile, "# REMARK         B :: %.2f \n", bfactor);
    fprintf(pFile, "# REMARK      QMIN :: %.6f \n", qvalues[0]);
    fprintf(pFile, "# REMARK      QMAX :: %.6f \n", qvalues[total-1]);
    fprintf(pFile, "# REMARK      LMAX :: %i \n", lmax);
    fprintf(pFile, "# REMARK      DMAX :: %.1f \n", model.getDmax());
    fprintf(pFile, "# REMARK  RESIDUES :: %i \n", model.getPDBModel().getTotalResidues());
    fprintf(pFile, "# REMARK     ATOMS :: %i \n", model.getPDBModel().getTotalCoordinates());
    fprintf(pFile, "# REMARK COLUMNS 1 :: %s \n", "index");
    fprintf(pFile, "# REMARK COLUMNS 2 :: %s \n", "momentum transfer (inverse Angstroms");
    fprintf(pFile, "# REMARK COLUMNS 3 :: %s \n", "Icalc");
    fprintf(pFile, "# REMARK COLUMNS 4 :: %s \n", "hydrated particle");
    fprintf(pFile, "# REMARK COLUMNS 5 :: %s \n", "excluded volume");
    fprintf(pFile, "# REMARK COLUMNS 6 :: %s \n", "cross-terms");


    float q_val, exp_term;

    for(unsigned int i=0; i < total; i++){
        q_val = qvalues[i];
        exp_term = expf(-(q_val*q_val)*bfactor*inv16Pi2)*cx;

        fprintf(pFile, "%5i %.3E %.5E %.5E %.5E %.5E\n", (i+1), qvalues[i], icalc[i], aP[i], aC[i]*exp_term*exp_term, aCXP[i]);
    }

    fclose(pFile);
}


void Fit::writeBestModelToFile(
                           std::vector<float> & qvalues,
                           std::string filename,
                           AtomisticModel & model,
                           const float * aP,  // particle amplitude
                           const float * aC,  // excluded volume
                           const float * aCXP // cross-term
){

    std::string name = model.getPDBModel().getFileStemName();

    char buffer [50];

//    sprintf (buffer, "B_%.2f", bestBfactor); // change the number to a string for name
//    std::string sbf = buffer;
//    std::replace(sbf.begin(), sbf.end(), '.', 'p');
//
//    sprintf (buffer, "CX_%.2f", bestCx); // change the number to a string for name
//    std::string cxs = buffer;
//    std::replace(cxs.begin(), cxs.end(), '.', 'p');
//
//    std::string nameOf = name + "_" + cxs + "_" + sbf + "_fit.dat";
    std::string nameOf = name + "_best_fit.dat";

    logger("OUTPUT FILE", nameOf);

    FILE * pFile = fopen(nameOf.c_str(), "w");
    unsigned int total = qvalues.size();

    // set header
    if (!header.empty()){
        fprintf(pFile,"%s", header.c_str());
    }

    fprintf(pFile, "# REMARK              :: %s\n", model.getPDBModel().getFilename().c_str());
    fprintf(pFile, "# REMARK         DATA :: %s \n", filename.c_str());
    fprintf(pFile, "# REMARK           CX :: %.4f \n", bestCx);
    fprintf(pFile, "# REMARK            B :: %.2f \n", bestBfactor);
    fprintf(pFile, "# REMARK        SCORE :: %.3f \n", bestScore);
    fprintf(pFile, "# REMARK        CHI^2 :: %.3f \n", bestChi);
    fprintf(pFile, "# REMARK           DW :: %.3f \n", bestDW);
    fprintf(pFile, "# REMARK       scalar :: %.2f \n", bestScaleCoef);
    fprintf(pFile, "# REMARK         QMIN :: %.6f \n", qvalues[0]);
    fprintf(pFile, "# REMARK         QMAX :: %.6f \n", qvalues[total-1]);
    fprintf(pFile, "# REMARK         LMAX :: %i \n", lmax);
    fprintf(pFile, "# REMARK         DMAX :: %.1f Angstroms\n", model.getDmax());
    fprintf(pFile, "# REMARK     RESIDUES :: %i \n", model.getPDBModel().getTotalResidues());
    fprintf(pFile, "# REMARK        ATOMS :: %i \n", model.getPDBModel().getTotalCoordinates());
    fprintf(pFile, "# REMARK (non-H) MASS :: %.1f \n", model.getPDBModel().getMW());
    fprintf(pFile, "# REMARK       VOLUME :: %.1f \n", model.getPDBModel().getVolume());
    fprintf(pFile, "# REMARK    COLUMNS 1 :: %s \n", "index");
    fprintf(pFile, "# REMARK    COLUMNS 2 :: %s \n", "momentum transfer (Angstroms^-1)");
    fprintf(pFile, "# REMARK    COLUMNS 3 :: %s \n", "I_observed");
    fprintf(pFile, "# REMARK    COLUMNS 4 :: %s \n", "I_model");
    fprintf(pFile, "# REMARK    COLUMNS 5 :: %s \n", "residual [I_obs - I_model]");
    fprintf(pFile, "# REMARK    COLUMNS 6 :: %s \n", "Particle Scattering");
    fprintf(pFile, "# REMARK    COLUMNS 7 :: %s \n", "Excluded Volume Scattering");
    fprintf(pFile, "# REMARK    COLUMNS 8 :: %s \n", "Cross-term Scattering");

    for(unsigned int i=0; i < total; i++){
        fprintf(pFile, "%5i %.3E %.5E %.5E %.5E %.5E %.5E %.5E\n", (i+1), qvalues[i], i_obs[i], i_best[i], (i_obs[i] - i_best[i]), aP[i], aC[i], aCXP[i]);
    }

    fclose(pFile);
}

/*
 * returns scale factor
 */
float Fit::calculateIofQfromModelWithVariances(unsigned int totalInWorkingSet,
                                               const float * pQ,
                                               const float * pApw,
                                               const float * pAc,
                                               const float * pAXW_cross_term,
                                               float cx,
                                               float b_factor,
                                               float * pIcalc
                                 ){
    float i_of_q, q_val;
    float exp_term;

    float c_obs_calc = 0.0f;
    float calc_var = 0.0f;

    for(unsigned int q =0; q < totalInWorkingSet; q++) {

        q_val = pQ[q]; // get qvalue to calculate b-factor
        exp_term = expf(-(q_val*q_val)*b_factor*inv16Pi2)*cx;

        i_of_q = pApw[q] + pAc[q]*exp_term*exp_term + exp_term*pAXW_cross_term[q];

        c_obs_calc += i_of_q*i_obs_over_variance[q];
        calc_var += i_of_q*i_of_q*inv_variance[q];

        pIcalc[q] = i_of_q;
    }

    return c_obs_calc/calc_var;
}


/*
 * returns scale factor
 */
float Fit::calculateIofQfromModelWithVariancesMooreBfactor(unsigned int totalInWorkingSet,
                                               const float * pQ,
                                               const float * pApw,
                                               const float * pNeighbors,
                                               const float * pNonNeighbors,
                                               const float * pAXW_cross_term,
                                               float cx,
                                               float b_factor,
                                               float * pIcalc
                                               ){

    float i_of_q, q_val;
    float exp_term;
    float constant = -0.5*b_factor;

    float c_obs_calc = 0.0f;
    float calc_var = 0.0f;

    for(unsigned int q =0; q < totalInWorkingSet; q++) {

        q_val = pQ[q]; // get qvalue to calculate b-factor
        exp_term = expf(constant*q_val*q_val)*cx;
        //exp_term = expf(-(q_val*q_val)*b_factor*inv16Pi2)*cx;

        //i_of_q = pApw[q] + pAc[q]*exp_term*exp_term + exp_term*pAXW_cross_term[q];
        //i_of_q = pApw[q] + pAc[q]*exp_term*exp_term + exp_term*pAXW_cross_term[q];
        i_of_q = pApw[q] + pNonNeighbors[q]*cx*cx + pNeighbors[q]*exp_term*exp_term + pAXW_cross_term[q];
        //i_of_q = pApw[q] + pNonNeighbors[q]*exp_term*exp_term + pNeighbors[q]*exp_term*exp_term + exp_term*pAXW_cross_term[q];

        // std::cout << q << " " << i_of_q << std::endl;
        c_obs_calc += i_of_q*i_obs_over_variance[q];
        calc_var += i_of_q*i_of_q*inv_variance[q];

        pIcalc[q] = i_of_q;
    }

    return c_obs_calc/calc_var;

}


/*
 * returns scale factor
 */
float Fit::calculateIofQfromModelWithBackground(unsigned int totalInWorkingSet,
                                               const float * pQ,
                                               const float * pApw,
                                               const float * pAc,
                                               const float * pAXW_cross_term,
                                               float cx,
                                               float b_factor,
                                               float * pIcalc
){
    float i_of_q, q_val;
    float exp_term;
    float y_val;

    std::vector<float> calculated(totalInWorkingSet);

    float sumxy=0, sumy = 0, sumx = 0, sumx2 = 0;

    for(unsigned int q =0; q < totalInWorkingSet; q++) {

        q_val = pQ[q]; // get qvalue to calculate b-factor
        exp_term = expf(-(q_val*q_val)*b_factor*inv16Pi2)*cx;

        i_of_q = pApw[q] + pAc[q]*exp_term*exp_term + exp_term*pAXW_cross_term[q];

        calculated[q] = i_of_q;
        // std::cout << q << " " << i_of_q << std::endl;
        y_val = i_obs[q];

        sumy += y_val;
        sumxy += i_of_q*y_val;
        sumx += i_of_q;
        sumx2 += i_of_q*i_of_q;
    }

    float scale = (totalInWorkingSet*sumxy - sumy*sumx)/(totalInWorkingSet*sumx2 - sumx*sumx);
    float background = (sumy - scale*sumx)/totalInWorkingSet;

    for(unsigned int q =0; q < totalInWorkingSet; q++) {
        pIcalc[q] = scale*calculated[q] + background;
    }

    return scale;
}


void Fit::calculateIofQfromModel(unsigned int totalInWorkingSet,
                                               const float * pQ,
                                               const float * pApw,
                                               const float * pAc,
                                               const float * pAXW_cross_term,
                                               float cx,
                                               float b_factor,
                                               float * pIcalc
                                               ){

    float q_val;
    float exp_term;

    for(unsigned int q =0; q < totalInWorkingSet; q++) {
        q_val = pQ[q];
        exp_term = expf(-(q_val*q_val)*b_factor*inv16Pi2)*cx;
        pIcalc[q] = pApw[q] + pAc[q]*exp_term*exp_term + exp_term*pAXW_cross_term[q];
    }
}

void Fit::simulate(float cx,
                   float bfactor,
                   std::vector<float> & qvalues,
                   AtomisticModel &model,
                   Waters &waterModel){


    unsigned int totalInWorkingSet = qvalues.size();

    std::vector<float> norm_asPs(totalInWorkingSet);
    std::vector<float> norm_asPWs(totalInWorkingSet);
    std::vector<float> norm_asWs(totalInWorkingSet);
    std::vector<float> norm_asCs(totalInWorkingSet);
    std::vector<float> asXW_cross_term(totalInWorkingSet);

    assemblePartials(model,
                     waterModel,
                     norm_asPs,
                     norm_asPWs,
                     norm_asWs,
                     norm_asCs,
                     asXW_cross_term,
                     totalInWorkingSet);

    float * pAc = norm_asCs.data();
    float * pApw = norm_asPWs.data();
    float * pAXW_cross_term = asXW_cross_term.data();

    std::vector<float> i_calc(totalInWorkingSet);
    auto pIcalc = i_calc.data();

    calculateIofQfromModel(totalInWorkingSet,
                                        qvalues.data(),
                                        pApw,
                                        pAc,
                                        pAXW_cross_term,
                                        cx,
                                        bfactor,
                                        pIcalc
    );

    writeModelToFile(cx, bfactor, qvalues, model, i_calc, pApw, pAc, pAXW_cross_term);
}



float Fit::calculateDurbinWatson(float *residuals, unsigned int total) {

    float val, diff, first = residuals[0];
    float squared_sum = first*first;
    float squared_diff = 0.0;

    for(unsigned int i=1; i<total; i++){
        val = residuals[i];
        diff = (val - residuals[i-1]);
        squared_diff += diff*diff;
        squared_sum += val*val;
    }

    return squared_diff/squared_sum;
}


// need a CV Chi-square test using fitted parameters from above
void Fit::chiFreeSearch(unsigned int totalRounds, IofQData &iofqdata, AtomisticModel &model, Waters &waterModel){

    std::vector<Result> results;

    if (!(totalRounds & 1)){
        totalRounds+=1;
    }

    for (unsigned int rnd=0; rnd<totalRounds; rnd++){

        SASTOOLS_UTILS_H::logger("CHI-free round", formatNumber(rnd+1));

        iofqdata.makeWorkingSet();
        std::vector<float> qvalues = iofqdata.getWorkingSetQvalues();

        waterModel.calculatePartialAmplitudes(lmax, iofqdata.getTotalInWorkingSet(), qvalues, (rnd > 1));
        model.calculatePartialAmplitudes(lmax, iofqdata.getTotalInWorkingSet(), qvalues, (rnd > 1));

        this->search(iofqdata, model, waterModel);
        results.emplace_back(Result(bestScore, bestScaleCoef, bestCx, bestBfactor, bestDW, bestChi, i_best, iofqdata.getSelectedIndices()));
    }

    // now sort and take the median of Results
    // Using lambda expressions in C++11
    std::sort(results.begin(), results.end(), [](const Result& lhs, const Result& rhs) {
        return lhs.score < rhs.score;
    });

    std::cout << "_______________________________________________________________________" << std::endl;
    std::cout << " SORTED SCORES " << std::endl;
    std::cout << "   SCORE ::   CHI   ::   DW  ::   CX  ::    B" << std::endl;

    char buffer [50];
    for(auto & res : results){
        sprintf (buffer, "%8.3f    %6.3f     %5.3f    %5.3f    %5.1f", res.score, res.getChi(), res.getDurbinWatson(), res.getCX(), res.getBfactor()); // change the number to a string for name
        std::cout << buffer << std::endl;
    }

    std::cout << "_______________________________________________________________________" << std::endl;
    // median
    int middle = (totalRounds-1)/2;

    auto & res = results[middle];
    bestScore = res.score;
    bestChi - res.getChi();
    bestDW = res.getDurbinWatson();
    bestCx = res.getCX();
    bestBfactor = res.getBfactor();
    i_best = res.getICalc();

    SASTOOLS_UTILS_H::logger("STATUS", "Writing median best-fit model");
    writeBestChiFreeModel(model, res.getSelectedIndices(), iofqdata, totalRounds);

    // calculate over all input data using best fitting median parameters
    SASTOOLS_UTILS_H::logger("STATUS", "FINAL FIT TO ALL DATA");

    iofqdata.setAllDataToWorkingSet();

    std::vector<float> qvalues = iofqdata.getWorkingSetQvalues();

//    // recalculate Partials for the initial models
    waterModel.calculatePartialAmplitudes(lmax, qvalues.size(), qvalues, true);
    model.calculatePartialAmplitudes(lmax, qvalues.size(), qvalues, true);
//
    fitFixedCxandBFactor(iofqdata, model, waterModel, bestBfactor, bestCx);
    writeBestModelToFile(qvalues, iofqdata.getFilename(), model, norm_aPWs.data(), norm_aCs.data(), aXW_cross_term.data());

}

//buffer,  qvalues, i_calc, diff
void Fit::writeIcalcAtCXToFile(std::string name, unsigned int total, const std::vector<float> & qvalues, std::vector<float> & icalc, std::vector<float> & diffs) {

    std::string nameOf = "cx_" + name+"_.txt";
    FILE * pFile = fopen(nameOf.c_str(), "w");

    fprintf(pFile, "CX: %s\n", name.c_str());
    for(unsigned int i=0; i < total; i++){
        fprintf(pFile, "%5i %.3E %.5E %.5E \n", i, qvalues[i], icalc[i], diffs[i]);
    }

    fclose(pFile);
}


void Fit::writeBestChiFreeModel(AtomisticModel &model,
                                const std::vector<unsigned int> &selectedIndices,
                                IofQData & iofqdata, unsigned int totalRNDS) {

    std::string name = model.getPDBModel().getFileStemName();

    char buffer [50];

    sprintf (buffer, "B_%.2f", bestBfactor); // change the number to a string for name
    std::string sbf = buffer;
    std::replace(sbf.begin(), sbf.end(), '.', 'p');

    sprintf (buffer, "CX_%.2f", bestCx); // change the number to a string for name
    std::string cxs = buffer;
    std::replace(cxs.begin(), cxs.end(), '.', 'p');

    std::string nameOf = name + "_" + cxs + "_" + sbf + "_median_fit.dat";

    logger("OUTPUT FILE", nameOf);

    FILE * pFile = fopen(nameOf.c_str(), "w");

    std::vector<float> qvalues;
    for(auto & ind : selectedIndices){
        qvalues.push_back(iofqdata.getXDataAt(ind));
    }
    unsigned int total = qvalues.size();

    header = iofqdata.getHeader();
    if (!header.empty()){
        fprintf(pFile,"%s", header.c_str());
    }

    fprintf(pFile, "# REMARK              :: %s\n", model.getPDBModel().getFilename().c_str());
    fprintf(pFile, "# REMARK     CHI-FREE :: BEST MEDIAN\n");
    fprintf(pFile, "# REMARK   TOTAL RNDS :: %i\n", totalRNDS);
    fprintf(pFile, "# REMARK           CX :: %.4f \n", bestCx);
    fprintf(pFile, "# REMARK            B :: %.2f \n", bestBfactor);
    fprintf(pFile, "# REMARK        SCORE :: %.3f \n", bestScore);
    fprintf(pFile, "# REMARK        CHI^2 :: %.3f \n", bestChi);
    fprintf(pFile, "# REMARK           DW :: %.3f \n", bestDW);
    fprintf(pFile, "# REMARK       scalar :: %.2f \n", bestScaleCoef);
    fprintf(pFile, "# REMARK         QMIN :: %.6f \n", qvalues[0]);
    fprintf(pFile, "# REMARK         QMAX :: %.6f \n", qvalues[total-1]);
    fprintf(pFile, "# REMARK         LMAX :: %i \n", lmax);
    fprintf(pFile, "# REMARK         DMAX :: %.1f Angstroms\n", model.getDmax());
    fprintf(pFile, "# REMARK     RESIDUES :: %i \n", model.getPDBModel().getTotalResidues());
    fprintf(pFile, "# REMARK        ATOMS :: %i \n", model.getPDBModel().getTotalCoordinates());
    fprintf(pFile, "# REMARK (non-H) MASS :: %.1f \n", model.getPDBModel().getMW());
    fprintf(pFile, "# REMARK       VOLUME :: %.1f \n", model.getPDBModel().getVolume());
    fprintf(pFile, "# REMARK    COLUMNS 1 :: %s \n", "index");
    fprintf(pFile, "# REMARK    COLUMNS 2 :: %s \n", "momentum transfer (Angstroms^-1)");
    fprintf(pFile, "# REMARK    COLUMNS 3 :: %s \n", "I_observed");
    fprintf(pFile, "# REMARK    COLUMNS 4 :: %s \n", "I_model");
    fprintf(pFile, "# REMARK    COLUMNS 5 :: %s \n", "residual [I_obs - I_model]");

    for(unsigned int i=0; i < total; i++){
        float calc = i_best[i];
        float obs = iofqdata.getYDataAt(selectedIndices[i]);
        fprintf(pFile, "%5i %.3E %.5E %.5E %.5E\n", (i+1), qvalues[i], obs, calc, obs-calc);
    }

    fclose(pFile);
}
