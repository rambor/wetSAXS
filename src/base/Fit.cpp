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


Fit::Fit(unsigned int lmax, unsigned int division, float lowcx, float highcx, float colo, float cohi) :
        lmax(lmax),
        total_cx(division),
        lowcx(lowcx),
        highcx(highcx),
        low_c_occ(colo),
        high_c_occ(cohi){

    deltaCX = (highcx-lowcx)/(float)division;
    bestScore = FLT_MAX;
    bestCx = lowcx;

}


void Fit::assemblePartials(AtomisticModel &model,
                           Waters &waterModel,
                           std::vector<float> & norm_Pamps,
                           std::vector<float> & norm_Wamps,
                           std::vector<float> & norm_Camps,
                           std::vector<float> & aPW_cross_terms,
                           std::vector<float> & aPX_cross_terms,
                           std::vector<float> & aXW_cross_terms,
                           unsigned int totalInWorkingSet){

    const float * almPR = model.getModelPartialsReal();
    const float * almPI = model.getModelPartialsImag();
    const float * xlmR = model.getExcludedVolumePartialsReal();
    const float * xlmI = model.getExcludedVolumePartialsImag();

    const float * wlmR = waterModel.getModelPartialsReal();
    const float * wlmI = waterModel.getModelPartialsImag();

    float * pAp = norm_Pamps.data();
    float * pAw = norm_Wamps.data();
    float * pAc = norm_Camps.data();

    // cross-terms
    float * pXW_cross_term = aXW_cross_terms.data();
    float * pPW_cross_term = aPW_cross_terms.data();
    float * pPX_cross_term = aPX_cross_terms.data();

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

            a_p_norm +=  *pAlmPR * *pAlmPR + *pAlmPI * *pAlmPI;     // sum of the A*A_transpose
            a_w_norm +=  *pWlmR * *pWlmR + *pWlmI * *pWlmI;         // water layer

            // add hydration layer waters to the excluded volume
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
        pAw[q] = a_w_norm; //
        pAc[q] = a_x_norm; // => a^x -> scaled by contrast term

        //pApw[q] = a_p_norm + a_w_norm*occ*occ + 2.0f*occ*(aPWR_term + aPWI_term); // hydratedParticle
        pPW_cross_term[q] = 2.0f*(aPWR_term + aPWI_term); // hydratedParticle

        // invert phase - this is a cross term that remains that reflects the cross term with surrounding solvent (not layer)
        pPX_cross_term[q] = -2.0f*(aPXR_term + aPXI_term); // * c_x
        pXW_cross_term[q] = -2.0f*(aXWR_term + aXWI_term); // * c_x x c_o

    } // end loop over q values
}

/*
 * |A_p|^2 - |A_x|^2 + (A_p - A_x)*A_m' + (A_p - A_x')*A_m
 * |A_p|^2 - cx^2*|A_x|^2 + cx*(A_p*A_m' - A_p'*A_m') - cx^2*(A_x*A_m' + A_x'*A_m)
 * the last term is the Ax-Am cross-term
 * the second to last term is the Ap-Am cross-term
 * partials I need are:
 * |A_p|^2 : norm_aPs
 * |A_x|^2 : norm_aCs
 * A_p*A_m' - A_p'*A_m' => approximate using sin(q*r_pm)/(q*r_pm) weighted : aXW_cross_term
 * A_x*A_m' + A_x'*A_m  => approximate using sin(q*r_xm)/(q*r_xm) weighted : aPW_cross_term
 */
void Fit::assemblePartialsWithMask(AtomisticModel &model,
                           Mask &mask,
                           unsigned int totalInWorkingSet){

    const float * almPR = model.getModelPartialsReal();
    const float * almPI = model.getModelPartialsImag();
    const float * xlmIntensity = model.getExcludedVolumePartialsReal();

    const float * cross_Ap_Am = mask.getMaskCross_Ap_Am();
    const float * cross_Ax_Am = mask.getMaskCross_Ax_Am();

    float * pAp = norm_aPs.data();
    float * pAc = norm_aCs.data();

    // cross-terms
    float * pAXM_cross_term = aXW_cross_term.data(); // exVol mask cross-term
    float * pAPM_cross_term = aPX_cross_term.data(); // particle mask cross-term

    float norm_a_p, norm_a_x;

    const float * pAlmPR, * pAlmPI;

    unsigned int totalpartialsperq = (lmax+1)*(lmax+1);
    // calculate squared norm terms
    for(unsigned int q = 0; q < totalInWorkingSet; q++) {

        norm_a_p = 0.0;
        norm_a_x = 0.0;

        unsigned int startAt = q*totalpartialsperq; // always positive
        // for each l, m take the complex norm and sum
        for(unsigned int lm=0; lm< totalpartialsperq; lm++){ // sums the L and M's per q value

            unsigned int index = startAt + lm;

            pAlmPR = &almPR[index];
            pAlmPI = &almPI[index];

            // complex norm
            norm_a_p +=  *pAlmPR * *pAlmPR + *pAlmPI * *pAlmPI;     // particle scattering
            norm_a_x += xlmIntensity[index]; // excluded volume
        }

        // PI constant gets absorbed into the scaling constant below during the chi2 fit
        pAp[q] = norm_a_p; //*pi_constant;
        pAc[q] = norm_a_x; // => a^x -> scaled by contrast term c_x^2
        //pApw[q] = a_p_norm + a_w_norm*occ*occ + 2.0f*occ*(aPWR_term + aPWI_term); // hydratedParticle

        // investigate here
        pAPM_cross_term[q] = cross_Ap_Am[q]; // *c_x
        pAXM_cross_term[q] = cross_Ax_Am[q]; // *c_x^2

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
    float * pAw = norm_aWs.data();
    float * pAc = norm_aCs.data();

    // cross-terms
    float * pAXW_cross_term = aXW_cross_term.data();
    float * pAPX_cross_term = aPX_cross_term.data();
    float * pAPW_cross_term = aPW_cross_term.data();

    float a_p_norm, a_x_norm, a_w_norm, aPWR_term, aPWI_term;
    float x_w_r, x_w_i;
    float aPXR_term, aPXI_term, aXWR_term, aXWI_term;
    const float * pAlmPR, * pAlmPI, * pXlmR, * pXlmI, * pWlmR, *pWlmI;

    unsigned int totalpartialsperq = (lmax+1)*(lmax+1);

    // partialRows maps a working-set row to the row of the partials it was computed for
    // (chiFreeSearch computes partials once for all data); empty means the rows coincide
    if (!partialRows.empty() && partialRows.size() != totalInWorkingSet)
        throw std::invalid_argument("Fit::assemblePartials: partialRows does not match the working set");

    // calculate squared norm terms
    for(unsigned int q = 0; q < totalInWorkingSet; q++) {

        a_p_norm = 0.0;
        a_x_norm = 0.0;
        a_w_norm = 0.0;
        aPWR_term = 0.0;
        aPWI_term = 0.0;
        aPXR_term = 0.0;
        aPXI_term = 0.0;
        aXWR_term = 0.0;
        aXWI_term = 0.0;

        const unsigned int row = partialRows.empty() ? q : partialRows[q];
        unsigned int startAt = row*totalpartialsperq; // always positive
        // for each l, m take the complex norm and sum
        for(unsigned int lm=0; lm< totalpartialsperq; lm++){ // sums the L and M's per q value

            unsigned int index = startAt + lm;

            pAlmPR = &almPR[index];
            pAlmPI = &almPI[index];
            pXlmR = &xlmR[index];
            pXlmI = &xlmI[index];
            pWlmR = &wlmR[index];
            pWlmI = &wlmI[index];

            // complex norm
            a_p_norm +=  *pAlmPR * *pAlmPR + *pAlmPI * *pAlmPI;     // particle scattering
            a_w_norm +=  *pWlmR * *pWlmR + *pWlmI * *pWlmI;         // water layer

            x_w_r = *pXlmR + *pWlmR; // combining hydrating waters into the excluded volume
            x_w_i = *pXlmI + *pWlmI;

            a_x_norm +=  x_w_r*x_w_r + x_w_i*x_w_i; // excluded volume, includes the water that form the scattered object

            // particle - water cross term
            aPWR_term += *pAlmPR * *pWlmR + *pAlmPI * *pWlmI;
            //aPWI_term += *pAlmPI * *pWlmI;

            // excluded volume and water cross-term
            // excluded volume and particle cross-term
            // aPXR_term += x_w_r * (*pAlmPR + occ* *pWlmR); // invert the phase to make it cross term between particle and solvent
            // aPXI_term += x_w_i * (*pAlmPI + occ* *pWlmI); //

            aPXR_term += x_w_r * *pAlmPR; // invert the phase to make it cross term between particle and solvent
            aPXI_term += x_w_i * *pAlmPI; // multiply by cx

            aXWR_term += x_w_r * *pWlmR; // multiply by cx and c_occupancy in search
            aXWI_term += x_w_i * *pWlmI; // invert the phase to make it cross term between particle and solvent
        }

        // PI constant gets absorbed into the scaling constant below during the chi2 fit
        pAp[q] = a_p_norm; //*pi_constant;
        pAw[q] = a_w_norm; //
        pAc[q] = a_x_norm; // => a^x -> scaled by contrast term

        //pApw[q] = a_p_norm + a_w_norm*occ*occ + 2.0f*occ*(aPWR_term + aPWI_term); // hydratedParticle
        pAPW_cross_term[q] = 2.0f*aPWR_term; // hydratedParticle

        // invert phase - negative - this is a cross term that remains that reflects the cross term with surrounding solvent (not layer)
        //pAXW_cross_term[q] = -2.0f*(aXWR_term + aXWI_term + aPXR_term + aPXI_term);
        pAPX_cross_term[q] = -2.0f*(aPXR_term + aPXI_term); // * c_x
        pAXW_cross_term[q] = -2.0f*(aXWR_term + aXWI_term); // * c_x x c_o

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
void Fit::setUpForFittingWithMask(IofQData &iofqdata, AtomisticModel &model, Mask &maskModel) {

    unsigned int totalInWorkingSet = iofqdata.getTotalInWorkingSet();

    // set the size of the vectors
    norm_aPs.resize(totalInWorkingSet); // particle scattering
    norm_aCs.resize(totalInWorkingSet); // excluded volume

    aXW_cross_term.resize(totalInWorkingSet); // excludedVol-mask
    aPX_cross_term.resize(totalInWorkingSet); // particle-mask

    assemblePartialsWithMask(model, maskModel, totalInWorkingSet);

    // populate experimental observation information
    const std::vector<Datum> & workingSet = iofqdata.getWorkingSet();

    i_obs.resize(totalInWorkingSet);
    i_obs_over_variance.resize(totalInWorkingSet);
    inv_variance.resize(totalInWorkingSet);
    i_best.resize(totalInWorkingSet);

    for(unsigned int q =0; q < totalInWorkingSet; q++) {
        const Datum *pW = &workingSet[q];
        i_obs[q] = pW->getI();
        i_obs_over_variance[q] = pW->getI() * pW->getInvVar();
        inv_variance[q] = pW->getInvVar();
    }
}

// set up the containers for holding complex norms for calculating intensities
void Fit::setUpForFitting(IofQData &iofqdata, AtomisticModel &model, Waters &waterModel) {

    unsigned int totalInWorkingSet = iofqdata.getTotalInWorkingSet();

    norm_aPs.resize(totalInWorkingSet);
    norm_aWs.resize(totalInWorkingSet);
    norm_aCs.resize(totalInWorkingSet);
    aXW_cross_term.resize(totalInWorkingSet);
    aPX_cross_term.resize(totalInWorkingSet);
    aPW_cross_term.resize(totalInWorkingSet);


    assemblePartials(model, waterModel, totalInWorkingSet);

    // populate experimental observation information
    //std::vector<Datum> & workingSet = const_cast<std::vector<Datum> &>(iofqdata.getWorkingSet());
    const std::vector<Datum> & workingSet = iofqdata.getWorkingSet();

    i_obs.resize(totalInWorkingSet);
    i_obs_over_variance.resize(totalInWorkingSet);
    inv_variance.resize(totalInWorkingSet);
    i_best.resize(totalInWorkingSet);

    for(unsigned int q =0; q < totalInWorkingSet; q++) {
        const Datum *pW = &workingSet[q];
        i_obs[q] = pW->getI();
        i_obs_over_variance[q] = pW->getI() * pW->getInvVar();
        inv_variance[q] = pW->getInvVar();
    }
}


void Fit::searchWithMask(IofQData &iofqdata, AtomisticModel &model, Waters &waterModel, Mask &maskModel) {
    // |Ap|^2 + cx*(ApAm) - cx^2 (AxAm) - |Ax|^2

    // setup function with select q-values, model, watermodel and mask
    setUpForFittingWithMask(iofqdata, model, maskModel); // takes < 0.001 seconds

    float * pAp = norm_aPs.data();
    float * pAc = norm_aCs.data();

    // cross-terms
    float * pAXW_cross_term = aXW_cross_term.data(); // exVol-mask cross-term
    float * pAPX_cross_term = aPX_cross_term.data(); // particle-mask cross-term

    // perform grid search
    float sum, diff, scalar;
    bestCx = 0.0f;
    bestOcc = 1.0;
    bestScore = FLT_MAX;

    auto & qvalues = iofqdata.getWorkingSetQvalues();
    auto pQ = qvalues.data();

    unsigned int totalInWorkingSet = qvalues.size();
    std::vector<float> residuals(totalInWorkingSet);
    std::vector<float> i_calc(totalInWorkingSet);
    auto pIcalc = i_calc.data();
    float chi, dw;

    // setting match point requires total waters in hydration model
    model.setMatchPointScale(waterModel.getTotalWaters());

    float contrastMatchScale = model.getMatchPoint(); // scale factor required to match forward scattering as if water
    SASTOOLS_UTILS_H::logger("Match Point", formatNumber(model.getMatchPoint(),3));

    float inv_n_s = 1.0f/(float)(totalInWorkingSet-2);
    float cx;
    float c_occ = 0.1;//low_c_occ;

    cx =  lowcx*contrastMatchScale;

    for(int incr_cx = 0; incr_cx < total_cx; incr_cx++){

        scalar = calculateIofQfromModelWithVariancesMask(totalInWorkingSet,
                                     pQ,
                                     pAp,
                                     pAc,
                                     pAXW_cross_term,
                                     pAPX_cross_term,
                                     cx,
                                     c_occ,
                                     pIcalc);

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
            bestCx = cx;
            bestScaleCoef = scalar;
            bestDW = dw;
            bestOcc = c_occ;

            // do point transfer and do the scalar after finishing
            for(unsigned int q =0; q < totalInWorkingSet; q++) {
                i_best[q] = scalar*pIcalc[q];
            }
        }
std::cout << "cx " << cx << std::endl;
        cx += deltaCX*contrastMatchScale;
    }

    std::cout << "_______________________________________________________________________" << std::endl;

    printBestScores();

    header = iofqdata.getHeader();
    std::vector<float> qrvalues = iofqdata.getWorkingSetQvalues();
    writeBestModelToFile(qrvalues, iofqdata.getFilename(), model);

}


void Fit::search(IofQData &iofqdata, AtomisticModel &model, Waters &waterModel) {

    setUpForFitting(iofqdata, model, waterModel); // takes < 0.001 seconds

    float * pAc = norm_aCs.data();
    float * pAp = norm_aPs.data();
    float * pAw = norm_aWs.data();

    // cross-terms
    float * pAXW_cross_term = aXW_cross_term.data();
    float * pAPX_cross_term = aPX_cross_term.data();
    float * pAPW_cross_term = aPW_cross_term.data();


    // perform grid search
    float sum, diff, scalar;
    bestCx = 0.0f;
    bestOcc = 1.0;
    bestScore = FLT_MAX;

    auto & qvalues = iofqdata.getWorkingSetQvalues();
    auto pQ = qvalues.data();

    unsigned int totalInWorkingSet = qvalues.size();
    std::vector<float> residuals(totalInWorkingSet);
    std::vector<float> i_calc(totalInWorkingSet);
    auto pIcalc = i_calc.data();

    float chi, dw;

    // setting match point requires total waters in hydration model
    model.setMatchPointScale(waterModel.getTotalWaters());
    float contrastMatchScale = model.getMatchPoint(); // scale factor required to match forward scattering as if water
    SASTOOLS_UTILS_H::logger("Match Point", formatNumber(model.getMatchPoint(),3));

    float inv_n_s = 1.0f/(float)(totalInWorkingSet-2);
    float cx;
    float c_occ = 0.0;//low_c_occ;
    float incr_occ = 0.12;

    //while (c_occ < high_c_occ) {
    while (c_occ < 10.01) {

        cx =  lowcx*contrastMatchScale;

        for(int incr_cx = 0; incr_cx < total_cx; incr_cx++){

                scalar = calculateIofQfromModelWithVariances(totalInWorkingSet,
                                             pQ,
                                             pAp,
                                             pAc,
                                             pAw,
                                             pAXW_cross_term,
                                             pAPX_cross_term,
                                             pAPW_cross_term,
                                             cx,
                                             c_occ,
                                             pIcalc);

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
                    bestCx = cx;
                    bestScaleCoef = scalar;
                    bestDW = dw;
                    bestOcc = c_occ;

                    // do point transfer and do the scalar after finishing
                    for(unsigned int q =0; q < totalInWorkingSet; q++) {
                        i_best[q] = scalar*pIcalc[q];
                    }
                }

            cx += deltaCX*contrastMatchScale;
        }
        c_occ += incr_occ;
    }


    SASTOOLS_UTILS_H::logger("CX RANGE", formatNumber(lowcx*contrastMatchScale,3) + " to " + formatNumber(cx,3));

    std::cout << "_______________________________________________________________________" << std::endl;

    printBestScores();

    // if (bestOcc > 1.05 && bestCx > 1.05) {
    //     for(unsigned int q =0; q < totalInWorkingSet; q++) {
    //         std::cout << q << " <-> " << qvalues[q] << " - - " <<  i_obs[q] << std::endl;
    //     }
    // }

    SASTOOLS_UTILS_H::logger("Effective Solvent Density", formatNumber(bestCx,3));
    SASTOOLS_UTILS_H::logger("Percent Adjustment > H2O", formatNumber(bestCx/contrastMatchScale,3));
//    SASTOOLS_UTILS_H::logger("Use this as CX", formatNumber(bestCx/contrastMatchScale,3));

    header = iofqdata.getHeader();

    std::vector<float> qrvalues = iofqdata.getWorkingSetQvalues();
    writeBestModelToFile(qrvalues, iofqdata.getFilename(), model);
    //writeBestModelToFile(const_cast<std::vector<float> &>(qvalues), model, pApw, pAc, pAXW_cross_term);

}

void Fit::printBestScores(){
    SASTOOLS_UTILS_H::logger("Best Score", formatNumber(bestScore,3));
    SASTOOLS_UTILS_H::logger("chi-squared", formatNumber(bestChi,3));
    SASTOOLS_UTILS_H::logger("Durbin-Watson", formatNumber(bestDW,3));
//    SASTOOLS_UTILS_H::logger("B-factor", formatNumber(bestBfactor,1));
    SASTOOLS_UTILS_H::logger("contrast(c_x)", formatNumber(bestCx,2));
    SASTOOLS_UTILS_H::logger("occupancy(c_o)", formatNumber(bestOcc,2));
}


void Fit::searchCOFactor(IofQData &iofqdata, AtomisticModel &model, Waters &waterModel, float fixedCx, float colo, float cohi, int totalCO) {

    setUpForFitting(iofqdata, model, waterModel);

    const float delta_co = (cohi - colo)/(float)totalCO; //
    unsigned int totalInWorkingSet = iofqdata.getTotalInWorkingSet();

    norm_aPs.resize(totalInWorkingSet);
    norm_aWs.resize(totalInWorkingSet);
    norm_aCs.resize(totalInWorkingSet);
    aXW_cross_term.resize(totalInWorkingSet);
    aPX_cross_term.resize(totalInWorkingSet);
    aPW_cross_term.resize(totalInWorkingSet);


    float * pAp = norm_aPs.data();
    float * pAw = norm_aWs.data();
    float * pAc = norm_aCs.data();

    // cross-terms
    float * pAxw_cross_term = aXW_cross_term.data();
    float * pApx_cross_term = aPX_cross_term.data();
    float * pApw_cross_term = aPW_cross_term.data();

    // perform grid search
    float sum, diff, scalar;
    float chi, dw;

    bestCx = 0.0f;
    bestOcc = colo;
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
    float co_factor = colo;

    float inv_n_s = 1.0f/(float)(totalInWorkingSet-2);

    for (int incr_co = 0; incr_co < totalCO; incr_co++){

        scalar = calculateIofQfromModelWithVariances(totalInWorkingSet,
            pQ,
            pAp,
            pAc,
            pAw,
            pAxw_cross_term,
            pApx_cross_term,
            pApw_cross_term,
            cx,
            co_factor,
            pIcalc);

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
            bestOcc = co_factor;
            bestScaleCoef = scalar;
            bestDW = dw;

            for(unsigned int q =0; q < totalInWorkingSet; q++) {
                i_best[q] = scalar * pIcalc[q];
            }
        }

        co_factor += delta_co;
    }
}


void Fit::searchCX(IofQData &iofqdata, AtomisticModel &model, Waters &waterModel, float fixedCO, float lowerCX, float upperCX, int divisions) {

    lowcx = lowerCX;
    highcx = upperCX;

    setUpForFitting(iofqdata, model, waterModel);

    deltaCX = (highcx-lowcx)/(float)divisions;

    unsigned int totalInWorkingSet = iofqdata.getTotalInWorkingSet();

    float * pAp = norm_aPs.data();
    float * pAc = norm_aCs.data();
    float * pAw = norm_aWs.data();

    // cross-terms
    float * pAXW_cross_term = aXW_cross_term.data();
    float * pAPW_cross_term = aPW_cross_term.data();
    float * pAPX_cross_term = aPX_cross_term.data();

    // perform grid search
    float sum, diff, scalar;
    float chi, dw;

    bestCx = 0.0f;
    bestOcc = fixedCO;
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
                                                     pAp,
                                                     pAc,
                                                     pAw,
                                                     pAXW_cross_term,
                                                     pAPX_cross_term,
                                                     pAPW_cross_term,
                                                     cx,
                                                     fixedCO,
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
void Fit::fitFixedCoAndCx(IofQData &iofqdata,
                               AtomisticModel &model,
                               Waters &waterModel,
                               float fixedCO,
                               float fixedCX) {

    setUpForFitting(iofqdata, model, waterModel);

    float * pAc = norm_aCs.data();
    float * pAp = norm_aPs.data();
    float * pAw = norm_aWs.data();

    // cross-terms
    float * pAxw_cross_term = aXW_cross_term.data();
    float * pApx_cross_term = aPX_cross_term.data();
    float * pApw_cross_term = aPW_cross_term.data();

    auto & qvalues = iofqdata.getWorkingSetQvalues();
    auto pQ = qvalues.data();

    unsigned int totalInWorkingSet = iofqdata.getTotalInWorkingSet();
    std::vector<float> residuals(totalInWorkingSet);
    std::vector<float> i_calc(totalInWorkingSet);
    auto pIcalc = i_calc.data();

    float inv_n_s = 1.0f/(float)(totalInWorkingSet-2);

    float scalar = calculateIofQfromModelWithVariances(totalInWorkingSet,
                             pQ,
                             pAp,
                             pAc,
                             pAw,
                             pAxw_cross_term,
                             pApx_cross_term,
                             pApw_cross_term,
                             fixedCX,
                             fixedCO,
                             pIcalc);

    // score the model
    float sum = 0.0f, diff;
    scalar = (bestScaleCoef < 0) ? scalar : bestScaleCoef;

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
    bestOcc = fixedCO;
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
                           float co,
                           std::vector<float> & qvalues,
                           AtomisticModel & model,
                           std::vector<float> & icalc,
                           const float * pAp,  // particle amplitude
                           const float * pAw,  // hydrated particle amplitude
                           const float * pAc,  // excluded volume
                           const float * pPW_cross_term, // cross-term
                           const float * pPX_cross_term, // cross-term
                           const float * pXW_cross_term // cross-term
                           ){

    std::string name = model.getPDBModel().getFileStemName();

    char buffer [50];

    sprintf (buffer, "CO_%.2f", co); // change the number to a string for name
    std::string sco = buffer;
    std::replace(sco.begin(), sco.end(), '.', 'p');

    sprintf (buffer, "CX_%.2f", cx); // change the number to a string for name
    std::string cxs = buffer;
    std::replace(cxs.begin(), cxs.end(), '.', 'p');

    std::string nameOf = name + "_" + cxs + "_" + sco + "_calc.dat";

    FILE * pFile = fopen(nameOf.c_str(), "w");
    unsigned int total = qvalues.size();

    fprintf(pFile, "# REMARK           :: %s\n", model.getPDBModel().getFilename().c_str());
    fprintf(pFile, "# REMARK        CX :: %.4f \n", cx);
    fprintf(pFile, "# REMARK        CO :: %.2f \n", co);
    fprintf(pFile, "# REMARK      QMIN :: %.6f \n", qvalues[0]);
    fprintf(pFile, "# REMARK      QMAX :: %.6f \n", qvalues[total-1]);
    fprintf(pFile, "# REMARK      LMAX :: %i \n", lmax);
    fprintf(pFile, "# REMARK      DMAX :: %.1f \n", model.getDmax());
    fprintf(pFile, "# REMARK  RESIDUES :: %i \n", model.getPDBModel().getTotalResidues());
    fprintf(pFile, "# REMARK     ATOMS :: %i \n", model.getPDBModel().getTotalCoordinates());
    fprintf(pFile, "# REMARK COLUMNS 1 :: %s \n", "index");
    fprintf(pFile, "# REMARK COLUMNS 2 :: %s \n", "momentum transfer (inverse Angstroms");
    fprintf(pFile, "# REMARK COLUMNS 3 :: %s \n", "Icalc");
    fprintf(pFile, "# REMARK COLUMNS 4 :: %s \n", "in vacuo particle");
    fprintf(pFile, "# REMARK COLUMNS 5 :: %s \n", "hydrated particle");
    fprintf(pFile, "# REMARK COLUMNS 6 :: %s \n", "excluded volume");
    fprintf(pFile, "# REMARK COLUMNS 7 :: %s \n", "cross-term Particle and hydrating waters");
    fprintf(pFile, "# REMARK COLUMNS 8 :: %s \n", "cross-term Particle and background solvent");
    fprintf(pFile, "# REMARK COLUMNS 9 :: %s \n", "cross-term Background solvent and hydrating waters");

    float co2 = co*co;
    float cx2 = cx*cx;
    float cocx = co*cx;

    for(unsigned int i=0; i < total; i++){

        fprintf(pFile, "%5i %.3E %.5E %.5E %.5E %.5E %.5E %.5E %.5E\n", (i+1), qvalues[i], icalc[i], pAp[i], co2*pAw[i],
            cx2*pAc[i], co*pPW_cross_term[i], cx*pPX_cross_term[i], cocx*pXW_cross_term[i]);
    }

    fclose(pFile);
}


void Fit::writeBestModelToFile(
                           std::vector<float> & qvalues,
                           std::string filename,
                           AtomisticModel & model
//                           const float * aP,  // particle amplitude
//                           const float * aC,  // excluded volume
//                           const float * aCXP // cross-term
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
    fprintf(pFile, "# REMARK           CO :: %.2f \n", bestOcc);
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

    const float * aP = norm_aPs.data();
    const float * aC = norm_aCs.data();
    //const float * aW = norm_aWs.data();

    for(unsigned int i=0; i < total; i++){
        fprintf(pFile, "%5i %.3E %.5E %.5E %.5E %.5E %.5E\n", (i+1), qvalues[i], i_obs[i], i_best[i], (i_obs[i] - i_best[i]), aP[i], aC[i]);
    }

    fclose(pFile);
}

/*
 * Calculate SAXS profile
 */
float Fit::calculateIofQfromModelWithVariances(unsigned int totalInWorkingSet,
                                             const float * pQ,
                                             const float * pAp,
                                             const float * pAc,
                                             const float * pAw,
                                             const float * pAXW_cross_term,
                                             const float * pAPX_cross_term,
                                             const float * pAPW_cross_term,
                                             float cx,
                                             float co,
                                             float * pIcalc) {
    float i_of_q;
    const float * q_val;

    float c_obs_calc = 0.0f;
    float calc_var = 0.0f;

    float cob2, cob;
    float co2 = co*co;
    float cx2 = cx*cx;
    float cocx = co*cx;
    float bf = 6;
    float inv = 1.0/(M_PI*M_PI*8);

    for(unsigned int q =0; q < totalInWorkingSet; q++) {
        q_val = &pQ[q]; // get qvalue to calculate b-factor
        cob2 = expf(-(*q_val * *q_val)*bf*inv);
        cob = expf(-(*q_val * *q_val)*bf*inv16Pi2);
        //pApw[q] = a_p_norm + a_w_norm*occ*occ + 2.0f*occ*(aPWR_term + aPWI_term); // hydratedParticle
         i_of_q = pAp[q] + co*pAPW_cross_term[q] + pAw[q]*co2
         + cocx * pAXW_cross_term[q]*cob + cx * pAPX_cross_term[q]*cob - 2*pAc[q]*cx2*cob2;


        // aPXR_term += x_w_r * (*pAlmPR + occ* *pWlmR); // split into two -> Wlm is modified occupancy (co)
        // aPXI_term += x_w_i * (*pAlmPI + occ* *pWlmI); // x_w is excluded volume and modified by cx (contrast)
        // pAPX_cross_term => x_w_r * *pAlmPR + x_w_i * *pAlmPI
        // pAXW_cross_term => x_w_r * *pAlmWR + x_w_i * *pWlmPI

        c_obs_calc += i_of_q*i_obs_over_variance[q];
        calc_var += i_of_q*i_of_q*inv_variance[q];

        pIcalc[q] = i_of_q;
    }

    return c_obs_calc/calc_var;
}

/*
 * Calculate SAXS profile
 */
float Fit::calculateIofQfromModelWithVariancesMask(unsigned int totalInWorkingSet,
                                             const float * pQ,
                                             const float * pAp,
                                             const float * pAc,
                                             const float * pAXW_cross_term,
                                             const float * pAPX_cross_term,
                                             float cx,
                                             float co,
                                             float * pIcalc) {
    float i_of_q;
    const float * q_val;

    float c_obs_calc = 0.0f;
    float calc_var = 0.0f;

    float cob2, cob;
    float co2 = co*co;
    float cx2 = cx*cx;
    float cocx = co*cx;
    float bf = 81;
    float inv = 1.0/(M_PI*M_PI*8);
    float factor = 1;//16*M_PI*M_PI;

    for(unsigned int q =0; q < totalInWorkingSet; q++) {
        q_val = &pQ[q]; // get qvalue to calculate b-factor
        cob2 = expf(-(*q_val * *q_val)*bf*inv);
        cob = expf(-(*q_val * *q_val)*bf*inv16Pi2);
        //pApw[q] = a_p_norm + a_w_norm*occ*occ + 2.0f*occ*(aPWR_term + aPWI_term); // hydratedParticle
        i_of_q = factor*pAp[q] + pAPX_cross_term[q]*cx - pAXW_cross_term[q]*cx2 +  -  pAc[q]*cx2*factor;
//std::cout << q << " " << (pAp[q] - 7.910 * pAc[q]) << std::endl;
        c_obs_calc += i_of_q*i_obs_over_variance[q];
        calc_var += i_of_q*i_of_q*inv_variance[q];

        pIcalc[q] = i_of_q;
    }

    return c_obs_calc/calc_var;
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
    float i_of_q;
    const float * q_val;
    float exp_term;

    float c_obs_calc = 0.0f;
    float calc_var = 0.0f;

    for(unsigned int q =0; q < totalInWorkingSet; q++) {

        q_val = &pQ[q]; // get qvalue to calculate b-factor
        exp_term = expf(-(*q_val * *q_val)*b_factor*inv16Pi2)*cx;

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
                                               const float * pAp,
                                               const float * pAw,
                                               const float * pAc,
                                               const float * pAPW_cross_term,
                                               const float * pAPX_cross_term,
                                               const float * pAXW_cross_term,
                                               float cx,
                                               float co,
                                               float * pIcalc
                                               ){

    float co2 = co*co;
    float cx2 = cx*cx;
    float cocx = co*cx;

    for(unsigned int q =0; q < totalInWorkingSet; q++) {

        // q_val = &pQ[q]; // get qvalue to calculate b-factor
        // cob = expf(-(*q_val * *q_val)*co*inv16Pi2);

        //pApw[q] = a_p_norm + a_w_norm*occ*occ + 2.0f*occ*(aPWR_term + aPWI_term); // hydratedParticle

        pIcalc[q] = pAp[q] + pAw[q]*co2 + co*pAPW_cross_term[q]
            + cocx * pAXW_cross_term[q] + cx * pAPX_cross_term[q] + pAc[q]*cx2;
    }


}

void Fit::simulate(float cx,
                   float co,
                   std::vector<float> & qvalues,
                   AtomisticModel &model,
                   Waters &waterModel){


    unsigned int totalInWorkingSet = qvalues.size();

    std::vector<float> norm_aPs(totalInWorkingSet);
    std::vector<float> norm_aWs(totalInWorkingSet);
    std::vector<float> norm_aCs(totalInWorkingSet);
    std::vector<float> aPW_cross_term(totalInWorkingSet);
    std::vector<float> aPX_cross_term(totalInWorkingSet);
    std::vector<float> aXW_cross_term(totalInWorkingSet);


    assemblePartials(model,
                     waterModel,
                     norm_aPs,
                     norm_aWs,
                     norm_aCs,
                     aPW_cross_term,
                     aPX_cross_term,
                     aXW_cross_term,
                     totalInWorkingSet);

    float * pAp = norm_aPs.data();
    float * pAc = norm_aCs.data();
    float * pAw = norm_aWs.data();

    float * pPW_cross_term = aPW_cross_term.data();
    float * pPX_cross_term = aPX_cross_term.data();
    float * pXW_cross_term = aXW_cross_term.data();

    std::vector<float> i_calc(totalInWorkingSet);
    auto pIcalc = i_calc.data();


    calculateIofQfromModel(totalInWorkingSet,
                           qvalues.data(),
                           pAp,
                           pAw,
                           pAc,
                           pPW_cross_term,
                           pPX_cross_term,
                           pXW_cross_term,
                           cx,
                           co,
                           pIcalc
                           );


    writeModelToFile(cx, co, qvalues, model, i_calc, pAp, pAw, pAc, pPW_cross_term, pPX_cross_term, pXW_cross_term);
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

    // Partial amplitudes of each q are independent of the others, so compute them ONCE for every
    // measured q. Each round then selects its working-set rows through partialRows.
    iofqdata.setAllDataToWorkingSet();
    std::vector<float> allQvalues = iofqdata.getWorkingSetQvalues();

    SASTOOLS_UTILS_H::logger("", "Calculating Amplitudes (all q, once)");
    model.calculatePartialAmplitudes(lmax, allQvalues.size(), allQvalues);
    waterModel.calculatePartialAmplitudes(lmax, allQvalues.size(), allQvalues);

    for (unsigned int rnd=0; rnd<totalRounds; rnd++){

        SASTOOLS_UTILS_H::logger("CHI-free round", formatNumber(rnd+1));

        iofqdata.makeWorkingSet();                   // sorted indices into the full data
        partialRows = iofqdata.getSelectedIndices(); // working-set row -> partials row

        this->search(iofqdata, model, waterModel);
        results.emplace_back(Result(bestScore, bestScaleCoef, bestCx, bestOcc, bestDW, bestChi, i_best, iofqdata.getSelectedIndices()));
    } // should converge to some stable number - after how many samplings?

    // now sort and take the median of Results
    // Using lambda expressions in C++11
    logger("Sorting", "::");
    std::sort(results.begin(), results.end(), [](const Result& lhs, const Result& rhs) {
        return lhs.score < rhs.score;
    });

    std::cout << "_______________________________________________________________________" << std::endl;
    std::cout << " SORTED SCORES " << std::endl;
    std::cout << "   SCORE ::   CHI   ::   DW  ::   CX  ::   CO " << std::endl;

    char buffer [50];
    for(auto & res : results){
        sprintf (buffer, "%8.3f    %6.3f     %5.3f    %5.3f    %5.2f",
            res.score,
            res.getChi(),
            res.getDurbinWatson(),
            res.getCX(),
            res.getCO()); // change the number to a string for name

        std::cout << buffer << std::endl;
    }

    std::cout << "_______________________________________________________________________" << std::endl;
    // median
    int middle = (totalRounds-1)/2;

    auto & res = results[middle];
    bestScore = res.score;
    bestChi = res.getChi();
    bestDW = res.getDurbinWatson();
    bestCx = res.getCX(); // percent * contrast_match_point
    bestOcc = res.getCO();
    i_best = res.getICalc();
    bestScaleCoef = res.getScale();

    SASTOOLS_UTILS_H::logger("STATUS", "Writing median best-fit model");
    writeBestChiFreeModel(model, res.getSelectedIndices(), iofqdata, totalRounds);
    SASTOOLS_UTILS_H::logger("CHI-free Parameters", "Median From Search");
    printBestScores();

    // calculate over all input data using best fitting median parameters
    SASTOOLS_UTILS_H::logger("STATUS", "FINAL FIT TO ALL DATA");
    iofqdata.setAllDataToWorkingSet();
    partialRows.clear(); // working set = all data, in the same order as the partials
    std::vector<float> qvalues = iofqdata.getWorkingSetQvalues();

    // partials for all q already exist: no recalculation needed
    fitFixedCoAndCx(iofqdata, model, waterModel, bestOcc, bestCx);
    SASTOOLS_UTILS_H::logger("CHI-free Parameters", "All Data");
    printBestScores();

    writeBestModelToFile(qvalues, iofqdata.getFilename(), model);

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

    sprintf (buffer, "CO_%.2f", bestOcc); // change the number to a string for name
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
    fprintf(pFile, "# REMARK           CO :: %.2f \n", bestOcc);
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
