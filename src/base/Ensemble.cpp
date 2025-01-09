// Copyright (c) 2023.
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
// Created by Robert Rambo on 21/09/2023.
//

#include <random>
#include "Ensemble.h"
#include "AtomisticModel.h"
#include "Waters.h"


Ensemble::Ensemble(std::vector<std::string> file_names, std::string directory) {

    // check files exits and directory is writable
    for (auto & filename : file_names){

        try{

            FileClass basefile = FileClass(filename);

            if (basefile.isPDB()){
                filenames.push_back(filename);
            }

        } catch (const std::invalid_argument& fl){
            std::cerr << " " << fl.what() << std::endl;
        }
    }


    // check if directory exists and is writeable
    if (!boost::filesystem::exists(directory)){
        throw std::invalid_argument("** ERROR FILE => DIRECTORY NOT FOUND : " + directory);
    } else {
        this->directory = directory;
    }
}

/*
 * Calculate the partials for each model using workingset.
 * Then run a similarity algorithm for all pairwise comparisons to make a clustering
 * sample from each cluster
 */
void Ensemble::search(unsigned int lmax, IofQData * iofqfile, bool forceRNA, bool keepWaters, int oversampling){

    // for each file in filenames, calculate intensities using workingset
    iofqfile->makeWorkingSet(oversampling);

    unsigned int totalInWorkingSet = iofqfile->getTotalInWorkingSet();
    auto qvalues = iofqfile->getWorkingSetQvalues();

    std::vector<Amplitudes> models = std::vector<Amplitudes>(); // Amplitudes is a struct defined in Ensemble.h, holds partials

    unsigned int total = 0;
    std::vector<unsigned int> indices;

    // read files and calculate amplitudes along with hydration
    for (auto & filename : filenames){

            try{
                AtomisticModel model = AtomisticModel(filename, forceRNA, keepWaters);

                model.calculatePartialAmplitudes(lmax, totalInWorkingSet, qvalues, false);
                if (model.getPDBModel().getTotalResidues() < 1){
                    throw (filename);
                }

                // create water model
                Waters waters = Waters();
                waters.hydrateAtomisticModel(model); // add waters to PDB
                waters.createSphericalCoordinateOfHydration();
                waters.calculatePartialAmplitudes(lmax, totalInWorkingSet, qvalues);

                models.emplace_back(Amplitudes(total, model.getPDBModel().getFileStemName()));
                Amplitudes * pCalc = &models[total]; //

                std::vector<float> aPs = std::vector<float>(totalInWorkingSet);
                std::vector<float> aPWs = std::vector<float>(totalInWorkingSet);
                std::vector<float> aCs = std::vector<float>(totalInWorkingSet);
                std::vector<float> aXW_cross_term = std::vector<float>(totalInWorkingSet);

                // assemble the partials and
                assemblePartials(lmax, model, waters, totalInWorkingSet, aPs, aPWs, aCs, aXW_cross_term);
                pCalc->addValues(aPs, aPWs, aCs, aXW_cross_term);

                indices.push_back(total);
                total++;
            } catch (std::string & fl){
                SASTOOLS_UTILS_H::logger("Improper PDB file", "skipping..." + fl);
            }
    }

    std::vector<float> basis(totalInWorkingSet, 0.0f); // holds the sum of the model particle amplitudes
    for(auto & amps : models){ // for each model, add the particle amplitudes
        auto & rAmps = amps.aPs;
        for(unsigned int q=0; q<totalInWorkingSet; q++){
            basis[q] += rAmps[q];
        }
    }

    // average the values within basis to act as normalization
    for(auto & val : basis){
        val /= (float)total;
    }

    // perform similarity calculation using Durbin Watson
    std::vector<float> residuals(totalInWorkingSet);
    std::vector<Score> scores;
    for (unsigned int i=0; i<(total-1); i++){
        Amplitudes * mod1 = &models[i];
        //auto & amp1 = mod1->pAPWs;
        auto & amp1 = mod1->aPs;
        unsigned int k=(i+1);

        for (; k<total; k++){

            Amplitudes * mod2 = &models[k];
            //auto & amp2 = mod2->pAPWs; // might need to scale?080989090999999
            auto & amp2 = mod2->aPs; // might need to scale?080989090999999

            for(unsigned int q=0; q<totalInWorkingSet; q++){
                residuals[q] = (amp1[q] - amp2[q])/basis[q]; // normalized residuals
            }
            scores.emplace_back(Score(i, k, calculateDurbinWatson(residuals.data(), totalInWorkingSet)));
        }
    }


    float max = 0;
    float min = 16.1;
    float sum=0;
    unsigned int max_row, max_col;
    for(auto & score : scores){ // find pair with maximum DW score
        sum += score.score;
        if (score.score > max ){
            max = score.score;
            max_row = score.row;
            max_col = score.col;
        }

        if(score.score < min){
            min = score.score;
        }

       // std::cout << score.row << " <-> " << score.col <<  " :: " << score.score << std::endl;
    }

    float average = sum/(float)scores.size();
    std::cout << " Average " << average << std::endl;
    std::cout << " MAX " << max << std::endl;
    std::cout << " MIN " << min << std::endl;

    /*
     * for two clusters, assign the remaining models to either max row or max col depending on their respective DW scores
     */
    unsigned int index;
    Score * pS1, * pS2;

    std::map<unsigned int, std::vector<unsigned int>> groups;
    groups.insert(std::pair<int, std::vector<unsigned int>>(max_row,std::vector<unsigned int>()));
    groups.insert(std::pair<int, std::vector<unsigned int>>(max_col,std::vector<unsigned int>()));

    for (unsigned int i=0; i<total; i++){

        if (i != max_row && i != max_col){

            if (i<max_row){
                index = i*total - (i*(i+1)/2) + max_row - (i + 1);
                pS1 = &scores[index];
            } else {
                index = max_row*total - (max_row*(max_row+1)/2) + i - (max_row + 1);
                pS1 = &scores[index];
            }

            if (i < max_col){
                index = i*total - (i*(i+1)/2) + max_col - (i + 1);
                pS2 = &scores[index];
            } else {
                index = max_col*total - (max_col*(max_col+1)/2) + i - (max_col + 1);
                pS2 = &scores[index];
            }

            // add model to the group with lowest score
            if (pS1->score < pS2->score){
                groups[max_row].push_back(i);
            } else {
                groups[max_col].push_back(i);
            }
        }
    }

    std::cout << groups[max_row].size() << " :: " << groups[max_col].size() << std::endl;

    for(auto & grp : groups){
        std::cout << " Parent " << grp.first << std::endl;
        for(auto & ind : grp.second){
            std::cout << "     :: " << ind << std::endl;
        }
    }


    // for a given model, find closest neighbor?
    // demonstrate development is beyond current state of the art

    unsigned int ensemble_limit = 5; // should be Lmax

    // pair with the largest distance is trivial
    for(unsigned int ensemble_size=2; ensemble_size<ensemble_limit; ensemble_size++) {

        float maxd = 0.0f;
        std::vector<unsigned int> result(ensemble_size);
        std::vector<unsigned int> best(ensemble_size);
        combinations(indices, scores, ensemble_size, 0, result, best, maxd);

//        std::cout << " Best " << ensemble_size << " " << maxd << " ";
//        for(auto & bb : best){
//            std::cout << bb << " ";
//        }
//        std::cout << std::endl;
    }

    // assign each remaining conformation to member of best.




}

void Ensemble::linearCombination(unsigned int lmax, IofQData * iofqfile,
                                 unsigned int total_cx,
                                 float lowcx,
                                 float highcx,
                                 float upper_Bfactor,
                                 bool forceRNA, bool keepWaters, int oversampling){

    unsigned int information_limit = lmax - 2;
    if (filenames.size() > information_limit){
        std::cerr << "Too many models for independent linear combination " << std::endl;
    }

    iofqfile->makeWorkingSet(oversampling > filenames.size() ? oversampling : filenames.size());

    unsigned int totalInWorkingSet = iofqfile->getTotalInWorkingSet();
    auto qvalues = iofqfile->getWorkingSetQvalues();
    auto pQ = qvalues.data();
    const std::vector<Datum> & workingSet = iofqfile->getWorkingSet();

    std::vector<float> i_obs_over_variance(totalInWorkingSet);
    std::vector<float> inv_variance(totalInWorkingSet);

    for(unsigned int q =0; q < totalInWorkingSet; q++) {
        const Datum *pW = &workingSet[q];
        i_obs_over_variance[q] = pW->getI() * pW->getInvVar();
        inv_variance[q] = pW->getInvVar();
    }

    std::vector<Amplitudes> models = std::vector<Amplitudes>();

    float contrastMatchScale=1.0;
    unsigned int total_models = 0;
    std::vector<unsigned int> indices;

    bool setMatch = true;
    for (auto & filename : filenames){

        try{
            AtomisticModel model = AtomisticModel(filename, forceRNA, keepWaters);

            model.calculatePartialAmplitudes(lmax, totalInWorkingSet, qvalues, false);

            if (model.getPDBModel().getTotalResidues() < 1){
                throw (filename);
            }

            // create water model
            Waters waters = Waters();
            waters.hydrateAtomisticModel(model); // add waters to PDB
            waters.createSphericalCoordinateOfHydration();
            waters.calculatePartialAmplitudes(lmax, totalInWorkingSet, qvalues);

            if (setMatch){
                model.setMatchPointScale(waters.getTotalWaters());
                contrastMatchScale = model.getMatchPoint();
                setMatch = false;
            }
            models.emplace_back(Amplitudes(total_models, model.getPDBModel().getFileStemName()));
            Amplitudes * pCalc = &models[total_models];

            std::vector<float> aPs = std::vector<float>(totalInWorkingSet);
            std::vector<float> aPWs = std::vector<float>(totalInWorkingSet);
            std::vector<float> aCs = std::vector<float>(totalInWorkingSet);
            std::vector<float> aXW_cross_term = std::vector<float>(totalInWorkingSet);

            // assemble the partials and
            assemblePartials(lmax, model, waters, totalInWorkingSet, aPs, aPWs, aCs, aXW_cross_term);
            pCalc->addValues(aPs, aPWs, aCs, aXW_cross_term);

            indices.push_back(total_models);
            total_models++;
        } catch (std::string & fl){
            SASTOOLS_UTILS_H::logger("Improper PDB file", "skipping..." + fl);
        }
    }

    boost::numeric::ublas::matrix<float> model_intensities(totalInWorkingSet, total_models);

    // perform grid search
    float deltaCX;
    const float delta_b = 0.331f; // actual equation is exp(-q^2*d_w/3), d_w goes to 100, divided by 3 is 33
    float chi, b_factor, dw;

    auto total_b = (int)std::ceil(upper_Bfactor/delta_b); // divided by 3

    float cx =  lowcx*contrastMatchScale, scalar;
    float inv_n_s = 1.0f/(float)(totalInWorkingSet-2);


    for(int incr_cx = 0; incr_cx < total_cx; incr_cx++){

        b_factor = 0.0f;

        for (int incr_b = 0; incr_b < total_b; incr_b++){

            scalar = populateIntensitiesMatrix(totalInWorkingSet,
                                                         pQ,
                                                         models,
                                                         cx,
                                                         b_factor,
                                                         model_intensities
            );

            // perform cross-entropy optimization




            // form pseudo inverse
            boost::numeric::ublas::matrix<float> um;
            boost::numeric::ublas::matrix<float> sm;
            boost::numeric::ublas::matrix<float> vm;
            // total_q must be > total_models
            svd(model_intensities, um, sm, vm);

            // perform fit
            boost::numeric::ublas::matrix<float> a_star;
            pseudo_inverse(lmax, a_star, um, sm, vm);

            // score the model
//            sum = 0.0f;
//            for(unsigned int q =0; q < totalInWorkingSet; q++) {
//                diff = i_obs[q] - scalar * pIcalc[q];
//                sum += diff*diff*inv_variance[q];
//                residuals[q] = diff;
//            }
//
//            sum *= inv_n_s;
//            chi = sum;
//            dw = fabs(2 - calculateDurbinWatson(residuals.data(), totalInWorkingSet)); // should be near zero
//            sum += dw;
//
//            if (sum < bestScore){
//
//                bestScore = sum;
//                bestChi = chi;
//                bestBfactor = b_factor;
//                bestCx = cx;
//                bestScaleCoef = scalar;
//                bestDW = dw;
//
//                // do point transfer and do the scalar after finishing
//                for(unsigned int q =0; q < totalInWorkingSet; q++) {
//                    i_best[q] = scalar * pIcalc[q];
//                }
//            }

            b_factor += delta_b;
        }

        cx += deltaCX*contrastMatchScale;
    }
}

void Ensemble::assemblePartials(unsigned int lmax, AtomisticModel &model,
                           Waters &waterModel,
                           unsigned int totalInWorkingSet,
                           std::vector<float> & aP,
                           std::vector<float> & pPW,
                           std::vector<float> & pCX,
                           std::vector<float> & pCrossterm
                           ){

    const float * almPR = model.getModelPartialsReal();
    const float * almPI = model.getModelPartialsImag();
    const float * xlmR = model.getExcludedVolumePartialsReal();
    const float * xlmI = model.getExcludedVolumePartialsImag();

    const float * wlmR = waterModel.getModelPartialsReal();
    const float * wlmI = waterModel.getModelPartialsImag();

    float * pAP = aP.data();
    float * pApw = pPW.data();
    float * pAc = pCX.data();

    // cross-terms
    float * pAXW_cross_term = pCrossterm.data();

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
        pAP[q] = a_p_norm;
        pAc[q] = a_x_norm; // excluded volume
        pApw[q] = a_p_norm + a_w_norm + 2.0f*(aPWR_term + aPWI_term); // hydratedParticle term
        pAXW_cross_term[q] = -2.0f*(aXWR_term + aXWI_term + aPXR_term + aPXI_term); // invert phase

    } // end loop over q values
}

// if DW score is 2, there is no autocorrelation
float Ensemble::calculateDurbinWatson(float *residuals, unsigned int total) {

    float val, diff, first = residuals[0];
    float squared_sum = first*first;
    float squared_diff = 0.0;

    for(unsigned int i=1; i<total; i++){
        val = residuals[i];
        diff = (val - residuals[i-1]);
        squared_diff += diff*diff;
        squared_sum += val*val;
    }

    val = 2.0f - squared_diff/squared_sum;
    return val*val;
}

float Ensemble::populateIntensitiesMatrix(unsigned int totalInWorkingSet, float *qvalues,
                                                    std::vector<Amplitudes> models,
                                                    float cx,
                                                    float b_factor,
                                                    boost::numeric::ublas::matrix<float> & intensities) {

    float q_val;
    float exp_term, exp_2;

    unsigned int total_models = models.size();

    for(unsigned int q =0; q < totalInWorkingSet; q++) {

        q_val = qvalues[q]; // get qvalue to calculate b-factor
        exp_term = expf(-(q_val*q_val)*b_factor*inv16Pi2)*cx;
        exp_2 = exp_term*exp_term;

        // update matrix column elements (fixed q row)
        for(unsigned int n=0; n<total_models; n++){
            auto & model = models[n];
            intensities(q,n) = model.pAPWs[q] + model.pExVols[q]*exp_2 + exp_term*model.pACrossTerms[q];
        }
    }

    return 0;
}

/*
 * adapted from SVD-ublas by Volodymyr Kysenko (vksnk)
 */
void Ensemble::svd(boost::numeric::ublas::matrix < float >&A,
                   boost::numeric::ublas::matrix < float >&QQL,
                   boost::numeric::ublas::matrix < float >&sigmas, // singular values
                   boost::numeric::ublas::matrix < float >&QQR) {
    unsigned int row_num = A.size1();
    unsigned int col_num = A.size2();

    QQL.resize(row_num, row_num);
    sigmas.resize(row_num, col_num);
    QQR.resize(col_num, col_num);

    identity(QQL);
    identity(QQR);

    unsigned int to = std::min(row_num, col_num);

    for (unsigned int i = 0; i < to; i++) {
        householder(A, QQL, i, i, true);
        householder(A, QQR, i, i + 1, false);
    }

    boost::numeric::ublas::vector < float > d(to, 0.0f);
    boost::numeric::ublas::vector < float > s(to + 1, 0.0f);

    for (unsigned int i = 0; i < to; i++) {
        d(i) = A(i, i);
        if (i < (to - 1))
            s(i + 1) = A(i, i + 1);
    }

    svd_qr_shift(QQL, QQR, d, s);

    sigmas.clear();

    for (unsigned int i = 0; i < to; i++)
        sigmas(i, i) = d(i); // singular values
}


void Ensemble::householder(boost::numeric::ublas::matrix < float > & A,
                           boost::numeric::ublas::matrix < float > & QQ,
            unsigned int row_start, unsigned int col_start, bool column) {

    unsigned int size = column ? A.size1() : A.size2();
    unsigned int start = column ? row_start : col_start;

    if (start >= size)
        return;

    boost::numeric::ublas::vector < float > x(size);
    for (unsigned int i = 0; i < size; i++) {
        if (i < start) {
            x(i) = 0;
        } else {
            if (column)
                x(i) = A(i, col_start);
            else
                x(i) = A(row_start, i);
        }
    }

    float x_norm = norm(x);

    float alpha = (x(start) >= 0.0f ? 1.0f : -1.0f) * x_norm;

    boost::numeric::ublas::vector<float> v = x;

    v(start) += alpha;
    normalize(v);

    boost::numeric::ublas::matrix < float > Q;

    if (column) {
        for (unsigned int i = 0; i < A.size2(); i++) {
            float sum_Av = 0.0f;

            for (unsigned int j = 0; j < A.size1(); j++)
                sum_Av = sum_Av + (v(j) * A(j, i));
            for (unsigned int j = 0; j < A.size1(); j++)
                A(j, i) = A(j, i) - 2 * v(j) * sum_Av;
        }

        for (unsigned int i = 0; i < A.size1(); i++) {
            float sum_Qv = 0.0f;

            for (unsigned int j = 0; j < A.size1(); j++)
                sum_Qv = sum_Qv + (v(j) * QQ(i, j));
            for (unsigned int j = 0; j < A.size1(); j++)
                QQ(i, j) = QQ(i, j) - 2 * v(j) * sum_Qv;
        }

    } else {
        for (unsigned int i = 0; i < A.size1(); i++) {
            float sum_Av = 0.0f;

            for (unsigned int j = 0; j < A.size2(); j++)
                sum_Av = sum_Av + (v(j) * A(i, j));
            for (unsigned int j = 0; j < A.size2(); j++)
                A(i, j) = A(i, j) - 2 * v(j) * sum_Av;
        }

        for (unsigned int i = 0; i < A.size2(); i++) {
            float sum_Qv = 0.0f;

            for (unsigned int j = 0; j < A.size2(); j++)
                sum_Qv = sum_Qv + (v(j) * QQ(i, j));
            for (unsigned int j = 0; j < A.size2(); j++)
                QQ(i, j) = QQ(i, j) - 2 * v(j) * sum_Qv;
        }
    }
}


void Ensemble::svd_qr_shift(boost::numeric::ublas::matrix < float >&u,
                            boost::numeric::ublas::matrix < float >&v,
                            boost::numeric::ublas::vector < float >&q,
                            boost::numeric::ublas::vector < float >&e)
{
    int n = (int)q.size();
    int m = (int)u.size1();

    bool goto_test_conv = false;

    for (int k = n - 1; k >= 0; k--) {
        //std::cout << "U = " << u << std::endl;

        for (int iter = 0; iter < ITER_MAX; iter++) {
            // test for split
            int l;
            for (l = k; k >= 0; l--) {
                goto_test_conv = false;
                if (fabs(e[l]) <= EPS) {
                    // set it
                    goto_test_conv = true;
                    break;
                }

                if (fabs(q[l - 1]) <= EPS) {
                    // goto
                    break;
                }
            }

            if (!goto_test_conv) {
                float c = 0.0;
                float s = 1.0;

                int l1 = l - 1;

                for (int i = l; i <= k; i++) {
                    float f = s * e[i];
                    e[i] = c * e[i];

                    if (fabs(f) <= EPS) {
                        break;
                    }

                    float g = q[i];
                    float h = pythag(f, g);
                    q[i] = h;
                    c = g / h;
                    s = -f / h;

                    for (int j = 0; j < m; j++) {
                        float y = u(j, l1);
                        float z = u(j, i);
                        u(j, l1) = y * c + z * s;
                        u(j, i) = -y * s + z * c;
                    }
                }
            }

            float z = q[k];

            if (l == k) {
                if (z < 0.0f) {
                    q[k] = -z;

                    for(int j = 0; j < n; j++)
                        v(j, k) = -v(j,k);
                }

                break;
            }

            if (iter >= ITER_MAX - 1) {
                break;
            }

            float x = q[l];
            float y = q[k - 1];
            float g = e[k - 1];
            float h = e[k];
            float f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0f * h * y);

            g = pythag(f, 1.0);

            if (f < 0) {
                f = ((x - z) * (x + z) + h * (y / (f - g) - h)) / x;
            } else {
                f = ((x - z) * (x + z) + h * (y / (f + g) - h)) / x;
            }

            float c = 1.0;
            float s = 1.0;

            for (int i = l + 1; i <= k; i++) {
                g = e[i];
                y = q[i];
                h = s * g;
                g = c * g;
                float zz = pythag(f, h);
                e[i - 1] = zz;
                c = f / zz;
                s = h / zz;
                f = x * c + g * s;
                g = -x * s + g * c;
                h = y * s;
                y = y * c;

                for (int j = 0; j < n; j++) {
                    x = v(j, i - 1);
                    zz = v(j, i);
                    v(j, i - 1) = x * c + zz * s;
                    v(j, i) = -x * s + zz * c;
                }

                zz = pythag(f, h);
                q[i - 1] = zz;
                c = f / zz;
                s = h / zz;
                f = c * g + s * y;
                x = -s * g + c * y;

                for (unsigned int j = 0; j < m; j++) {
                    float yy = u(j, i - 1);
                    float zzz = u(j, i);
                    u(j, i - 1) = yy * c + zzz * s;
                    u(j, i) = -yy * s + zzz * c;
                }
            }
            e[l] = 0.0;
            e[k] = f;
            q[k] = x;
        }
    }
}

// new call for AMR, treat and prevent


void Ensemble::pseudo_inverse(
        unsigned int limit,
        boost::numeric::ublas::matrix < float >&A_star,
        boost::numeric::ublas::matrix < float >& um,
        boost::numeric::ublas::matrix < float >& sm, // singular values
        boost::numeric::ublas::matrix < float >& vm) {

    boost::numeric::ublas::matrix<float> inv_sm(sm.size1(), sm.size2(), 0.0f);
    // form pseudo inverse up to the lmax limit
    // not sure if I should be doing reduced SVD or regular?
    for(unsigned int i=0; i<limit; i++){
        inv_sm(i,i) = 1.0f/sm(i,i);
    }

    // A_star = vm*inv_sm_T*um_T
    A_star = boost::numeric::ublas::prod(boost::numeric::ublas::trans(inv_sm), boost::numeric::ublas::trans(um));
    A_star = boost::numeric::ublas::prod(vm, A_star);
}


void Ensemble::standardizeData(nsNNLS::vector * vec, nsNNLS::vector * out_vec){

    unsigned int window = 31;
    double sum=0.0;
    out_vec->setSize(vec->length());

    for(int i=0; i<window; i++){
        sum += vec->get(vec->length() - window + i);
    }

    standardizedMin = sum/(float) window;//nonData.getMinY();
    standardizedScale = std::abs(vec->get(3) - standardizedMin);

    double invstdev = 1.0/standardizedScale;

    for(int r=0; r<vec->length(); r++){
        out_vec->set(r, (vec->get(r) - standardizedMin)*invstdev);
    }

}


float Ensemble::crossEntropyOptimization(std::vector<float> & i_obs_over_variance,
                                        std::vector<float> & inv_variance,
                                        boost::numeric::ublas::matrix<float> & coeffs,
                                        boost::numeric::ublas::matrix<float> & model_iqs,
                                        boost::numeric::ublas::matrix<float> & i_obs){

    std::random_device rd;
    std::mt19937 gen(rd());

    unsigned int totalInWorkingSet = model_iqs.size1();
    float inv_trial_size = 1.0f/(float)trial_size;

    auto topN = (unsigned int)(frac*totalPerRound);

    unsigned int totalModels = model_iqs.size2();
    std::vector<unsigned int> model_counts(totalModels); // should be one for every model
    boost::numeric::ublas::matrix<float> b_vec(totalModels,1);
    boost::numeric::ublas::matrix<float> i_calc(model_iqs.size1(),1);

    std::uniform_real_distribution<float> distributions(0,1.0);
    std::uniform_int_distribution<unsigned int> indices(0, totalModels-1);
    std::vector<Probability> probabilities(totalModels);

    Probability * pProb;

    std::vector<Trial> topTrials;
    std::vector<unsigned int> trial(trial_size);

    topTrials.reserve(topN);
    topTrials.resize(topN);

    const unsigned int last = topN-1;
    float inv_topN = 1.0f/(float)topN;
    float scalar, score, chi, dw, diff, i_of_q;

    std::vector<float> residuals(totalInWorkingSet);
    std::vector<float> selections(totalModels, 0.0);
    std::vector<float> selectionsSquared(totalModels, 0.0);

    float inv_n_s = 1.0f/(float)(totalInWorkingSet-2);
    float best_score = FLT_MAX;

    for(unsigned int t = 0; t<totalTrials; t++){

        unsigned int topAdded=0;

        for(unsigned int n=0; n<totalPerRound; n++) {

            b_vec.resize(totalModels,1);
            std::fill(model_counts.begin(), model_counts.end(), 0);

            // sample distribution
            for (unsigned ts=0; ts<trial_size; ts++) {
                bool doIt = true;
                while(doIt){
                    unsigned int index = indices(gen);
                    if (distributions(gen) <= probabilities[index].value){
                        model_counts[index] += 1;
                        doIt = false;
                        break;
                    }
                }
            }

            // convert to weights
            for (unsigned m=0; m<totalModels; m++) {
                b_vec(m,0) = (float)model_counts[m] * inv_trial_size;
            }

            i_calc = boost::numeric::ublas::prod(model_iqs, b_vec);
            // perform fit
            float c_obs_calc = 0.0f;
            float calc_var = 0.0f;

            for(unsigned int q =0; q < totalInWorkingSet; q++) {
                i_of_q = i_calc(q,0);
                c_obs_calc += i_of_q*i_obs_over_variance[q];
                calc_var += i_of_q*i_of_q*inv_variance[q];
            }

            scalar = c_obs_calc/calc_var;

            score = 0.0f;
            for(unsigned int q =0; q < totalInWorkingSet; q++) {
                diff = i_obs(q,0) - scalar * i_calc(q,0);
                score += diff*diff*inv_variance[q];
                residuals[q] = diff;
            }

            score *= inv_n_s;
            chi = score;
            dw = fabs(2 - calculateDurbinWatson(residuals.data(), totalInWorkingSet)); // should be near zero
            score += 0.1*dw; // could be weighted?

            //update best list
            if (topAdded < topN){
                Trial * pTrial = &topTrials[topAdded];
                pTrial->value = score;
                pTrial->model_weights.swap(b_vec);
                pTrial->chi = chi;
                pTrial->dw = dw;
                topAdded++;
                std::sort(topTrials.begin(), topTrials.begin()+topAdded);
            } else if (score < topTrials[last].value){
                /*
                 * replace last entry and sort
                 */
                Trial * pTrial = &topTrials[last];
                pTrial->value = score;
                pTrial->model_weights.swap(b_vec);
                pTrial->chi = chi;
                pTrial->dw = dw;
                std::sort(topTrials.begin(), topTrials.end());
            }
        }

        /*
         * for each bead in topN, count occurrences
         */
        std::fill(selections.begin(), selections.end(), 0.0);
        float bval;
        for(auto & trial : topTrials){
            for(unsigned int b=0; b<totalModels; b++){
                selections[b] += trial.model_weights(b,0);
            }
        }

        float oldprob;
        for(unsigned int b=0; b<totalModels; b++){
            pProb = &probabilities[b];
            oldprob = (1.0f - updateAlpha)*pProb->value;
            (pProb->value) = updateAlpha*selections[b]*inv_topN + oldprob;
        }


        if (topTrials[0].value < best_score){
            best_score = topTrials[0].value;
            std::cout << " Best Score " << t << " " << best_score << " " << topTrials[0].dw << " " << topTrials[0].model_weights << std::endl;
            for(unsigned int i=0; i<totalModels; i++){
                coeffs(i,0) = topTrials[0].model_weights(i,0);
            }
        } // stopping criteria?
    }

    std::vector<float> averages(totalModels,0.0);
    for(unsigned int i=0; i<5; i++){
        for(unsigned int b=0; b<totalModels; b++) {
            averages[b] += topTrials[i].model_weights(b, 0);
        }
    }

    for(unsigned int b=0; b<totalModels; b++) {
        std::cout << b << " " << averages[b]/5 << std::endl;
    }

    return best_score;
}