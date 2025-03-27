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

/*
 * CE optimization method akin to knapsack problem
 * adjustable parameters :: totalTrials
 * coeffs are output values for model weights
 */
float Ensemble::crossEntropyOptimization(std::vector<float> & i_obs_over_variance,
                                        std::vector<float> & inv_variance,
                                        boost::numeric::ublas::matrix<float> & coeffs,
                                        boost::numeric::ublas::matrix<float> & model_iqs,
                                         std::vector<float> & qvalues,
                                         boost::numeric::ublas::matrix<float> & i_calc_all,
                                        boost::numeric::ublas::matrix<float> & i_obs){

    std::random_device rd;
    std::mt19937 gen(rd());

    unsigned int totalInWorkingSet = model_iqs.size1();
    float inv_ensemble_size = 1.0f/(float)ensemble_size;

    auto topN = (unsigned int)(rarity_param * totalPerRound);

    unsigned int totalModels = model_iqs.size2();
    std::vector<unsigned int> model_counts(totalModels); // should be one for every model
    boost::numeric::ublas::matrix<float> b_vec(totalModels,1);
    boost::numeric::ublas::matrix<float> i_calc(model_iqs.size1(),1);

    std::uniform_real_distribution<float> distributions(0,1.0);
    std::uniform_int_distribution<unsigned int> indices(0, totalModels-1);

    std::vector<Probability> probabilities(totalModels); // each element initialized as 1/2
    Probability * pProb;

    std::vector<Trial> topTrials;
    std::vector<unsigned int> trial(ensemble_size);

    topTrials.reserve(topN);
    topTrials.resize(topN);

    const unsigned int last = topN-1;
    float inv_topN = 1.0f/(float)topN;
    float scalar, score, chi, dw, diff, i_of_q, temp_residuals;

    std::vector<float> residuals(totalInWorkingSet);
    std::vector<float> selections(totalModels, 0.0);
    std::vector<float> selectionsSquared(totalModels, 0.0);
    std::vector<unsigned int> selectedIndices(ensemble_size);


    float inv_n_s = 1.0f/(float)(totalInWorkingSet-2);
    float best_score = FLT_MAX, best_scale_factor = 1;

    char buffer [50];

    logger("Total Models in Search", std::to_string(totalModels));
    logger("ENSEMBLE SIZE", std::to_string(ensemble_size));
    logger("Top N", std::to_string(topN));
    logger("Total Trials", std::to_string(totalTrials));
    logger("Total Per Trial", std::to_string(totalPerRound));
    logger("------------------------------","------------------------------");

    for (unsigned m=0; m<totalModels; m++) {
        b_vec(m,0) = 1.0/(float)totalModels;
    }

    i_calc_all = boost::numeric::ublas::prod(model_iqs, b_vec);
    // calculate scalar
    scalar = calculate_scalar(totalInWorkingSet, i_calc_all, i_obs_over_variance, inv_variance);
    // print out the entire ensemble composite scattering curve
    //i_calc_all *= scalar;

    std::string nameOf = "all.dat";
    FILE * pFile = fopen(nameOf.c_str(), "w");
    for(unsigned int q =0; q < totalInWorkingSet; q++) {
        float q_val = qvalues[q];
        fprintf(pFile, "%5i %.5E %.5E\n", q+1,  qvalues[q], i_calc_all(q,0));
    }
    fclose(pFile);

    for(unsigned int t = 0; t<totalTrials; t++){ // upper limit - can terminate early if changes < epsilon

        unsigned int topAdded=0;

        for(unsigned int n=0; n<totalPerRound; n++) {

            b_vec.resize(totalModels,1);
            selectedIndices.resize(ensemble_size);

            std::fill(model_counts.begin(), model_counts.end(), 0); // should be changed to map<key,value> pair

            // pick members of the ensemble, can do with replacement meaning same members can be picked multiple times
            unsigned int ens = 0;
            for (unsigned ts=0; ts < ensemble_size; ts++) {
                bool doIt = true;
                while(doIt){
                    unsigned int index = indices(gen);
                    if (distributions(gen) <= probabilities[index].value){
                        selectedIndices[ens] = index;
                        model_counts[index] += 1;
                        ens++;
                        doIt = false;
                        break;
                    }
                }
            }

            std::sort(selectedIndices.begin(), selectedIndices.end());

            // convert to weights
            for (unsigned m=0; m<totalModels; m++) {
                b_vec(m,0) = (float)model_counts[m] * inv_ensemble_size;
            }

            i_calc = boost::numeric::ublas::prod(model_iqs, b_vec);
            // calculate scalar
            scalar = calculate_scalar(totalInWorkingSet, i_calc, i_obs_over_variance, inv_variance);

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
                pTrial->indices.swap(selectedIndices);
                pTrial->chi = chi;
                pTrial->dw = dw;
                pTrial->scale_factor = scalar;
                topAdded++;
                std::sort(topTrials.begin(), topTrials.begin()+topAdded);
            } else if (score < topTrials[last].value){
                // is the new trial difference or similar to existing trial
                // perform cosine similarity for selected indices
                // if similar and lower energy, replace the configuration
                // if different than all of them, and lower energy, replace the last model and sort the list
                /*
                 * replace last entry and sort
                 */
                Trial * pTrial = &topTrials[last];
                pTrial->value = score;
                pTrial->model_weights.swap(b_vec);
                pTrial->indices.swap(selectedIndices);
                pTrial->chi = chi;
                pTrial->dw = dw;
                pTrial->scale_factor = scalar;
                std::sort(topTrials.begin(), topTrials.end());
            }
        }

        /*
         * for each in topN, count occurrences
         */
        std::fill(selections.begin(), selections.end(), 0.0);
        float bval;
        for(auto & trial : topTrials){
            for(unsigned int b=0; b<totalModels; b++){
                selections[b] += trial.model_weights(b,0); // weights are normalized to ensemble size
            }
        }

        // updating is naive, could do something more sophisticated using expectation-maximization
        if (t < 13){ // smooth updating model
            for(unsigned int b=0; b<totalModels; b++){
                pProb = &probabilities[b];
                (pProb->value) = updateAlpha*selections[b]*inv_topN + (1.0f - updateAlpha)*pProb->value;
            }
        } else { // switch to coarse
            for(unsigned int b=0; b<totalModels; b++){
                probabilities[b].value = selections[b]*inv_topN;
            }
        }

        // copy coefficients of the best performances
        if (topTrials[0].value < best_score){
            best_score = topTrials[0].value;
            best_scale_factor = topTrials[0].scale_factor;
            for(unsigned int i=0; i<totalModels; i++){
                coeffs(i,0) = topTrials[0].model_weights(i,0)*best_scale_factor;
            }
        }

        sprintf (buffer, "TRIAL %3i", t); // change the number to a string for name
        std::string strial = buffer;

        sprintf (buffer, "BEST SCORE %7.4f DELTA %.4f", best_score, std::abs(best_score-topTrials[last].value)); // change the number to a string for name
        logger(strial, buffer);
    }

    logger("------------------------------","------------------------------");
    // find smallest probability in set greater than zero
    float smallest = FLT_MAX;
    for(unsigned int b=0; b<totalModels; b++) {
        float * p = & probabilities[b].value;
        if (*p < smallest && *p > 0) {
            smallest = *p;
        }
    }

    // get the top 4 models and fit them to the data
    std::vector<unsigned int> top4(4);
    float maxp = 0;
    for(unsigned int b=0; b<totalModels; b++) {
        float * p = & probabilities[b].value;
        if (*p > maxp ) {
            maxp = *p;
            top4[0] = b;
        }
    }
    float maxlimit = maxp;
    maxp = 0;
    for(unsigned int b=0; b<totalModels; b++) {
        float * p = & probabilities[b].value;
        if (*p > maxp && maxp <= maxlimit) {
            maxp = *p;
            top4[1] = b;
        }
    }

    maxlimit = maxp;
    maxp = 0;
    for(unsigned int b=0; b<totalModels; b++) {
        float * p = & probabilities[b].value;
        if (*p > maxp && maxp <= maxlimit) {
            maxp = *p;
            top4[2] = b;
        }
    }

    maxlimit = maxp;
    maxp = 0;
    for(unsigned int b=0; b<totalModels; b++) {
        float * p = & probabilities[b].value;
        if (*p > maxp && maxp <= maxlimit) {
            maxp = *p;
            top4[3] = b;
        }
    }

    // perform linear fits of top four identified conformers
//    boost::numeric::ublas::matrix<float> model_iqs_4(model_iqs.size1(), 4);
//    unsigned int col = 0;
//    for (auto & ind : top4) {
//        // combine and add noise
//        for(unsigned int q=0; q<totalInWorkingSet; q++){
//            model_iqs_4(q, col) = model_iqs(q, ind);
//        }
//        col+=1;
//    }

    std::vector<std::vector<unsigned int>> comb;
    comb.emplace_back(std::vector<unsigned int>(1));
    comb.emplace_back(std::vector<unsigned int>(2));
    comb.emplace_back(std::vector<unsigned int>(3));
    comb.emplace_back(std::vector<unsigned int>(4));

    combinatorial(
                  inv_variance,
                  top4,
                  model_iqs,
                  i_obs,
                  totalInWorkingSet,
                  4);

    nameOf = "selected_ensemble_summary.txt";
    pFile = fopen(nameOf.c_str(), "w");

    unsigned int totalSelected = 0;
    for(unsigned int b=0; b<totalModels; b++) {
        if (probabilities[b].value > 0) {
            unsigned int bbt = std::ceil(probabilities[b].value/smallest);
            Param *pr = &params[b];
            totalSelected+=1;
            selected_files.push_back(pr->filename);
            for(unsigned int j=0; j<bbt; j++){
                fprintf(pFile, "%5i %s %7.5f %3i %8.4f %9.4f\n", b, pr->filename.c_str(), probabilities[b].value, bbt, pr->rg, pr->dmax);
            }
        }
    }
    fclose(pFile);


    nameOf = "selected_ensemble_for_alignment.txt";
    pFile = fopen(nameOf.c_str(), "w");

    for(unsigned int b=0; b<totalModels; b++) {
        if (probabilities[b].value > 0) {
            unsigned int bbt = std::ceil(probabilities[b].value/smallest);
            Param *pr = &params[b];
            totalSelected+=1;
            selected_files.push_back(pr->filename);
            fprintf(pFile, "%s %5i \n", pr->filename.c_str(), bbt);
        }
    }
    fclose(pFile);


    logger("Total Unique Models", std::to_string(totalSelected));

    nameOf = "all_summary.txt";
    pFile = fopen(nameOf.c_str(), "w");

    for(unsigned int b=0; b<totalModels; b++) {
        Param *pr = &params[b];
        fprintf(pFile, "%5i %s %7.5f %8.4f\n", b, pr->filename.c_str(), pr->rg, pr->dmax);
    }
    fclose(pFile);

    return best_score;
}

/*
 *
 * model_iqs need to be normalized!  Meaning, the values in the column should scale from 1.0 as max
 *
 * n is the total number of models in the search
 * k is the linear combination of models to use
 * k < n, maxk < n
 *
 * need an Lmax check on the number of columns
 * unsigned int limitSm = (lmax < trial_models.size2()) ? lmax : trial_models.size2();
 *
 * using SVD, scaling constants that are negative are nonsensical as a solution (can't have negative concentrations)
 * We will only report solutions that are all positive and infer that poor solutions not reported have negative scaling constants
 * NNLS was tried, results take longer but similar.
 *
 */
void Ensemble::combinatorial (
                              std::vector<float> &inv_variance,
                              std::vector<unsigned int> selected_indices_of_model_iqs,
                              boost::numeric::ublas::matrix<float> &model_iqs,
                              boost::numeric::ublas::matrix<float> &i_obs,
                              unsigned int totalQInWorkingSet,
                              int maxk) {

    // must refer to last index of total models
    unsigned int n_models = selected_indices_of_model_iqs.size()-1;

    if (maxk >= n_models){ // throw exception

    }

    logger("Total Models", std::to_string(n_models));
    logger("Starting", "pairwise SVD search");

    boost::numeric::ublas::matrix<float> um;
    boost::numeric::ublas::matrix<float> sm;
    boost::numeric::ublas::matrix<float> vm;
    boost::numeric::ublas::matrix<float> astar;

    std::vector<float> residuals(totalQInWorkingSet);
    float diff, score, chi, dw;
    float inv_n_s = 1.0f/(float)(totalQInWorkingSet-2);

    std::vector<Trial> topTrials;
    topTrials.reserve(maxk);
    topTrials.resize(maxk);

    for (unsigned int k=1; k<=maxk; k++){

        boost::numeric::ublas::matrix<float> trial_models(totalQInWorkingSet,k);

        float best_score = FLT_MAX;

        auto & bestTrial = topTrials[k-1];
        // initialize values in numbers
        std::vector<unsigned int> numbers(k);
        for (unsigned int i = 0; i < k; ++i)
            numbers.at(i) = i;

        do { // iterate through all combinations of k
            /*
             * load the trial_models matrix
             */
            for (unsigned int i=0; i<k; i++){
                unsigned int model_index = selected_indices_of_model_iqs[numbers[i]];
                for(unsigned int q=0; q < totalQInWorkingSet; q++){
                    trial_models(q, i) = model_iqs(q, model_index) ;
                }
            }

            boost::numeric::ublas::matrix<float> a_matrix(trial_models); // make copy of matrix
            this->svd(trial_models, um, sm, vm);
            // lmax might be 23 and trial_models.size2 might be 2, so looking for the best pair
            // unsigned int limitSm = (lmax < trial_models.size2()) ? lmax : trial_models.size2();
            unsigned int limitSm = trial_models.size2();
            this->pseudo_inverse(limitSm, astar, um, sm, vm);

            boost::numeric::ublas::matrix<float> constants = boost::numeric::ublas::prod(astar, i_obs);

            // all scaling constants must be greater than zero
            bool skip = true;
            for(int c=0; c<constants.size1(); c++){
                if (constants(c,0) < 0){
                    skip = false;
                }
            }

            if (skip){ // score the answer
                boost::numeric::ublas::matrix<float> i_calc = boost::numeric::ublas::prod(a_matrix, constants);
                // calculate residuals
                for(unsigned int q =0; q < totalQInWorkingSet; q++) {
                    diff = i_obs(q,0) - i_calc(q,0);
                    score += diff*diff*inv_variance[q];
                    residuals[q] = diff;
                }

                score *= inv_n_s;

                chi = score;
                dw = fabs(2 - calculateDurbinWatson(residuals.data(), totalQInWorkingSet)); // should be near zero
                score += 0.1*dw;

                if (score < best_score){
                    best_score = score;
                    std::cout << "Best " << best_score << " " << k << std::endl;

                    float sumIt = 0;
                    for(int c=0; c<constants.size1(); c++){
                        sumIt += constants(c,0);
                    }

                    for(int c=0; c<constants.size1(); c++){
                        std::cout << " c" << numbers[c] << " : " << constants(c,0)/sumIt << " ";
                    }
                    std::cout << "\n";

                    bestTrial.chi  = chi; // record each top trial for combination
                    bestTrial.dw = dw;
                    bestTrial.indices = numbers;
                    bestTrial.model_weights.swap(constants);
                }
            }

        } while (this->next_combination((int)n_models, (int)k, numbers));

        // write out the best trials, report the best fit and the score difference between first and last
//        std::cout << "k " << k << " " <<  topTrials[0].chi << " " << topTrials[0].dw << std::endl;
//        auto & besty = topTrials[0];
//        for(unsigned int ind =0; ind < k; ind++){
//            std::cout << ind << " "  << filenames[besty.indices[ind]] << " " << besty.model_weights(ind,0) << std::endl;
//        }

//        for(unsigned q=0 ; q<totalQInWorkingSet; q++){
//            std::cout << q <<  " ";
//            for(unsigned int ind =0; ind < k; ind++){
//                std::cout << model_iqs(q, besty.indices[ind])*besty.model_weights(ind,0) << " ";
//            }
//            std::cout << i_obs(q,0) << std::endl;
//        }


        // write out best performing models.
        for (auto & best : topTrials){

        }


    }
}