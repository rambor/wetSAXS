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

#ifndef WETSAXS_ENSEMBLE_H
#define WETSAXS_ENSEMBLE_H


#include <sastools/FileClass.h>
#include <boost/filesystem/path.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <vector>
#include <string>
#include <sastools/IofQData.h>
#include "AtomisticModel.h"
#include "Waters.h"
#include "../fsrc/nnls.h"


class Ensemble {

    struct Amplitudes {
        unsigned int index;
        std::string name;
        std::vector<float> aPs;
        std::vector<float> pAPWs;
        std::vector<float> pACrossTerms;
        std::vector<float> pExVols;

        Amplitudes(unsigned int index, std::string name) : index(index), name(name) {}

        void addValues(std::vector<float> aps, std::vector<float> aPWs, std::vector<float> pACs, std::vector<float> crossterm){
            aPs = std::move(aps);
            pAPWs = std::move(aPWs);
            pExVols = std::move(pACs);
            pACrossTerms = std::move(crossterm);
        }
    };

    struct Param {
        std::string filename;
        float rg;
        float dmax;

        Param(std::string s, float f, float d) : filename(std::move(s)), rg(f), dmax(d) {}

    };

    public:
    struct Score {
        unsigned int row, col;
        float score;

        Score(unsigned int row, unsigned int col, float score) : row(row), col(col), score(score) {}
    };


    struct Probability {
        float value;
        Probability() : value(0.50) {} // initialize
    };

    struct Trial{
        double value;
        double chi, dw, scale_factor;
        boost::numeric::ublas::matrix<float> model_weights;
        std::vector<unsigned int> indices;
        Trial() = default;

        Trial(double val, boost::numeric::ublas::matrix<float> & vec) : value(val) {
            model_weights = std::move(vec);
        }

        bool operator<(const Trial & a) const{
            return value < a.value;
        }


        /*
         * cosine similarity
         * both arrays should be sorted
         */
        float score(std::vector<unsigned int> & test){
            unsigned int * pvalA, * pvalB;
            unsigned sumA = 0;
            unsigned sumB = 0;
            unsigned sumAB = 0;
            for(unsigned int ind=0; ind<indices.size(); ind++){
                pvalA = &indices[ind];
                pvalB = &test[ind];
                sumAB += *pvalA* *pvalB;
                sumA += *pvalA* *pvalA;
                sumB += *pvalB* *pvalB;
            }

            return (float)sumAB/(sqrtf(float(sumA)) * sqrtf(float(sumB)));
        }
    };

    std::vector<std::string> filenames;
    std::vector<Param> params;
    std::string directory;
    const float inv16Pi2 = 1.0/(16.0*M_PI_2);
    const float EPS = 0.00001;
    const int ITER_MAX = 50;
    float standardizedMin;
    float standardizedScale;

    const float updateAlpha=0.33333367f;
    float rarity_param = 0.0171;
    unsigned int totalTrials = 41;
    unsigned int ensemble_size = 21;
    unsigned int totalPerRound = 3331;
    std::vector<std::string> selected_files;

public:
    Ensemble(std::vector<std::string> filenames, std::string directory);

    unsigned int getTotalFiles(){ return filenames.size();}

    void search(unsigned int lmax, IofQData * file, bool forceRNA, bool keepWaters = false, int oversampling=5);

    void assemblePartials(unsigned int lmax, AtomisticModel &model, Waters &waterModel, unsigned int totalInWorkingSet,
                          std::vector<float> & aP,
                          std::vector<float> & pPW,
                          std::vector<float> & pCX,
                          std::vector<float> & pCrossterm);

    std::string getDirectory(){ return directory;}

    float calculateDurbinWatson(float *residuals, unsigned int total);


    /*
     *
     * find the pairs with the largest summed difference
     * for 2, find largest difference
     * for 3, find the largest summed difference (maximising area of triangle)
     * for 4, find the largest summed differences
     *
     * subset_len is the ensemble size
     * maxd is the maximum calculated difference between initialized as 0
     *
     */
    void combinations(std::vector<unsigned int> &indices, std::vector<Score> &scores,
                      int subset_len,
                      int startPosition,
                      std::vector<unsigned int> & result,
                      std::vector<unsigned int> &best,
                      float & maxd){

        if (subset_len == 0){
            unsigned int index, total = indices.size();
            float sum = 0;
            // calculate the distance sum and return
            // results has an array of sorted indices
            for(unsigned int row = 0; row<(result.size()-1); row++){
                unsigned int rr = result[row];
                unsigned int col = row+1;
                for(; col<result.size(); col++){
                    index = rr * total - (rr * (rr + 1) / 2) + result[col] - (rr + 1);
                    sum += scores[index].score;
                }
            }

            if (sum > maxd){
                maxd = sum;
                std::copy(result.begin(), result.end(), best.begin());
            }

            return;
        }

        for (int i = startPosition; i <= indices.size()-subset_len; i++){
            result[result.size() - subset_len] = indices[i];
            combinations(indices, scores, subset_len-1, i+1, result, best, maxd);
        }
    }


    float populateIntensitiesMatrix(unsigned int totalQ,
                                              float *qvalues,
                                              std::vector<Amplitudes> models,
                                              float cx,
                                              float bfactor,
                                              boost::numeric::ublas::matrix<float> & intensities);

    /*
    * adapted from SVD-ublas by Volodymyr Kysenko (vksnk)
    */
    float pythag(float a, float b) {
        float absa = fabs(a);
        float absb = fabs(b);

        if (absa > absb) {
            return absa * sqrtf(1.0f + simde_math_powf(absb / absa, 2));
        } else {
            return absb * sqrtf(1.0f + simde_math_powf(absa / absb, 2));
        }
    }

    /*
     * adapted from SVD-ublas by Volodymyr Kysenko (vksnk)
     */
    void identity(boost::numeric::ublas::matrix < float >&m){
        for (unsigned int i = 0; i < m.size1(); i++)
            for (unsigned int j = 0; j < m.size2(); j++)
                m(i, j) = (i == j) ? 1.0f : 0.0f;
    }

    /*
     * adapted from SVD-ublas by Volodymyr Kysenko (vksnk)
     */
    float norm(boost::numeric::ublas::vector < float >&x){
        float x_norm = 0.0;

        for (unsigned int i = 0; i < x.size(); i++)
            x_norm += simde_math_powf(x(i),2);

        return sqrtf(x_norm);
    }

    /*
     * adapted from SVD-ublas by Volodymyr Kysenko (vksnk)
     */
    void normalize(boost::numeric::ublas::vector < float >&x){
        float invx_norm = 1.0f/norm(x);
        for (unsigned int i = 0; i < x.size(); i++) {
            x(i) *= invx_norm;
        }
    }

    void svd(boost::numeric::ublas::matrix<float> &A,
             boost::numeric::ublas::matrix<float> &QQL,
             boost::numeric::ublas::matrix<float> &QQW,
             boost::numeric::ublas::matrix<float> &QQR);


    void householder(boost::numeric::ublas::matrix<float> &A,
                     boost::numeric::ublas::matrix<float> &QQ,
                     unsigned int row_start,
                     unsigned int col_start,
                     bool column);

    void svd_qr_shift(boost::numeric::ublas::matrix<float> &u,
                      boost::numeric::ublas::matrix<float> &v,
                      boost::numeric::ublas::vector<float> &q,
                      boost::numeric::ublas::vector<float> &e);


    void pseudo_inverse(unsigned int lmax,
                        boost::numeric::ublas::matrix<float> &A_star,
                        boost::numeric::ublas::matrix<float> &um,
                        boost::numeric::ublas::matrix<float> &sm,
                        boost::numeric::ublas::matrix<float> &vm);

//    void nnls(boost::numeric::ublas::matrix<float> &y_m, boost::numeric::ublas::matrix<float> &x_m,
//              boost::numeric::ublas::matrix<float> &A_m);

    unsigned int get_max_index(boost::numeric::ublas::matrix<float> &vector);

    void standardizeData(nsNNLS::vector *vec, nsNNLS::vector *out_vec);

    float getStandardizedScale(){ return standardizedScale;}
    float getStandardizedMin(){ return standardizedMin;}


    /**
     *
     * @param iofqdata
     * @param cx
     * @param bfactor
     *
     * dmax cutoff
     */
    void ce_search(IofQData & iofqdata, float cx, float bfactor, bool dmaxFilter){

        iofqdata.makeWorkingSet(5); // meaningless number
        unsigned int totalInWorkingSet = iofqdata.getTotalInWorkingSet();
        auto qvalues = iofqdata.getWorkingSetQvalues();
        const std::vector<Datum> & workingSet = iofqdata.getWorkingSet();

        unsigned lmax = (unsigned int)((iofqdata.getQmax()*iofqdata.getDmax())/SIMDE_MATH_PI) + 1;
        SASTOOLS_UTILS_H::logger("LMAX", std::to_string(lmax));

        std::vector<float> i_obs_over_variance(totalInWorkingSet);
        std::vector<float> inv_variance(totalInWorkingSet);
        boost::numeric::ublas::matrix<float> y_vector(totalInWorkingSet, 1);

        for(unsigned int q =0; q < totalInWorkingSet; q++) {
            const Datum *pW = &workingSet[q];
            i_obs_over_variance[q] = pW->getI() * pW->getInvVar();
            inv_variance[q] = pW->getInvVar();
            y_vector(q, 0) = pW->getI();
        }

        boost::numeric::ublas::matrix<float> coeffs(filenames.size(),1);
        boost::numeric::ublas::matrix<float> models(totalInWorkingSet,filenames.size());
        unsigned int mdl_cnt=0;

        for (auto & filename : filenames){

            try{
                AtomisticModel model = AtomisticModel(filename, false, false);

                model.calculatePartialAmplitudes(lmax, totalInWorkingSet, qvalues, false);

                if (model.getPDBModel().getTotalResidues() < 1){
                    throw ("skipping..." + filename);
                }

                if (dmaxFilter && model.getDmax() > iofqdata.getDmax()){
                    throw(filename.append(" dmax exceeded ").append(std::to_string(model.getDmax())));
                }

                // create water model
                Waters waters = Waters();
                waters.hydrateAtomisticModel(model); // add waters to PDB
                waters.createSphericalCoordinateOfHydration();
                waters.calculatePartialAmplitudes(lmax, totalInWorkingSet, qvalues);

                std::vector<float> aPs = std::vector<float>(totalInWorkingSet);
                std::vector<float> aPWs = std::vector<float>(totalInWorkingSet);
                std::vector<float> aCs = std::vector<float>(totalInWorkingSet);
                std::vector<float> aXW_cross_term = std::vector<float>(totalInWorkingSet);

                // assemble the partials and
                this->assemblePartials(lmax, model, waters, totalInWorkingSet, aPs, aPWs, aCs, aXW_cross_term);

                // combine
                for(unsigned int q=0; q<totalInWorkingSet; q++){
                    float q_val = qvalues[q]; // get qvalue to calculate b-factor
                    float exp_term = expf(-(q_val*q_val)*bfactor*inv16Pi2)*cx;
                    // from assemblePartials in Ensemble
                    //        pAP[q] = a_p_norm;
                    //        pAc[q] = a_x_norm; // excluded volume
                    //        pApw[q] = a_p_norm + a_w_norm + 2.0f*(aPWR_term + aPWI_term); // hydratedParticle term
                    //        pAXW_cross_term[q] = -2.0f*(aXWR_term + aXWI_term + aPXR_term + aPXI_term); // invert phase
                    models(q, mdl_cnt) = aPWs[q] + aCs[q]*exp_term*exp_term + exp_term*aXW_cross_term[q];
                }

                params.emplace_back(Param(filename, model.calculateRg(), model.getDmax()));

                mdl_cnt++;
            } catch (std::string & fl){
                SASTOOLS_UTILS_H::logger("Improper PDB file", fl);
            }
        }

        if (dmaxFilter){
            coeffs.resize(mdl_cnt,1, true);
            models.resize(totalInWorkingSet,mdl_cnt, true);
        }

        SASTOOLS_UTILS_H::logger("Starting", "CE SEARCH");
        rarity_param = 0.071371;
        ensemble_size = 37;
        totalPerRound = 7*models.size2(); // should be optimized but for now will do 3xtimes

        boost::numeric::ublas::matrix<float> i_calc_all(models.size1(),1);

        float val = this->crossEntropyOptimization(i_obs_over_variance,
                                                   inv_variance,
                                                   coeffs,
                                                   models,
                                                   qvalues,
                                                   i_calc_all,
                                                   y_vector);

        // calculate for all q-values
        //
        // iofqdata.setAllDataToWorkingSet();
        // coeffs relates only to the
        auto i_calc = boost::numeric::ublas::prod(models, coeffs);

        std::string nameOf = "top_fit.dat";
        FILE * pFile = fopen(nameOf.c_str(), "w");
        for(unsigned int q =0; q < totalInWorkingSet; q++) {
            float p = y_vector(q,0);
            float c = i_calc(q,0);
            float q_val = qvalues[q];
            fprintf(pFile, "%5i %.6E %.5E %.5E %.5E %.5E\n", q+1, q_val, p, c, (p-c), i_calc_all(q,0));
        }
        fclose(pFile);

//        iofqdata.setAllDataToWorkingSet();
//
//        for(auto filename : selected_files){
//            AtomisticModel model = AtomisticModel(filename, false, false);
//
//            model.calculatePartialAmplitudes(lmax, totalInWorkingSet, qvalues, false);
//
//            // create water model
//            Waters waters = Waters();
//            waters.hydrateAtomisticModel(model); // add waters to PDB
//            waters.createSphericalCoordinateOfHydration();
//            waters.calculatePartialAmplitudes(lmax, totalInWorkingSet, qvalues);
//
//            std::vector<float> aPs = std::vector<float>(totalInWorkingSet);
//            std::vector<float> aPWs = std::vector<float>(totalInWorkingSet);
//            std::vector<float> aCs = std::vector<float>(totalInWorkingSet);
//            std::vector<float> aXW_cross_term = std::vector<float>(totalInWorkingSet);
//
//            // assemble the partials and
//            this->assemblePartials(lmax, model, waters, totalInWorkingSet, aPs, aPWs, aCs, aXW_cross_term);
//
//            // combine and add noise
//            for(unsigned int q=0; q<totalInWorkingSet; q++){
//                float q_val = qvalues[q]; // get qvalue to calculate b-factor
//                float exp_term = expf(-(q_val*q_val)*bfactor*inv16Pi2)*cx;
//                // from assemblePartials in Ensemble
//                //        pAP[q] = a_p_norm;
//                //        pAc[q] = a_x_norm; // excluded volume
//                //        pApw[q] = a_p_norm + a_w_norm + 2.0f*(aPWR_term + aPWI_term); // hydratedParticle term
//                //        pAXW_cross_term[q] = -2.0f*(aXWR_term + aXWI_term + aPXR_term + aPXI_term); // invert phase
//                models(q, mdl_cnt) = aPWs[q] + aCs[q]*exp_term*exp_term + exp_term*aXW_cross_term[q];
//            }
//            params.emplace_back(Param(filename, model.calculateRg(), model.getDmax()));
//
//            mdl_cnt++;
//        }

    }


    /*
     * use for mixture models, given 10 models, try all pairs, up to all 10
     * Not to be used for ensemble selection
     *
     */
    void combinatorial_search(IofQData & iofqdata, float cx, float bfactor){

        //iofqdata.makeWorkingSet(11); // meaningless number
        iofqdata.setAllDataToWorkingSet();
        unsigned int totalInWorkingSet = iofqdata.getTotalInWorkingSet();
        auto qvalues = iofqdata.getWorkingSetQvalues();
        const std::vector<Datum> & workingSet = iofqdata.getWorkingSet();

        unsigned lmax = (unsigned int)((iofqdata.getQmax()*iofqdata.getDmax())/SIMDE_MATH_PI) + 1;
        SASTOOLS_UTILS_H::logger("LMAX", std::to_string(lmax));

        std::vector<float> inv_variance(totalInWorkingSet);
        boost::numeric::ublas::matrix<float> i_obs(totalInWorkingSet, 1);

        // rescale data to max at 1
        float max_I = -FLT_MIN;
        for(unsigned int q =0; q < totalInWorkingSet; q++) {
            const Datum *pW = &workingSet[q];
            if (max_I < pW->getI()){
                max_I = pW->getI();
            }
        }

        // rescale the observed datasets
        for(unsigned int q =0; q < totalInWorkingSet; q++) {
            const Datum *pW = &workingSet[q];
            inv_variance[q] = pW->getInvVar()/max_I/max_I;
            i_obs(q, 0) = pW->getI()/max_I;
        }

        boost::numeric::ublas::matrix<float> coeffs(filenames.size(),1);
        boost::numeric::ublas::matrix<float> models(totalInWorkingSet,filenames.size());

        std::vector<unsigned int> selected_indices_of_model_iqs;
        unsigned int mdl_cnt=0;

        for (auto & filename : filenames){

            try{
                AtomisticModel model = AtomisticModel(filename, false, false);

                model.calculatePartialAmplitudes(lmax, totalInWorkingSet, qvalues, false);

                if (model.getPDBModel().getTotalResidues() < 1){
                    throw (filename);
                }

                // create water model
                Waters waters = Waters();
                waters.hydrateAtomisticModel(model); // add waters to PDB
                waters.createSphericalCoordinateOfHydration();
                waters.calculatePartialAmplitudes(lmax, totalInWorkingSet, qvalues);

                std::vector<float> aPs = std::vector<float>(totalInWorkingSet);
                std::vector<float> aPWs = std::vector<float>(totalInWorkingSet);
                std::vector<float> aCs = std::vector<float>(totalInWorkingSet);
                std::vector<float> aXW_cross_term = std::vector<float>(totalInWorkingSet);

                // assemble the partials and
                this->assemblePartials(lmax, model, waters, totalInWorkingSet, aPs, aPWs, aCs, aXW_cross_term);

                // combine and add noise
                for(unsigned int q=0; q<totalInWorkingSet; q++){
                    float q_val = qvalues[q]; // get qvalue to calculate b-factor
                    float exp_term = expf(-(q_val*q_val)*bfactor*inv16Pi2)*cx;
                    // from assemblePartials in Ensemble
                    //        pAP[q] = a_p_norm;
                    //        pAc[q] = a_x_norm; // excluded volume
                    //        pApw[q] = a_p_norm + a_w_norm + 2.0f*(aPWR_term + aPWI_term); // hydratedParticle term
                    //        pAXW_cross_term[q] = -2.0f*(aXWR_term + aXWI_term + aPXR_term + aPXI_term); // invert phase
                    models(q, mdl_cnt) = aPWs[q] + aCs[q]*exp_term*exp_term + exp_term*aXW_cross_term[q];
                }

                params.emplace_back(Param(filename, model.calculateRg(), model.getDmax()));
                selected_indices_of_model_iqs.push_back(mdl_cnt);
                mdl_cnt++;
            } catch (std::string & fl){
                SASTOOLS_UTILS_H::logger("Improper PDB file", "skipping..." + fl);
            }
        }

        if (selected_indices_of_model_iqs.size() > 10){
            throw std::invalid_argument("Too many models for Mixture Fitting > 10");
        }

        SASTOOLS_UTILS_H::logger("Starting", "Combinatorial SEARCH");

        // need to rescale each model to go from 1
        unsigned int total_models = models.size2();
        for(unsigned m=0; m<total_models; m++){
            max_I = -FLT_MIN;
            for(unsigned int q=0; q<totalInWorkingSet; q++){
                if (max_I < models(q, m)){
                    max_I = models(q, m);
                }
            }

            for(unsigned int q=0; q<totalInWorkingSet; q++){
                models(q, m) /= max_I;
            }
        }

        unsigned int maxk = (selected_indices_of_model_iqs.size()) < 7 ? selected_indices_of_model_iqs.size() : selected_indices_of_model_iqs.size() - 1;

        this->combinatorial(
                            inv_variance,
                            selected_indices_of_model_iqs,
                            models,
                            i_obs,
                            totalInWorkingSet,
                            maxk);
    }

    float crossEntropyOptimization(std::vector<float> &i_obs_over_variance,
                                   std::vector<float> &inv_variance,
                                   boost::numeric::ublas::matrix<float> &coeffs,
                                   boost::numeric::ublas::matrix<float> &model_iqs,
                                   std::vector<float> & qvalues,
                                   boost::numeric::ublas::matrix<float> &i_calc_all,
                                   boost::numeric::ublas::matrix<float> &i_obs);


    // ens.combinations(values, 1, 50, 1, 2);
    // prints all pairs => maxk = 2
    // prints all triples => maxk = 3
    void combinations (std::vector<int> vec, int start, int n, int k, int maxk) {

        int i;

        /* k here counts through positions in the maxk-element v.
         * if k > maxk, then the v is complete and we can use it.
         */
        if (k > maxk) {
            /* insert code here to use combinations as you please */
            for (i=1; i<=maxk; i++){
                printf ("%i ", vec[i]);
            }

            printf ("\n");
            return;
        }

        /* for this k'th element of the v, try all start..n
         * elements in that position
         */
        for (i=start; i<=n; i++) {

            vec[k] = i;
            /* recursively generate combinations of integers
             * from i+1..n
             */
            combinations (vec, i+1, n, k+1, maxk);
        }
    }


    void combinatorial (
                        std::vector<float> &inv_variance,
                        std::vector<unsigned int> selected_indices_of_model_iqs,
                        boost::numeric::ublas::matrix<float> &model_iqs,
                        boost::numeric::ublas::matrix<float> &i_obs,
                        unsigned int totalQInWorkingSet,
                        int maxk);


    void increment_neighbor(std::vector<unsigned int> &arr, int idx) {

        for (int i = idx; i < arr.size(); ++i)
            arr.at(i) = arr.at(i - 1) + 1;
    }

    bool _next_combination(int n,
                           int r,
                           std::vector<unsigned int> &curr_comb,
                           int idx) {

        if (idx < 0)
            return false;

        if (curr_comb.at(idx) < n - r + idx + 1) {
            ++curr_comb.at(idx);
            increment_neighbor(curr_comb, idx + 1);
            return true;
        }
        else
            return _next_combination(n, r, curr_comb, idx - 1);
    }

    bool next_combination(int n,
                          int r,
                          std::vector<unsigned int> &curr_comb){

        return _next_combination(n, r, curr_comb, r - 1);
    }

    float calculate_scalar(unsigned int totalInWorkingSet, boost::numeric::ublas::matrix<float> & i_calc,
                           std::vector<float> & i_obs_over_variance,
                           std::vector<float> & inv_variance){

        float c_obs_calc = 0.0f;
        float calc_var = 0.0f;
        float * i_of_q;

        for(unsigned int q =0; q < totalInWorkingSet; q++) {
            i_of_q = &i_calc(q,0);
            c_obs_calc += *i_of_q * i_obs_over_variance[q];
            calc_var += *i_of_q * *i_of_q * inv_variance[q];
        }

        return c_obs_calc/calc_var;
    }
};


#endif //WETSAXS_ENSEMBLE_H
