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
        double chi, dw;
        boost::numeric::ublas::matrix<float> model_weights;
        Trial() = default;

        Trial(double val, boost::numeric::ublas::matrix<float> & vec) : value(val) {
            model_weights = std::move(vec);
        }

        bool operator<(const Trial & a) const{
            return value < a.value;
        }
    };

    std::vector<std::string> filenames;
    std::string directory;
    const float inv16Pi2 = 1.0/(16.0*M_PI_2);
    const float EPS = 0.00001;
    const int ITER_MAX = 50;
    float standardizedMin;
    float standardizedScale;

    const float updateAlpha=0.67f;
    float frac = 0.03171;
    unsigned int totalTrials = 21;
    unsigned int trial_size = 100;
    unsigned int totalPerRound = 3331;

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


    void linearCombination(unsigned int lmax, IofQData *iofqfile,
                           unsigned int total_cx,
                           float lowcx,
                           float highcx,
                           float upperB,
                           bool forceRNA=false, bool keepWaters=false,
                           int oversampling=3);

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
    void identity(boost::numeric::ublas::matrix < float >&m)
    {
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

    void svd(boost::numeric::ublas::matrix<float> &A, boost::numeric::ublas::matrix<float> &QQL,
             boost::numeric::ublas::matrix<float> &QQW, boost::numeric::ublas::matrix<float> &QQR);


    void householder(boost::numeric::ublas::matrix<float> &A, boost::numeric::ublas::matrix<float> &QQ,
                     unsigned int row_start, unsigned int col_start,
                     bool column);

    void svd_qr_shift(boost::numeric::ublas::matrix<float> &u, boost::numeric::ublas::matrix<float> &v,
                      boost::numeric::ublas::vector<float> &q, boost::numeric::ublas::vector<float> &e);


    void pseudo_inverse(unsigned int lmax,
                        boost::numeric::ublas::matrix<float> &A_star, boost::numeric::ublas::matrix<float> &um,
                        boost::numeric::ublas::matrix<float> &sm, boost::numeric::ublas::matrix<float> &vm);

//    void nnls(boost::numeric::ublas::matrix<float> &y_m, boost::numeric::ublas::matrix<float> &x_m,
//              boost::numeric::ublas::matrix<float> &A_m);

    unsigned int get_max_index(boost::numeric::ublas::matrix<float> &vector);

    void standardizeData(nsNNLS::vector *vec, nsNNLS::vector *out_vec);

    float getStandardizedScale(){ return standardizedScale;}
    float getStandardizedMin(){ return standardizedMin;}

    float crossEntropyOptimization(std::vector<float> &i_obs_over_variance,
                                   std::vector<float> &inv_variance,
                                   boost::numeric::ublas::matrix<float> &coeffs,
                                   boost::numeric::ublas::matrix<float> &model_iqs,
                                   boost::numeric::ublas::matrix<float> &i_obs);
};


#endif //WETSAXS_ENSEMBLE_H
