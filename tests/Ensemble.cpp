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

#include "gtest/gtest.h"
#include "support.hpp"
#include <Ensemble.h>
#include <random>

class EnsembleTests : public ::testing::Test {

public:
    IofQData iofqdata;
    // setup
    // BSA_0p7.dat
    // BSA_high.dat
    EnsembleTests() : ::testing::Test(),
                 iofqdata(IofQData(fixture(BSA_b_sx.dat), false)){}


    void addNoise(unsigned int totalInWorkingSet, float scaleIt, std::vector<float> &i_obs,
                  std::vector<float> & i_obs_over_variance,
                  const std::vector<Datum> & workingSet,
                  boost::numeric::ublas::matrix<float> & y_vector){

        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<float>distribution(0, 1.0);

        for(unsigned int qq=0; qq<totalInWorkingSet; qq++){

            float * iobsq = &i_obs[qq];

            float p=0;
            float snr = 0;
            float noise = 0;
            float val2;

            if (qq < 5){
                unsigned int start_index = (qq < 2) ? 0 : qq - 2;
                for(unsigned int i=0; i< 5; i++){
                    Datum * pW = const_cast<Datum *>(&workingSet[start_index + i]);
                    p += (i_obs[start_index + i] * i_obs[start_index + i]); // calculated
                    snr += (pW->getI()*pW->getI())/(pW->getSigma()*pW->getSigma());
                    val2 = distribution(gen);
                    noise += val2*val2;
                }

            } else if (qq < totalInWorkingSet - 5){
                for(unsigned int i=0; i<5; i++){
                    Datum * pW = const_cast<Datum *>(&workingSet[qq + i - 2]);
                    p += (i_obs[qq + i - 2] * i_obs[qq + i - 2]);
                    snr += (pW->getI()*pW->getI())/(pW->getSigma()*pW->getSigma());
                    val2 = distribution(gen);
                    noise += val2*val2;
                }

            } else {
                unsigned int start_index = (qq >= totalInWorkingSet - 2) ? totalInWorkingSet - 5 : qq - 2;
                for(unsigned int i=0; i<5; i++){
                    Datum * pW = const_cast<Datum *>(&workingSet[start_index + i]);
                    p += (i_obs[start_index + i] * i_obs[start_index + i]);
                    snr += (pW->getI()*pW->getI())/(pW->getSigma()*pW->getSigma());
                    val2 = distribution(gen);
                    noise += val2*val2;
                }
            }

            snr /= 5.0;

            //float p = *iobsq * *iobsq;// *scaleIt *scaleIt;//; power is sum of squared values divided by length
            float np = p/(snr*snr);//10.0f*std::log10(p)/snr; // p/(snr*snr);
            float sigma = sqrtf(np/2);
            float alpha = sqrtf(p/(snr*noise));
            float ori = *iobsq;
            *iobsq = *iobsq + alpha*distribution(gen);
            y_vector(qq, 0) = *iobsq*scaleIt;

            Datum * pW = const_cast<Datum *>(&workingSet[qq]);
            i_obs_over_variance[qq] =  y_vector(qq, 0) * pW->getInvVar();
        }
    }
};



TEST_F(EnsembleTests, ConstructorTest){

    std::vector<std::string> files;

    files.emplace_back(fixture(bsa_md_1_7.pdb));

    boost::filesystem::path full_path(boost::filesystem::current_path());
    std::cout << "Current path is : " << full_path << std::endl;

    Ensemble ens = Ensemble(files, full_path.string());

    ASSERT_EQ(ens.getTotalFiles(), 1);

}

TEST_F(EnsembleTests, ConstructorDirDoesNotExistExceptionTest){

    std::vector<std::string> files;
    files.emplace_back("doesnotExist.txt");

    EXPECT_ANY_THROW(Ensemble(files, "ohno"));
}


TEST_F(EnsembleTests, ConstructorDirExistsExceptionTest){

    std::vector<std::string> files;
    files.emplace_back("doesnotExist.txt");

    std::string dir = boost::filesystem::current_path().string();

    EXPECT_NO_THROW(Ensemble(files, dir));
}

TEST_F(EnsembleTests, ConstructorDirAndFilesExistsExceptionTest){
    std::vector<std::string> files;

    //files.emplace_back(fixture(bsa.pdb));
    files.emplace_back(fixture(bsa_md_1_7.pdb));
    files.emplace_back(fixture(bsa_md_1_997.pdb));
    files.emplace_back(fixture(bsa_md_1_9.pdb));
    files.emplace_back(fixture(bsa_md_1_70.pdb));
    files.emplace_back(fixture(bsa_md_1_100.pdb));
    files.emplace_back(fixture(bsa_md_1_200.pdb));
    files.emplace_back(fixture(bsa_md_1_300.pdb)); // 6
    files.emplace_back(fixture(bsa_md_1_301.pdb));
    files.emplace_back(fixture(bsa_md_1_302.pdb));
    files.emplace_back(fixture(bsa_md_1_400.pdb));
    files.emplace_back(fixture(bsa_md_1_700.pdb)); // 10
    files.emplace_back(fixture(bsa_md_1_701.pdb));
    files.emplace_back(fixture(bsa_md_1_999.pdb));

    Ensemble ens = Ensemble(files, test_dir());

    iofqdata.extractData();
    iofqdata.setDmax(92);

    unsigned lmax = (unsigned int)((iofqdata.getQmax()*iofqdata.getDmax())/SIMDE_MATH_PI) + 1;
    ens.search(lmax, &iofqdata, false, false, 4);

}


TEST_F(EnsembleTests, SingularValueDeterminationTest){

    boost::numeric::ublas::matrix<float> um;
    boost::numeric::ublas::matrix<float> sm;
    boost::numeric::ublas::matrix<float> vm;

    std::vector<std::string> files;
    files.emplace_back(fixture(bsa_md_1_7.pdb));
    files.emplace_back(fixture(bsa_md_1_997.pdb));
    files.emplace_back(fixture(bsa_md_1_9.pdb));
    files.emplace_back(fixture(bsa_md_1_70.pdb));
    files.emplace_back(fixture(bsa_md_1_100.pdb));
    files.emplace_back(fixture(bsa_md_1_200.pdb));
    files.emplace_back(fixture(bsa_md_1_300.pdb)); // 6
    files.emplace_back(fixture(bsa_md_1_301.pdb));
    files.emplace_back(fixture(bsa_md_1_302.pdb));
    files.emplace_back(fixture(bsa_md_1_400.pdb));
    files.emplace_back(fixture(bsa_md_1_700.pdb)); // 10
    files.emplace_back(fixture(bsa_md_1_701.pdb));
    files.emplace_back(fixture(bsa_md_1_999.pdb));

    // m > n for Householder
    std::vector<float> qvalues;
    float minq = 0.001, maxq = 0.31;
    unsigned int divisions = 601;
    float deltaq = (maxq-minq)/(float)divisions;
    for(unsigned int q=0; q<divisions; q++){
        qvalues.emplace_back(minq + q*deltaq);
    }

    unsigned int cnt=0, totalInWorkingSet = qvalues.size();
    unsigned lmax = (unsigned int)((maxq*91)/SIMDE_MATH_PI) + 1;

    boost::numeric::ublas::matrix<float> exam(totalInWorkingSet,files.size());


    Ensemble ens = Ensemble(files, test_dir());

    for (auto & filename : files){

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
            ens.assemblePartials(lmax, model, waters, totalInWorkingSet, aPs, aPWs, aCs, aXW_cross_term);

            for(unsigned int q=0; q<totalInWorkingSet; q++){
                exam(q, cnt) = aPs[q];
            }
            cnt++;
        } catch (std::string & fl){
            SASTOOLS_UTILS_H::logger("Improper PDB file", "skipping..." + fl);
        }
    }

    ens.svd(exam, um, sm, vm);

//    boost::numeric::ublas::matrix<float> astar;
//    ens.pseudo_inverse(lmax, astar, um, sm, vm);

    // based on condition number, the number of independent columns should be Lmax
    unsigned int limit = (sm.size1() < sm.size2()) ? sm.size1() : sm.size2();
    std::cout << " LMAX " << lmax << " row " << sm.size1() << " cols " << sm.size2() << std::endl;
    for(unsigned int i=0; i < limit; i++){
        std::cout << i << "  " << sm(i,i) << " " << sm(0,0)/sm(i,i) << " " << ((i > 0) ? (sm(i-1,i-1) - sm(i,i))/sm(i-1,i-1) : 0) << std::endl;
    }

}


TEST_F(EnsembleTests, CEDecompositionTest){

    iofqdata.extractData();
    iofqdata.setDmax(92);
    //iofqdata.makeWorkingSet(11);
    iofqdata.setAllDataToWorkingSet();
    unsigned int totalInWorkingSet = iofqdata.getTotalInWorkingSet();
    auto qvalues = iofqdata.getWorkingSetQvalues();
    const std::vector<Datum> & workingSet = iofqdata.getWorkingSet();

    std::vector<float> i_obs_over_variance(totalInWorkingSet);
    std::vector<float> inv_variance(totalInWorkingSet);

    for(unsigned int q =0; q < totalInWorkingSet; q++) {
        const Datum *pW = &workingSet[q];
        i_obs_over_variance[q] = pW->getI() * pW->getInvVar();
        inv_variance[q] = pW->getInvVar();
    }


    std::vector<std::string> files;
    files.emplace_back(fixture(1REF_single.pdb));
    files.emplace_back(fixture(6yyt.pdb));
    files.emplace_back(fixture(bsa_md_1_300.pdb)); // 6
    files.emplace_back(fixture(bent_6.pdb)); // 10

// all same protein
//    files.emplace_back(fixture(bsa_md_1_7.pdb));
//    files.emplace_back(fixture(bsa_md_1_997.pdb));
//    files.emplace_back(fixture(bsa_md_1_300.pdb)); // 6
//    files.emplace_back(fixture(bsa_md_1_700.pdb)); // 10

    std::vector<float> ratios = {0.17, 0.23, 0.53, 0.07};

    std::vector<float> i_obs(qvalues.size());
    std::fill(i_obs.begin(), i_obs.end(),0.0f);

    unsigned int cnt=0;
    float maxq = qvalues[totalInWorkingSet-1];

    unsigned lmax = (unsigned int)((maxq*92)/SIMDE_MATH_PI) + 1;

    boost::numeric::ublas::matrix<float> models(totalInWorkingSet,files.size());
    boost::numeric::ublas::matrix<float> y_vector(totalInWorkingSet, 1);

    Ensemble ens = Ensemble(files, test_dir());

    for (auto & filename : files){

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
            ens.assemblePartials(lmax, model, waters, totalInWorkingSet, aPs, aPWs, aCs, aXW_cross_term);

            // combine and add noise
            for(unsigned int q=0; q<totalInWorkingSet; q++){
                models(q, cnt) = aPs[q];
                i_obs[q] += ratios[cnt]*aPs[q];
            }

            cnt++;
        } catch (std::string & fl){
            SASTOOLS_UTILS_H::logger("Improper PDB file", "skipping..." + fl);
        }
    }

    // scale it
    float scaleIt = 0.0;
    for(unsigned int qq=0; qq<7; qq++){
        Datum * pW = const_cast<Datum *>(&workingSet[qq]);
        scaleIt += pW->getI()/i_obs[qq];
    }
    scaleIt /= 7.0;

    addNoise(totalInWorkingSet, scaleIt, i_obs, i_obs_over_variance, workingSet, y_vector);

    boost::numeric::ublas::matrix<float> coeffs(models.size1(),1);
    boost::numeric::ublas::matrix<float> icalcall(totalInWorkingSet,1);

    // crossEntropyOptimization
    ens.crossEntropyOptimization(i_obs_over_variance, inv_variance, coeffs, models, qvalues, icalcall, y_vector);
}

TEST_F(EnsembleTests, SVDLinearDecompositionWithNoiseTest){

    iofqdata.extractData();
    iofqdata.setDmax(92);
    iofqdata.setAllDataToWorkingSet();

    // get real noise from a dataset
    unsigned int totalInWorkingSet = iofqdata.getTotalInWorkingSet();
    auto qvalues = iofqdata.getWorkingSetQvalues();
    const std::vector<Datum> & workingSet = iofqdata.getWorkingSet();

    std::vector<float> i_obs_over_variance(totalInWorkingSet);
    std::vector<float> inv_variance(totalInWorkingSet);

    for(unsigned int q =0; q < totalInWorkingSet; q++) {
        const Datum *pW = &workingSet[q];
        i_obs_over_variance[q] = pW->getI() * pW->getInvVar();
        inv_variance[q] = pW->getInvVar();
    }

    std::vector<std::string> files;
    files.emplace_back(fixture(bsa_md_1_7.pdb));
    files.emplace_back(fixture(bsa_md_1_997.pdb));
    files.emplace_back(fixture(bsa_md_1_300.pdb)); // 6
    files.emplace_back(fixture(bsa_md_1_700.pdb)); // 10

    std::vector<float> ratios = {0.17, 0.23, 0.53, 0.07};
    boost::numeric::ublas::matrix<float> factors(4,1);
    factors(0,0) = 0.17;
    factors(1,0) = 0.23;
    factors(2,0) = 0.53;
    factors(3,0) = 0.07;

    // m > n for Householder
    std::vector<float> i_obs(qvalues.size());
    std::fill(i_obs.begin(), i_obs.end(), 0.0f);

    unsigned int cnt=0;
    float maxq = qvalues[totalInWorkingSet-1];

    unsigned lmax = (unsigned int)((maxq*91)/SIMDE_MATH_PI) + 1;

    boost::numeric::ublas::matrix<float> models(totalInWorkingSet,files.size());
    Ensemble ens = Ensemble(files, test_dir());

    float max_Iq = -FLT_MIN;

    for (auto & filename : files){

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
            ens.assemblePartials(lmax, model, waters, totalInWorkingSet, aPs, aPWs, aCs, aXW_cross_term);

            // combine and add noise
            float * aP;
            for(unsigned int q=0; q<totalInWorkingSet; q++){
                aP = &aPs[q];
                models(q, cnt) = *aP;
                i_obs[q] += ratios[cnt]* *aP;
                if (*aP > max_Iq){
                    max_Iq = *aP;
                }
            }

            cnt++;
        } catch (std::string & fl){
            SASTOOLS_UTILS_H::logger("Improper PDB file", "skipping..." + fl);
        }
    }

    for(unsigned int col =0; col<files.size(); col++){
        float sum = 0;
        for(unsigned int q=0; q<5; q++){
            sum += models(q, col); // max_Iq; // rescale all the models to be between < 1
        }

        sum /= (float)5;
        for(unsigned int q=0; q<totalInWorkingSet; q++){
            models(q, col) /= sum; // scale data to 1
        }
    }

    float sum = 0;
    for(unsigned int q=0; q<5; q++){
        sum += i_obs[q];
    }

    /*
     * scale the data to 1
     */
    float scaleIt = 0.0;
    for(unsigned int qq=0; qq<7; qq++){
        Datum * pW = const_cast<Datum *>(&workingSet[qq]);
        scaleIt += pW->getI()/i_obs[qq];
    }

    scaleIt = 1; // i_obs[0];

    boost::numeric::ublas::matrix<float> y_vector(totalInWorkingSet, 1);
    sum /= (float)5;

    //   for(unsigned int q=0; q<totalInWorkingSet; q++){
//        i_obs[q] /= sum;
//        y_vector(q,0) = i_obs[q];
//    }

    addNoise(totalInWorkingSet, scaleIt, i_obs, i_obs_over_variance, workingSet, y_vector);

    for(unsigned int q=0; q<totalInWorkingSet; q++){
        y_vector(q,0) /= sum;
    }

    for(unsigned int i=0; i<i_obs.size(); i++){
        std::cout << i << " " << qvalues[i] << " " << y_vector(i,0) << std::endl;
    }

    // perform SVD and pseudoinverse
    boost::numeric::ublas::matrix<float> um;
    boost::numeric::ublas::matrix<float> sm;
    boost::numeric::ublas::matrix<float> vm;
    boost::numeric::ublas::matrix<float> astar;

    ens.svd(models, um, sm, vm);
    unsigned int limitSm = (lmax < models.size2()) ? lmax : models.size2();
    ens.pseudo_inverse(limitSm, astar, um, sm, vm);

    boost::numeric::ublas::matrix<float> constants = boost::numeric::ublas::prod(astar, y_vector);

    std::cout << constants << std::endl;

    sum=0.0;
    for (int i=0; i<constants.size1(); i++){
        sum += constants(i,0);
    }

    for (int i=0; i<constants.size1(); i++){
        std::cout << i << " " << constants(i,0)/sum << std::endl;
    }

    //std::cout << constants << std::endl;
    // test seems to fail with noise

    ASSERT_NEAR(constants(0,0), factors(0,0), 0.03);
    ASSERT_NEAR(constants(1,0), factors(1,0), 0.03);
    ASSERT_NEAR(constants(2,0), factors(2,0), 0.03);
    ASSERT_NEAR(constants(3,0), factors(3,0), 0.03);

}

TEST_F(EnsembleTests, SVDLinearDecompositionTest){

    std::vector<std::string> files;
    files.emplace_back(fixture(bsa_md_1_7.pdb));
    files.emplace_back(fixture(bsa_md_1_997.pdb));
    files.emplace_back(fixture(bsa_md_1_300.pdb)); // 6
    files.emplace_back(fixture(bsa_md_1_700.pdb)); // 10

    std::vector<float> ratios = {0.17, 0.23, 0.53, 0.07};
    boost::numeric::ublas::matrix<float> factors(4,1);
    factors(0,0) = 0.17;
    factors(1,0) = 0.23;
    factors(2,0) = 0.53;
    factors(3,0) = 0.07;

    // m > n for Householder
    std::vector<float> qvalues;
    float minq = 0.001, maxq = 0.31;
    unsigned int divisions = 601;
    float deltaq = (maxq-minq)/(float)divisions;
    for(unsigned int q=0; q<divisions; q++){
        qvalues.emplace_back(minq + q*deltaq);
    }

    std::vector<float> i_obs(qvalues.size());
    std::fill(i_obs.begin(), i_obs.end(),0.0f);

    unsigned int cnt=0, totalInWorkingSet = qvalues.size();
    maxq = qvalues[totalInWorkingSet-1];

    unsigned lmax = (unsigned int)((maxq*91)/SIMDE_MATH_PI) + 1;

    boost::numeric::ublas::matrix<float> models(totalInWorkingSet,files.size());

    Ensemble ens = Ensemble(files, test_dir());

    for (auto & filename : files){

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
            ens.assemblePartials(lmax, model, waters, totalInWorkingSet, aPs, aPWs, aCs, aXW_cross_term);

            // combine and add noise
            for(unsigned int q=0; q<totalInWorkingSet; q++){
                models(q, cnt) = aPs[q];
                i_obs[q] += ratios[cnt]*aPs[q];
            }

            cnt++;
        } catch (std::string & fl){
            SASTOOLS_UTILS_H::logger("Improper PDB file", "skipping..." + fl);
        }
    }


    for(unsigned m=0; m<models.size2(); m++){
        float max_I = -FLT_MIN;
        for(unsigned int q=0; q<totalInWorkingSet; q++){
            if (max_I < models(q, m)){
                max_I = models(q, m);
            }
        }

        for(unsigned int q=0; q<totalInWorkingSet; q++){
            models(q, m) /= max_I;
        }
    }



    boost::numeric::ublas::matrix<float>i_obs_syn = boost::numeric::ublas::prod(models, factors);

    // perform SVD and pseudoinverse
    boost::numeric::ublas::matrix<float> um;
    boost::numeric::ublas::matrix<float> sm;
    boost::numeric::ublas::matrix<float> vm;
    boost::numeric::ublas::matrix<float> astar;

    //The models matrix will be modified in the SVD by Househoulder transformation
//    boost::numeric::ublas::matrix<float> ori(models);

    ens.svd(models, um, sm, vm);
    unsigned int limitSm = (lmax < models.size2()) ? lmax : models.size2();
    ens.pseudo_inverse(limitSm, astar, um, sm, vm);

    // multiple obs vector by pseudo inverse
    boost::numeric::ublas::matrix<float> constants = boost::numeric::ublas::prod(astar, i_obs_syn);

//    boost::numeric::ublas::matrix<float> calc = boost::numeric::ublas::prod(ori, constants);
//    float rsum = 0;
//
//    for(unsigned int i=0; i<models.size1(); i++){
//        float residual = (calc(i,0) - i_obs_syn(i,0));
//        std::cout << i << " " << residual << " " << calc(i,0) << " " << i_obs_syn(i,0) << std::endl;
//        rsum += residual*residual;
//    }


    std::cout << constants << std::endl;
    ASSERT_NEAR(constants(0,0), factors(0,0), 0.001);
    ASSERT_NEAR(constants(1,0), factors(1,0), 0.001);
    ASSERT_NEAR(constants(2,0), factors(2,0), 0.001);
    ASSERT_NEAR(constants(3,0), factors(3,0), 0.001);

}

TEST_F(EnsembleTests, NNLSLinearDecompositionTest){

    std::vector<std::string> files;
    files.emplace_back(fixture(bsa_md_1_7.pdb));
    files.emplace_back(fixture(bsa_md_1_997.pdb));
    files.emplace_back(fixture(bsa_md_1_300.pdb)); // 6
    files.emplace_back(fixture(bsa_md_1_700.pdb)); // 10

    std::vector<float> ratios = {0.17, 0.23, 0.53, 0.07};
    boost::numeric::ublas::matrix<float> factors(4,1);
    factors(0,0) = ratios[0];
    factors(1,0) = ratios[1];
    factors(2,0) = ratios[2];
    factors(3,0) = ratios[3];

    // m > n for Householder
    std::vector<float> qvalues;
    float minq = 0.001, maxq = 0.31;
    unsigned int divisions = 601;
    float deltaq = (maxq-minq)/(float)divisions;
    for(unsigned int q=0; q<divisions; q++){
        qvalues.emplace_back(minq + q*deltaq);
    }

    std::vector<float> i_obs(qvalues.size());
    std::fill(i_obs.begin(), i_obs.end(),0.0f);

    unsigned int cnt=0, totalInWorkingSet = qvalues.size();
    maxq = qvalues[totalInWorkingSet-1];

    unsigned lmax = (unsigned int)((maxq*91)/SIMDE_MATH_PI) + 1;

    nsNNLS::denseMatrix modelsNNLS(totalInWorkingSet, files.size());
    nsNNLS::vector y_vector(totalInWorkingSet);

    Ensemble ens = Ensemble(files, test_dir());

    for (auto & filename : files){

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
            ens.assemblePartials(lmax, model, waters, totalInWorkingSet, aPs, aPWs, aCs, aXW_cross_term);

            // combine and add noise
            for(unsigned int q=0; q<totalInWorkingSet; q++){
                modelsNNLS.set(q, cnt, aPs[q]);
                i_obs[q] += ratios[cnt]*aPs[q];
            }

            cnt++;
        } catch (std::string & fl){
            SASTOOLS_UTILS_H::logger("Improper PDB file", "skipping..." + fl);
        }
    }

    /*
     * scale the data to 1
     */
//    std::vector<float> scalars(modelsNNLS.ncols(), 0.0);
//    for(unsigned int m=0; m<modelsNNLS.ncols(); m++){
//        float maxv = 0;
//        for(unsigned int q=0; q<totalInWorkingSet; q++){
//            float fval = modelsNNLS.get(q,m);
//            if (maxv < fval){
//                maxv = fval;
//            }
//        }
//        scalars[m] = maxv;
//        for(unsigned int q=0; q<totalInWorkingSet; q++){
//            float fval = modelsNNLS.get(q,m)/maxv;
//            modelsNNLS.set(q,m, fval);
//        }
//    }
//
//    float maxv = -FLT_MAX;
//    for(unsigned int q=0; q<totalInWorkingSet; q++){
//        if (i_obs[q] > maxv){
//            maxv = i_obs[q];
//        }
//    }


    for(unsigned int i=0; i<totalInWorkingSet; i++){
        y_vector.set(i, 1000*i_obs[i]);
    }

    nsNNLS::nnls solver = nsNNLS::nnls( &modelsNNLS, &y_vector, 20000);
    int flag = solver.optimize();
    auto vecb = solver.getSolution();


    float bsum = 0;
    for(unsigned int i=0; i<modelsNNLS.ncols(); i++){
        bsum += vecb->get(i);
    }

    for(unsigned int i=0; i<modelsNNLS.ncols(); i++){
        std::cout << " b " << i << " " << vecb->get(i)/bsum << std::endl;
    }
}


TEST_F(EnsembleTests, NNLSLinearDecompositionWithNoiseTest){

    iofqdata.extractData();
    iofqdata.setDmax(92);
    iofqdata.setAllDataToWorkingSet();

    unsigned int totalInWorkingSet = iofqdata.getTotalInWorkingSet();
    auto qvalues = iofqdata.getWorkingSetQvalues();
    const std::vector<Datum> & workingSet = iofqdata.getWorkingSet();

    std::vector<float> i_obs_over_variance(totalInWorkingSet);
    std::vector<float> inv_variance(totalInWorkingSet);

    for(unsigned int q =0; q < totalInWorkingSet; q++) {
        const Datum *pW = &workingSet[q];
        i_obs_over_variance[q] = pW->getI() * pW->getInvVar();
        inv_variance[q] = pW->getInvVar();
    }

    std::vector<std::string> files;
    files.emplace_back(fixture(1REF_single.pdb));
    files.emplace_back(fixture(6yyt.pdb));
    files.emplace_back(fixture(bsa_md_1_300.pdb)); // 6
    files.emplace_back(fixture(bent_6.pdb)); // 10
//    files.emplace_back(fixture(bsa_md_1_7.pdb));
//    files.emplace_back(fixture(bsa_md_1_997.pdb));
//    files.emplace_back(fixture(bsa_md_1_300.pdb)); // 6
//    files.emplace_back(fixture(bsa_md_1_700.pdb)); // 10

    std::vector<float> ratios = {0.17, 0.23, 0.53, 0.07};
    boost::numeric::ublas::matrix<float> factors(files.size(),1);
    factors(0,0) = 0.17;
    factors(1,0) = 0.23;
    factors(2,0) = 0.53;
    factors(3,0) = 0.07;

    // m > n for Householder
    unsigned int cnt=0;
    float maxq = qvalues[totalInWorkingSet-1];

    unsigned lmax = (unsigned int)((maxq*91)/SIMDE_MATH_PI) + 1;

    nsNNLS::denseMatrix modelsNNLS(totalInWorkingSet, files.size());
    Ensemble ens = Ensemble(files, test_dir());

    for (auto & filename : files){

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
            ens.assemblePartials(lmax, model, waters, totalInWorkingSet, aPs, aPWs, aCs, aXW_cross_term);

            // combine and add noise
            for(unsigned int q=0; q<totalInWorkingSet; q++){
                modelsNNLS.set(q, cnt, aPs[q]);
            }

            cnt++;
        } catch (std::string & fl){
            SASTOOLS_UTILS_H::logger("Improper PDB file", "skipping..." + fl);
        }
    }

    std::vector<float> i_obs(qvalues.size());
    std::fill(i_obs.begin(), i_obs.end(),0.0f);
    /*
     * scale the data to 1
     */
    std::vector<float> scalars(modelsNNLS.ncols(), 0.0);
    for(unsigned int m=0; m<modelsNNLS.ncols(); m++){
        float maxv = 0;
        for(unsigned int q=0; q<totalInWorkingSet; q++){
            float fval = modelsNNLS.get(q,m);
            if (maxv < fval){
                maxv = fval;
            }
        }
        scalars[m] = maxv;
        float fraction = ratios[m];
        for(unsigned int q=0; q<totalInWorkingSet; q++){
            float fval = modelsNNLS.get(q,m)/maxv;
            modelsNNLS.set(q,m, fval);
            i_obs[q] += fraction*fval;
        }
    }

    boost::numeric::ublas::matrix<float> y_vector(totalInWorkingSet, 1);
    addNoise(totalInWorkingSet, 1.0, i_obs, i_obs_over_variance, workingSet, y_vector);

    nsNNLS::vector y_vectorNNLS(totalInWorkingSet);
    for(int i=0; i<totalInWorkingSet; i++){
        y_vectorNNLS.set(i, y_vector(i,0));
    }

    nsNNLS::nnls solver = nsNNLS::nnls( &modelsNNLS, &y_vectorNNLS, 200000);
    int flag = solver.optimize();
    auto vecb = solver.getSolution();

    float bsum = 0;
    for(unsigned int i=0; i<modelsNNLS.ncols(); i++){
        bsum += vecb->get(i);
    }

    for(unsigned int i=0; i<modelsNNLS.ncols(); i++){
        std::cout << " b " << i << " " << vecb->get(i)/bsum << std::endl;
    }

    nsNNLS::vector ax = nsNNLS::vector(modelsNNLS.nrows());
    modelsNNLS.dot(false, vecb, &ax);

//    for(unsigned int i=0;i<totalInWorkingSet; i++){
//        std::cout << qvalues[i] << " " << y_vector(i,0) << " " <<  ax.get(i) << std::endl;
//    }

    ASSERT_NEAR(vecb->get(0)/bsum, factors(0,0), 0.05);
    ASSERT_NEAR(vecb->get(1)/bsum, factors(1,0), 0.05);
    ASSERT_NEAR(vecb->get(2)/bsum, factors(2,0), 0.05);
    ASSERT_NEAR(vecb->get(3)/bsum, factors(3,0), 0.03);
}

TEST_F(EnsembleTests, SVDTest){

    boost::numeric::ublas::matrix<float> um;
    boost::numeric::ublas::matrix<float> sm;
    boost::numeric::ublas::matrix<float> vm;

    // m > n for Householder
//    boost::numeric::ublas::matrix<float> exam(4,2);
//    exam(0,0)=2.0f;
//    exam(0,1)=4.0f;
//    exam(1,0)=1.0f;
//    exam(1,1)=3.0f;
//    exam(2,0)=0.0f;
//    exam(2,1)=0.0f;
//    exam(3,0)=0.0f;
//    exam(3,1)=0.0f;
    // m > n for Householder
    // Sigmas: [2](3.23607,1.23607)
    // 3 -1
    // 1  3
    // 1  1
    boost::numeric::ublas::matrix<float> exam(3,2);
    exam(0,0)=3.0f;
    exam(1,0)=1.0f;
    exam(2,0)=1.0f;
    exam(0,1)=-1.0f;
    exam(1,1)=3.0f;
    exam(2,1)=1.0f;

    std::vector<std::string> files;
    files.emplace_back(fixture(bsa_md_1_302.pdb));
    files.emplace_back(fixture(bsa_md_1_400.pdb));
    files.emplace_back(fixture(bsa_md_1_700.pdb)); // 10
    files.emplace_back(fixture(bsa_md_1_701.pdb));
    files.emplace_back(fixture(bsa_md_1_999.pdb));

    Ensemble ens = Ensemble(files, test_dir());

    ens.svd(exam, um, sm, vm);

    std::cout << um << std::endl; // U-matrix
    std::cout << sm << std::endl; // singular values
    std::cout << vm << std::endl; // V-matrix (needs to be transposed)
    //boost::numeric::ublas::trans(vm);
    ASSERT_NEAR(3.16228, sm(0,0), 0.00001);
    ASSERT_NEAR(3.4641, sm(1,1), 0.0001);
}

TEST_F(EnsembleTests, CombinationTestsExceptionTest){

    float bfactor = 7.9;
    float cx = 0.817;
    const float inv16Pi2 = 1.0/(16.0*M_PI_2);

    std::vector<std::string> files;

    iofqdata.extractData();
    iofqdata.setDmax(92);
    //iofqdata.setAllDataToWorkingSet();
    iofqdata.makeWorkingSet(11);

    unsigned int totalInWorkingSet = iofqdata.getTotalInWorkingSet();
    auto qvalues = iofqdata.getWorkingSetQvalues();
    const std::vector<Datum> & workingSet = iofqdata.getWorkingSet();

    unsigned lmax = (unsigned int)((iofqdata.getQmax()*iofqdata.getDmax())/SIMDE_MATH_PI) + 1;

    files.emplace_back(fixture(bsa.pdb));           //0
    files.emplace_back(fixture(bsa_md_1_7.pdb));    //1
    files.emplace_back(fixture(bsa_md_1_997.pdb));  //2
    files.emplace_back(fixture(bsa_md_1_9.pdb));    //3
    files.emplace_back(fixture(bsa_md_1_70.pdb));   //4
    files.emplace_back(fixture(bsa_md_1_400.pdb));  //5
    files.emplace_back(fixture(bsa_md_1_700.pdb));  //6
    files.emplace_back(fixture(bsa_md_1_701.pdb));  //7

    Ensemble ens = Ensemble(files, test_dir());

    boost::numeric::ublas::matrix<float> models(totalInWorkingSet,files.size());

    unsigned int mdl_cnt=0;

    for (auto & filename : files){

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
            ens.assemblePartials(lmax, model, waters, totalInWorkingSet, aPs, aPWs, aCs, aXW_cross_term);

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

            mdl_cnt++;
        } catch (std::string & fl){
            SASTOOLS_UTILS_H::logger("Improper PDB file", "skipping..." + fl);
        }
    }

//    int n = 5, r = 2;
//    unsigned int maxk = 3;
//
//    for (unsigned int k=2; k<=maxk; k++) {
//
//        std::vector<unsigned int> numbers(k);
//        // initialize values in numbers
//        for (unsigned int i = 0; i < k; ++i)
//            numbers.at(i) = i;
//
//        do {
//            for (auto x: numbers)
//                std::cout << x << " ";
//            std::cout << "\n";
//
//        } while (ens.next_combination(n, (int)k, numbers));
//    }

    std::vector<float> inv_variance(totalInWorkingSet);
    boost::numeric::ublas::matrix<float> i_obs(totalInWorkingSet,1);

    for(unsigned int q =0; q < totalInWorkingSet; q++) {
        const Datum *pW = &workingSet[q];
        inv_variance[q] = pW->getInvVar();
        i_obs(q,0) = pW->getI();
    }

    std::vector<unsigned int> selected_indices_of_model_iqs = {0,1,2,3,4,5,6,7};

    ens.combinatorial(
                      inv_variance,
                      selected_indices_of_model_iqs,
                      models,
                      i_obs,
                      totalInWorkingSet,
                      4);
}