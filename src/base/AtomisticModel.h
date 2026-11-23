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

#ifndef WETSAXS_ATOMISTICMODEL_H
#define WETSAXS_ATOMISTICMODEL_H

#include <vector>
#include <string>
#include <random>
#include <sastools/vector3.h>
#include <sastools/PDBModel.h>
#include <sastools/Bead.h>
#include <set>
#include "SphericalHarmonics.h"
#include "SphericalBessels.h"
#include "functions_WETSAXS.h"
#include "Waters.h"

class Waters;

class AtomisticModel {
    struct HashPair {
        template <typename T>
        std::size_t operator()(const std::pair<T, T>& p) const noexcept {
            static_assert(std::is_unsigned_v<T>, "Use unsigned types for best results");

            // SplitMix64 / better mixing for 64-bit
            auto hash1 = splitmix64(p.first);
            auto hash2 = splitmix64(p.second);

            // Combine with good avalanche
            hash1 ^= hash2 + 0x9e3779b97f4a7c15ULL + (hash1 << 6) + (hash1 >> 2);
            return hash1;
        }

    private:
        static uint64_t splitmix64(uint64_t x) noexcept {
            x += 0x9e3779b97f4a7c15ULL;
            x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
            x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
            return x ^ (x >> 31);
        }
    };


    std::string filename;
    bool forceRNA = false;
    bool keepWaters = false;

    PDBModel pdbModel;

    unsigned int totalAtoms;  // maximum number of waters per residue, should match declaration in waters
    unsigned int totalDistances;
    unsigned int lmax;

    float matchPoint = 1.0;
    float forwardScatteringParticle, forwardScatteringExVolume;

    // setup for dynamically allocated arrays
    std::vector<float> rvalues;
    std::vector<float> thetas; // could be aligned memory
    std::vector<float> phis;
    std::vector<float> bead_rvalues;
    std::vector<float> bead_thetas; // could be aligned memory
    std::vector<float> bead_phis;

    std::vector<float> ex_volume_weights;
    std::vector<float> almRealSum; // partial sums
    std::vector<float> almImagSum;
    std::vector<float> almRealExcludedVolSum; // partial sums
    std::vector<float> almImagExcludedVolSum;

    std::vector<float> binnedPartials;

    std::set <int> atomList;

    std::vector<float> fractionOccAtoms;
    std::vector<Bead> beads;

    SphericalHarmonics she, she_x;
    SphericalBessels sbj, sbj_x;

    float rg = 0;

public:
    AtomisticModel()=default;

    AtomisticModel(std::string filename, bool forceRNA, bool keepWaters);

    // copy constructor - prevent copying
    AtomisticModel(const AtomisticModel & model)= delete;
    // copy assignment - prevent copying
    AtomisticModel & operator=(const AtomisticModel & model)= delete;

    // move assignment operator
    AtomisticModel & operator=(AtomisticModel && model) noexcept {

        if (&model == this){
            return *this;
        }

        filename = std::move(model.filename);
        forceRNA = model.forceRNA;
        keepWaters = model.keepWaters;

        pdbModel = std::move(model.pdbModel);
        totalAtoms = model.totalAtoms;
        totalDistances = model.totalDistances;

        lmax = model.lmax;

        matchPoint = model.matchPoint;
        forwardScatteringParticle = model.forwardScatteringParticle;
        forwardScatteringExVolume = model.forwardScatteringExVolume;

        rvalues = std::move(model.rvalues);

        thetas = std::move(model.thetas);
        phis = std::move(model.phis);
        bead_rvalues = std::move(model.bead_rvalues);
        bead_thetas = std::move(model.bead_thetas);
        bead_phis = std::move(model.bead_phis);

        ex_volume_weights = std::move(model.ex_volume_weights);
        almRealSum = std::move(model.almRealSum);
        almImagSum = std::move(model.almImagSum);

        almRealExcludedVolSum = std::move(model.almRealExcludedVolSum); // partial sums
        almImagExcludedVolSum = std::move(model.almImagExcludedVolSum);

        binnedPartials = std::move(model.binnedPartials);
        atomList = std::move(model.atomList);

        fractionOccAtoms = std::move(model.fractionOccAtoms);
        beads = std::move(model.beads);

        she = std::move(model.she);
        she_x = std::move(model.she_x);
        sbj = std::move(model.sbj);
        sbj_x = std::move(model.sbj_x);

        rg = std::move(model.rg);

        return *this;
    }

    AtomisticModel (AtomisticModel && model) noexcept {
        *this = std::move(model);
    }

    ~AtomisticModel()=default;

    /*
     * total atoms in calculation should equal size of rvalues, thetas and phis
     */
    unsigned int getTotalAtoms(){ return totalAtoms;}

    /**
     * total number of unique atoms by atom type or atomic number
     * @return
     */
    unsigned int getTotalUniqueAtoms(){ return atomList.size();}
    std::set <int> & getAtomList(){ return atomList; }

    const float * getThetas(){ return thetas.data();}
    const float * getPhis(){ return phis.data();}
    const float * getRValues() { return rvalues.data();}

    const float * getModelPartialsReal() { return almRealSum.data(); } // partial sums
    const float * getModelPartialsImag() { return almImagSum.data();}

    const float * getExcludedVolumePartialsReal() { return almRealExcludedVolSum.data();} // partial sums
    const float * getExcludedVolumePartialsImag(){ return almImagExcludedVolSum.data();}

    const float * getBinnedPartials() { return binnedPartials.data();}

    std::vector<float> * getRValuesVector() { return &rvalues;}

    void createSphericalHarmonics(unsigned int lmax);

    void createSphericalBessels(unsigned int totalqvalues, std::vector<float> &qvalues);

    const SphericalBessels &getSbj() const;

    PDBModel & getPDBModel() { return pdbModel; }
    //const PDBModel * getPdbModel() const { return &pdbModel; }

    float getDmax(){ return pdbModel.getDmax(); }

    void calculatePartialAmplitudes(unsigned int lmax,
                                    unsigned int totalqvalues,
                                    std::vector<float> &qvalues,
                                    bool recalculate = false);

    void writeInVacuoModel(unsigned int lmax, unsigned int totalqvalues, std::vector<float> &qvalues);

    void calculatePartialAmplitudesNoExcludedVolume(unsigned int lmax,
                                               unsigned int totalqvalues,
                                               std::vector<float> &qvalues);

    void addWatersToAtomicModel(Waters &waterModel);

    void setMatchPointScale(int * pAtomicNumbers);
    void setMatchPointScale(unsigned int totalWatersInHydration);

    float getMatchPoint(){ return matchPoint;}

    void createHCPGridModel(float bin_width);

    //float addUp(unsigned int total_width, float *vector_to_sum);

    void populateSphBesselsArray(int l, std::vector<float> &qr, float *bessel_alq);

    // void setFractionalVolume(float volume){
    //
    //     fractionOccAtoms.resize(totalAtoms);
    //     auto pAtomicVolumes = pdbModel.getAtomicVolumeVec();
    //
    //     float inv = 1.0f/volume;
    //
    //     for(unsigned int i=0; i<totalAtoms; i++){
    //         fractionOccAtoms[i] =  pAtomicVolumes[i]*inv; // normalized to volume of water, looking at fractions of water volume
    //     }
    // }

    /*
     * Convert the Gerstein-Voronoi atomic radii into Gaussian radii
     *
     */
    //    void setGaussianRadii(){
    //        auto pAtomicVolumes = pdbModel.getAtomicVolumeVec();
    //
    //        gaussian_radii.resize(totalAtoms);
    //        float * pvol;
    //        float invpi3halves = 1.0f/std::sqrtf(M_PI*M_PI*M_PI);
    //
    //        for(unsigned int i=0; i<totalAtoms; i++){
    //            pvol = &pAtomicVolumes[i];
    //            gaussian_radii[i] =  std::cbrt(*pvol * invpi3halves); // normalized to volume of water, looking at fractions of water volume
    //        }
    //    }

    /*
     * Gaussian volumes where v = r^3 * pi^3/2
     */
    void setFractionalVolumeAtQ(float q_val){

        fractionOccAtoms.resize(totalAtoms);

        float * gaussian_radii = pdbModel.getAtomicGaussianRadii();

        float factor, diff, inv = 1.0f/(4*M_PI);
        float pw_vol = 29.9;
        float qval2 = q_val*q_val;

        for(unsigned int i=0; i<totalAtoms; i++){

            diff = gaussian_radii[i]/1.75f;
            factor = pw_vol*(1 - diff*diff*diff);

            fractionOccAtoms[i] = 0.334f*factor * std::expf(-std::cbrt(factor*factor)*inv*qval2); // normalized to volume of water, looking at fractions of water volume
        }
    }


    const float * getFractionalVolume() {
        // if (fractionOccAtoms.empty()){
        //     setFractionalVolume(29.9); // 29.9 volume for H2O
        // }
        return fractionOccAtoms.data();
    }


    SphericalHarmonics & getSHE(){ return she;}

    SphericalBessels & getSBJ() { return sbj;}

    float calculateRg(){
        const int * pAtomicNumbers = pdbModel.getAtomicNumberVec();
        rg = 0;
        for(unsigned int n=0; n < totalAtoms; n++){

            float rval = rvalues[n];
            rg += rval*rval; // unique atomic numbers
        }
        rg *= 1/(float)totalAtoms;
        rg = std::sqrtf(rg);
        return rg;
    }

    float getRg(){ return (rg ==0) ? this->calculateRg(): rg ;}

    float calculateStructureFactorWaterAtQ(float qval) {

        float delta_r = 0.02;
        float sum = 0;

        /*
         * f_0(2.27) = 0
         * f_1(2.29) + ... +
         * numerical integration using trapezoid rule
         */
        float rval = 0.0f; // there are no distances less than 2.29 (water to water)
        float endR = 18.2;

        float sincf;

        rval += delta_r;
        unsigned int ind = 1;

        while (rval < 2.3) {
            float qr = qval * rval;
            if (std::abs(qr) < 1e-10) {
                sincf = 1.0f;                 // lim x→0 sin(x)/x = 1
            } else {
                sincf = sinc(qr);
            }

            sum += ((ind % 2 == 0) ? 2.0f : 4.0f)*rval*rval*(-1.0f)*sincf;
            //std::cout << rval << " " << 0 << std::endl;
            rval += delta_r;
            ind++;
        }

        while (rval < endR) {
            float qr = qval * rval;
            if (std::abs(qr) < 1e-10) {
                sincf = 1.0f;                 // lim x→0 sin(x)/x = 1
            } else {
                sincf = sinc(qr);
            }

            // std::cout << rval << " " << rval*rval*(g_OO_SQ_matched(rval)-1.0f) << std::endl;
            // std::cout << rval << " " << g_OO_SQ_matched(rval) << std::endl;

            //std::cout << rval << " " << rval*rval*( calculateWaterWeight(rval)-1.0f ) << std::endl;
            sum += ((ind % 2 == 0) ? 2.0f : 4.0f)*rval*rval*( g_OO_SQ_matched(rval) - 1.0f ) * sincf;
            //sum += ((ind % 2 == 0) ? 2.0f : 4.0f)*rval*rval*( calculateWaterWeight(rval)-1.0f ) * sincf;
            rval += delta_r;
            ind++;
        }


        // for (unsigned int i=0; i<totalN; i++) {// integrate
        //     float qr = qval*rval;
        //     sum += calculateWaterWeight(rval)*sinf(qr)/qr;
        //    //     std::cout << rval << " " << calculateWaterWeight(rval) << std::endl;
        //
        //     rval += delta_r;
        // }

        sum *= delta_r/3.0f;
        return 1.0f + 4.0f*M_PI*0.03342*sum;
    }

    /*
     * returns P(r) distribution for water
     * r_min is 2.3
     */
    float calculateWaterWeight(const float dis) {

        float a_1 = 2.4 ;
        float r_1 = 2.8 ;
        float s_1 = 0.5/(0.28*0.28);
        float a_2 = 0.5 ;
        float r_2 = 4.5 ;
        float s_2 = 0.5/(0.55*0.55);
        float a_3 = 0.2 ;
        float r_3 = 7.1 ;
        float s_3 = 0.5/(0.9*0.9) ;
        // unsigned int num;
        // float minLimit = 2.49*2.49;
        // float maxLimit = 10.1*10.1;
        // P(r) distribution of water => r^2*( g(r)-1 )
        // maximum value at P(2.8) is 1, normalized to peak value
        //float value = ((1 + a_1*expf(-std::powf(dis-r_1,2)*s_1) + a_2*expf(-std::powf(dis-r_2,2)*s_2) + a_3*expf(-std::powf(dis-r_3,2)*s_3))*(1-expf(-8*std::powf(dis-2.3,2))) - 1)
        float g_r = 1 +
         (a_1*expf(-std::powf(dis-r_1,2)*s_1) +
         a_2*expf(-std::powf(dis-r_2,2)*s_2) +
         a_3*expf(-std::powf(dis-r_3,2)*s_3))*
         (1.0f-expf(-8*std::powf( dis - 2.3, 2)));
        return std::max(0.0f, g_r);
    }

    // g_OO(r) plus the long-range correction fitted to match real S(Q) data.
    float g_OO_SQ_matched(float rval) {

        // SHOULD return 0.0f for values less than cutoff
        float diff = rval-7.62330643;
        float g = g_oo(rval) +
            -0.01531494*expf(diff*diff*-2.0926947568f);
        return g;
    }

    float g_oo(const float rval) {

        constexpr double r_min_exclusion = 2.5;
        if (rval < r_min_exclusion) {
            return 0.0;
        }

        float a_1 = 1.6022288043;
        float r_1 = 2.7992424956 ;
        float s_1 = -23.8463006847317;
        float a_2 = -0.5044352646 ;
        float r_2 = 2.512087913132 ;
        float s_2 = -82.5754762463057;
        float a_3 = 2.835902295 ;
        float r_3 = 4.5106758451 ;
        float s_3 = -1.16992522513;
        float a_4 = -2.713383622 ;
        float r_4 = 4.506014238 ;
        float s_4 = -0.999784342623;
        float a_5 = 0.07666115 ;
        float r_5 = 6.689911741 ;
        float s_5 = -2.979445942709;
        float a_6 = -0.881242444 ;
        float r_6 = 2.287563421 ;
        float s_6 = -17.1019985978;

        // unsigned int num;
        // float minLimit = 2.49*2.49;
        // float maxLimit = 10.1*10.1;
        // P(r) distribution of water => r^2*( g(r)-1 )
        // maximum value at P(2.8) is 1, normalized to peak value
        //float value = ((1 + a_1*expf(-std::powf(dis-r_1,2)*s_1) + a_2*expf(-std::powf(dis-r_2,2)*s_2) + a_3*expf(-std::powf(dis-r_3,2)*s_3))*(1-expf(-8*std::powf(dis-2.3,2))) - 1)
        float g_r = 1 +
         a_1*expf(std::powf(rval-r_1,2)*s_1) +
         a_2*expf(std::powf(rval-r_2,2)*s_2) +
         a_3*expf(std::powf(rval-r_3,2)*s_3) +
         a_4*expf(std::powf(rval-r_4,2)*s_4) +
         a_5*expf(std::powf(rval-r_5,2)*s_5) +
         a_6*expf(std::powf(rval-r_6,2)*s_6);
        return g_r;
        //return std::max(0.0f, g_r);
    }

    float sinc(float x) {
        if (std::fabs(x) < 1e-8) {
            return 1.0f - x * x / 6.0f;
        }
        return std::sinf(x) / x;
    }


    float excluded_volume_scattering(float probe_radius=2.1, float q_value = 0.01);

    template <class F>
    void forNear(const vector3 & p,
        int nx,
        int ny,
        int nz,
        const std::vector<std::vector<int>> & cells,
        F&& f) const {

        int cx = std::clamp(static_cast<int>(std::floor(p.x)),0,nx-1);
        int cy = std::clamp(static_cast<int>(std::floor(p.y)),0,ny-1);
        int cz = std::clamp(static_cast<int>(std::floor(p.z)),0,nz-1);

        for (int dx=-1;dx<=1;++dx) for (int dy=-1;dy<=1;++dy) for (int dz=-1;dz<=1;++dz) {
            int nx_=cx+dx, ny_=cy+dy, nz_=cz+dz;
            if (nx_<0||ny_<0||nz_<0||nx_>=nx||ny_>=ny||nz_>=nz) continue;
            for (int idx : cells[(nx_*ny+ny_)*nz+nz_]) f(idx);
        }
    }


    int flat(int nx, int ny, int nz, float cell, const float * pX, const float * pY, const float * pZ, const vector3& lo) const {
        int cx = std::clamp(static_cast<int>(std::floor((*pX-lo.x)/cell)),0,nx-1);
        int cy = std::clamp(static_cast<int>(std::floor((*pY -lo.y)/cell)),0,ny-1);
        int cz = std::clamp(static_cast<int>(std::floor((*pZ-lo.z)/cell)),0,nz-1);
        return (cx*ny+cy)*nz+cz;
    }

    // Generate one uniform random point inside the sphere
    vector3 randomPointInSphere(vector3 & vec,
        float radius,
        std::mt19937& rng,
        std::uniform_real_distribution<float>& uni,
        std::normal_distribution<float>& normal, float & distance)
    {

        // 1. Radius (cube-root of uniform → correct volume density)
        float u = uni(rng);
        distance = radius * std::cbrtf(u);          // or std::pow(u, 1.0/3.0)

        // 2. Random direction (isotropic)
        vector3 n_vec(normal(rng), normal(rng), normal(rng));

        float norm = n_vec.sqlength();
        // Extremely unlikely to be zero, but protect just in case
        if (norm < 1e-12) {
            return randomPointInSphere(vec, radius, rng, uni, normal, distance); // retry
        }

        // 3. Scale + translate
        return {
            vec + n_vec.normalize()*distance
        };
    }
};


#endif //WETSAXS_ATOMISTICMODEL_H
