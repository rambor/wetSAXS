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

#include "AtomisticModel.h"

#include "Waters.h"
#include "water_shell.h"

AtomisticModel::AtomisticModel(std::string filename, bool forceRNA, bool keepWaters) :
            filename(filename),
            forceRNA(forceRNA),
            keepWaters(keepWaters) {

    /*
     * spherical harmonic calculations use centered coordinates
     */
    pdbModel = PDBModel(filename, keepWaters, forceRNA);

//    std::clock_t startTime = std::clock();
//    float mass = pdbModel.calculateMW();
//    double    runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;
//    std::cout << "total time :: " << runtime << " " << mass << std::endl;

    totalAtoms = pdbModel.getTotalCoordinates();

    rvalues.resize(totalAtoms);
    thetas.resize(totalAtoms);
    phis.resize(totalAtoms);

    // totalDistances = totalAtoms*(totalAtoms-1)/2;
    // ex_volume_weights.resize(totalDistances);

    const float * pXv = pdbModel.getCenteredXVec();
    const float * pYv = pdbModel.getCenteredYVec();
    const float * pZv = pdbModel.getCenteredZVec();

    //const int * pAtomicNumbers = pdbModel.getAtomicNumberVec();
    const int * pAtomicNumbers = pdbModel.getAdjustedASFNumbers();

    float xval, yval, zval, rval;

    // calculate spherical coordinates of atoms - must be centered
    for(unsigned int n=0; n < totalAtoms; n++){
        xval = pXv[n];
        yval = pYv[n];
        zval = pZv[n];

        rval = sqrtf(xval*xval + yval*yval + zval*zval);

        rvalues[n] = rval;
        thetas[n] = acosf(zval/rval);
        phis[n] = (xval == 0 || xval == 0.0f) ? 0.0f : atan2f(yval, xval); // need to check if this is a good way to test for 0

        atomList.insert(pAtomicNumbers[n]); // unique atomic numbers
    }

    // need to reassign atomic numbers to include hydrogens

    atomList.insert(99);  // 99 is unique for water
}

/*
 * Total partials per q value is (lmax+1)^2
 * Total partials im array is totalqvalues*(lmax+1)*(lmax+1)
 * A_lm(q_i) index is given by (l*l + (l-m)*q_index
 *
 * Excluded volume is calculated at the same coordinates of the particle's atoms but uses atomic scattering form factor for water
 *
 */
void AtomisticModel::calculatePartialAmplitudes(unsigned int lmax,
                                                unsigned int totalqvalues,
                                                std::vector < float > & qvalues,
                                                bool recalculate){

    pdbModel.writeCenteredCoordinatesToFile("bsa_centered");
    auto totalPDBAtoms = pdbModel.getTotalCoordinates();
    this->lmax = lmax;
    unsigned int bessel_index, q_index=0, lmax_plus_1 = lmax + 1, index;

    if (recalculate){

        sbj.recalculate(totalqvalues, qvalues, &rvalues);
        //sbj_x.recalculate(totalqvalues, qvalues, &bead_rvalues);

    } else {
//        createHCPGridModel(this->getDmax()/(float)lmax_plus_1);
        // SHE depends on atomic coordinates (r, theta, phi)
        this->createSphericalHarmonics(lmax);
        //this->createSphericalHarmonicsEx();

        // SBJ depends on qvalues and r
        this->createSphericalBessels(totalqvalues, qvalues);
        //this->createSphericalBesselsEx(totalqvalues, qvalues);
    }

    // if (neighbors.empty()){
         this->createNeighborhoods(10.1*10.1);
    // }

    // for each q-value calculate atomic scattering form factor
    // need to assemble into partial sums at each q-value?

    std::map<int, float> asf; // for each q-value, this is calculated over the number of unique atoms

    for(auto & atom : atomList){ // water is element 99
        asf.insert(std::pair<int, float>(atom, 0.0f));
    }

    // total partial amplitudes is (lmax+1)^2 per q-value
    almRealSum.resize((lmax+1)*(lmax+1)*totalqvalues);
    almImagSum.resize((lmax+1)*(lmax+1)*totalqvalues);

    almRealExcludedVolSum.resize((lmax+1)*(lmax+1)*totalqvalues); // partial sums
    almImagExcludedVolSum.resize((lmax+1)*(lmax+1)*totalqvalues);

    float * const ptrRAlm = (totalAtoms != 0) ? almRealSum.data() : nullptr;
    float * const ptrIAlm = (totalAtoms != 0) ? almImagSum.data() : nullptr;

    float * const ptrREx = (totalAtoms != 0) ? almRealExcludedVolSum.data() : nullptr;
    float * const ptrIEx = (totalAtoms != 0) ? almImagExcludedVolSum.data() : nullptr;

    float *pYLM, *pYLMI, bessel_product, bessel_product_Ex, waterASF, * pBessel;
    float *pYLMEx, *pYLMExI, *pBesselEx;

    float partialSumRAp, partialSumIAp, partialSumREx, partialSumIEx;
    float sumExCrossTerm, sumExSelfTerm, weight;

    auto pSBJ = sbj.getpSphericalBessels();
    //auto pAtomicNumbers = pdbModel.getAtomicNumberVec();
    auto pAtomicNumbers = pdbModel.getAdjustedASFNumbers(); // does not include the hydration

    //this->setFractionalVolume(29.9); // set fractional volume of the waters in excluded volume

    auto pReal = she.getpDataYLMReal();
    auto pImag = she.getpDataYLMImag();

    std::clock_t startTime = std::clock();
    //double runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;
    unsigned int lmax_plus_1_totalAtoms = lmax_plus_1*totalAtoms;
    unsigned int long ylm_index, base_index, temp_index;

    unsigned int totalWeights = neighbor_weights.size();
    std::vector<float> bessel_weights(totalWeights);
    auto * pRealYLMEx = (float *)_aligned_malloc(sizeof(float)*totalWeights, 16); // aligned memory
    auto * pImagYLMEx = (float *)_aligned_malloc(sizeof(float)*totalWeights, 16); // aligned memory
    auto * pWeightedBessels = (float *)_aligned_malloc(sizeof(float)*totalWeights, 16); // aligned memory
    auto * pQRvalues = (float *)_aligned_malloc(sizeof(float)*totalWeights, 16); // aligned memory
    float qr;

    for(auto & qvalue : qvalues){

        for(auto & ff : asf){ // calculate atomic scattering form factors for given qvalue
            ff.second = functions_WETSAXS::asf(ff.first, qvalue);
        }

        waterASF = asf[99];//functions_WETSAXS::asf(99, qvalue);

        //  this->setFractionalVolumeAtQ(qvalue);

        // populate qr values
        unsigned int ll_index = (lmax+1)*(lmax+1)*q_index;// + l*l + l + m;

        float intensity = calculateStructureFactorWaterAtQ(qvalue);
        float temp, bessel_products, bessel_products_real, bessel_products_imag, * pTemp;

        for(int l=0; l<=lmax; l++){

            bessel_index = lmax_plus_1_totalAtoms*q_index + totalAtoms*l;
            base_index = (l*l)*totalAtoms; // always positive

            // calculate spherical bessels for given q, l over all atoms
            // populateSphBesselsArray(l, qr, bessel_at_ql);
            /*
             * for m = 0
             */
            // for given (l,m) calculate Spherical Harmonic
            partialSumRAp = 0.0f;

            // excluded volume term based on atomic coordinates
            partialSumREx = 0.0f;
            partialSumIEx = 0.0f;

            sumExCrossTerm = 0.0f;
            sumExSelfTerm = 0.0f;

            ylm_index = base_index + l * totalAtoms; // for m = 0

            /*
             * for m = 0
             */
            startTime = std::clock();
            // neighbor_indices and neighbor_weights are long vectors > totalAtoms
            for (unsigned int n = 0; n < totalWeights; n++) { // over all atoms including hydration waters
                    auto *pIndex = &neighbor_indices[n];
                    pWeightedBessels[n] = neighbor_weights[n] * pSBJ[bessel_index + *pIndex];
                    pRealYLMEx[n] = pWeightedBessels[n] * pReal[ylm_index + *pIndex];
            }

            // two cross term sums - one is Ap-Am and the other is Ax-Am where real and imaginary terms need to be kept separate
            // Am cross-terms will be generated using a different algorithm (Mask class)
            // Ax and Ap can be summed directly

            for (unsigned int n = 0; n < totalAtoms; n++) { // over all atoms
                pBessel = &pSBJ[bessel_index + n];
                //bessel_product = asf[pAtomicNumbers[n]] * *pBessel; // for each m value, using same bessel values at l
                //bessel_product_Ex = (waterASF - fractionOccAtoms[n]) * *pBessel; // for estimating excluded volume using water ASF

                // real term for spherical harmonic
                pYLM = &pReal[ylm_index];
                //partialSumRAp += bessel_product * *pYLM;
                partialSumRAp += ( n < totalPDBAtoms  ? asf[pAtomicNumbers[n]] : waterASF)* *pBessel * *pYLM;
                //imaginary term is zero since sin(m*theta) with m = 0 is always zero.
                //partialSumREx += bessel_product_Ex * *pYLM;  // excluded volume, real

                //pYLM = &pImag[ylm_index];
                //partialSumIAp += bessel_product * *pYLM;    // particle, imaginary
                //partialSumIEx += bessel_product_Ex * *pYLM; // excluded volume, imaginary

                // excluded volume must be a separate calculation that uses the weights and waterASF
                temp = *pBessel * *pYLM; // shall be scaled by Water ASF
                sumExSelfTerm += temp*temp; // self term
                // now do the weighted cross-term
                bessel_products = 0;
                for (unsigned int v =neighbors_offset[n]; v < neighbors_offset[n+1]; v++) { // over all atoms including hydration waters
                    bessel_products += pRealYLMEx[v];
                }

                sumExCrossTerm += temp * bessel_products;
                ++ylm_index;
            }
            //std::cout << qvalue << " " << l << " " << partialSumAEx << std::endl;
            // 0.000015 sec without internal loop
            // 0.0061
            // weighted_sum = addUp(totalWeights, pRealYLMEx);
            // std::cout << " RUN TIME " << (std::clock() - startTime)/(double) CLOCKS_PER_SEC << std::endl;

            int partials_index = ll_index + l*l + l;

            // particle scattering for m=0
            ptrRAlm[partials_index] = partialSumRAp; // need to keep the partials for fitting
            ptrIAlm[partials_index] = 0;// sum is zero due to sin(0) = 0 where m = 0;

            // excluded volume scattering for m=0
            // sum f^2 + 2 * sum f_i*f_j for j > i
            sumExCrossTerm *= 2;
            sumExCrossTerm += sumExSelfTerm; // excluded volume I(q) for m = 0 at given l
            ptrREx[partials_index] = waterASF*waterASF*sumExCrossTerm;// partialSumREx; // uses volume fraction with assigned waters

            //intensity += waterASF*waterASF*sumExCrossTerm;
            ptrIEx[partials_index] = 0;; // final should be negative of the amplitude (Babinet's)

            // biggest time saving will be here, need to exploit the symmetry of spherical hamormonics (for m = -m)
            for(int m=1; m<=l; m++) {

                // for given (l,m) calculate Spherical Harmonic
                partialSumRAp = 0.0f;
                partialSumIAp = 0.0f;

                // excluded volume term based on atomic coordinates
                partialSumREx = 0.0f;
                partialSumIEx = 0.0f;

                sumExCrossTerm = 0.0f;
                sumExSelfTerm = 0.0f;

                ylm_index = base_index + (l + m) * totalAtoms; // incorrect

                for (unsigned int n = 0; n < totalWeights; n++) { // over all atoms including hydration waters
                    temp_index = ylm_index + neighbor_indices[n];
                    pTemp = &pWeightedBessels[n];

                    pRealYLMEx[n] = *pTemp * pReal[temp_index]; // each value is specified by l and m
                    pImagYLMEx[n] = *pTemp * pImag[temp_index];
                }


                // for each (l,m) calculate partial sum at given q value
                float real_product = 0, imaginary_product = 0;
                for (unsigned int n = 0; n < totalAtoms; n++) { // over all atoms

                    pBessel = &pSBJ[bessel_index + n];
                    bessel_product = ( n < totalPDBAtoms  ? asf[pAtomicNumbers[n]] : waterASF)* *pBessel;
                    //bessel_product_Ex = (waterASF - fractionOccAtoms[n]) * *pBessel; // for estimating excluded volume using water ASF
                    //bessel_product_Ex = waterASF * *pBessel; // for estimating excluded volume using water ASF

                    // in using the radial distribution function, we need to discount the waters based on cross-term correspondences, so waters far from each other have no correlation
                    // real term
                    pYLM = &pReal[ylm_index];
                    //real_product = bessel_product * *pYLM;     // particle, real
                    //partialSumRAp += real_product;     // particle, real
                    partialSumRAp += bessel_product* *pYLM;
                    //partialSumREx += bessel_product_Ex * *pYLM;  // excluded volume, real

                    // imaginary term
                    pYLMI = &pImag[ylm_index];
                    //imaginary_product = bessel_product * *pYLMI;
                    partialSumIAp += bessel_product* *pYLMI;    // particle, imaginary
                    //partialSumIEx += bessel_product_Ex * *pYLMI; // excluded volume, imaginary

                    // r = 0, self term for each position
                    sumExSelfTerm += (*pYLM * *pYLM + *pYLMI * *pYLMI)* *pBessel* *pBessel;

                    bessel_products_real = 0;
                    bessel_products_imag = 0;
                    for (unsigned int v =neighbors_offset[n]; v < neighbors_offset[n+1]; v++) { // over all atoms including hydration waters
                        bessel_products_real += pRealYLMEx[v];
                        bessel_products_imag += pImagYLMEx[v];
                    }
                    //complete the cross-term products and add to running total
                    real_product += bessel_products_real * *pBessel * *pYLM;
                    imaginary_product += bessel_products_imag * *pBessel * *pYLMI;

                    ++ylm_index;
                }

                int update_index = partials_index + m;

                // particle scattering
                // for m > 0
                ptrRAlm[update_index] = partialSumRAp;
                ptrIAlm[update_index] = partialSumIAp;

                // excluded volume
                // intensity is the sum of the squared amplitudes at each l,
                ptrREx[update_index] = waterASF*waterASF*(sumExSelfTerm + 2*(real_product + imaginary_product)); // update at l, m
                ptrIEx[update_index] = 0;//partialSumIEx; //

                //intensity += 2*ptrREx[update_index]; // 2x is for m and -m

                // partials_index = ll_index + l*l + l;
                // for m < 0
                update_index = partials_index - m;
                ptrREx[update_index] = waterASF*waterASF*(sumExSelfTerm + 2*(real_product + imaginary_product)); // The squared sum is same for -m and m

                if (m & 1) { // if m is odd, flip sign of the real part
                    // std::cout << count << " ( " << (partials_index -m) <<  " ) "  << ptrRAlm[partials_index -m] << " " << (-1*partialSumRAp) << std::endl;
                    // particle scattering
                    ptrRAlm[update_index] = -partialSumRAp;
                    ptrIAlm[update_index] = partialSumIAp;

                    // excluded volume scattering  - remember to flip sign for Babinet's
                    // ptrREx[update_index] = -partialSumREx; // uses volume fraction with assigned waters
                    ptrIEx[update_index] = partialSumIEx;  //

                } else { // if m is even, flip sign of the Imaginary part
                    // std::cout << count << " ( " << (partials_index -m) <<  " ) "  << ptrIAlm[partials_index -m] << " " << (-1*partialSumIAp) << std::endl;
                    // particle scattering
                    ptrRAlm[update_index] = partialSumRAp;
                    ptrIAlm[update_index] = -partialSumIAp;

                    // excluded volume scattering - remember to flip sign for Babinet's
                    //ptrREx[update_index] = partialSumREx;  // uses volume fraction with assigned water
                    ptrIEx[update_index] = -partialSumIEx; //
                }
            }
        }

        //intensity *= waterASF*waterASF;//*(1 + 4*M_PI*0.0334*addUp(totalWeights, pQRvalues));

        // print out
       std::cout << qvalue << " " << intensity << std::endl;
        ++q_index;
    }

exit(0);
    //calculateExVolPartialAmplitudes(lmax, totalqvalues, qvalues, recalculate);
    logger("------------------------------","------------------------------");

    // calculate I(0) for protein and excluded volume

    try {
       // pAtomicNumbers = pdbModel.getAtomicNumberVec();
        setMatchPointScale(pAtomicNumbers);
    } catch (std::overflow_error &e ){
        std::cout << e.what() << " :: ";
    }

    double runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;

    //    this->writeInVacuoModel(lmax, totalqvalues, qvalues);
    logger("Total TIME (Model)", formatNumber((float)runtime,8));
}


void AtomisticModel::calculateExVolPartialAmplitudes(unsigned int lmax,
                                                     unsigned int totalqvalues,
                                                     std::vector < float > & qvalues,
                                                     bool recalculate){

    unsigned int bessel_index, q_index=0, lmax_plus_1 = lmax + 1;
    float *pYLM, bessel_product_Ex, waterASF, * pBessel;
    float partialSumREx, partialSumIEx;

    almRealExcludedVolSum.resize((lmax+1)*(lmax+1)*totalqvalues);; // partial sums
    almImagExcludedVolSum.resize((lmax+1)*(lmax+1)*totalqvalues);

    float * const ptrREx = (beads.size() != 0) ? almRealExcludedVolSum.data() : nullptr;
    float * const ptrIEx = (beads.size() != 0) ? almImagExcludedVolSum.data() : nullptr;

    unsigned int totalBeads = beads.size();
    unsigned int lmax_plus_1_totalBeads = lmax_plus_1*totalBeads;
    unsigned int long ylm_index, base_index;

    auto pSBJ = sbj_x.getpSphericalBessels();

    auto pReal = she_x.getpDataYLMReal();
    auto pImag = she_x.getpDataYLMImag();

//    if (recalculate){
//        sbj_x.recalculate(totalqvalues, qvalues, &bead_rvalues);
//    } else {
//        // SHE depend on atomic coordinates (r, theta, phi)
//        this->createSphericalHarmonicsEx();
//
//        // SBJ depends on qvalues and r
//        this->createSphericalBesselsEx(totalqvalues, qvalues);
//    }

    float inv = 1.0f/(4*M_PI);
    float diff = 0.6366/1.75;
    float factor = 29.9f*(1 - diff*diff*diff);

    for(auto & qvalue : qvalues){

        float qval2 = qvalue*qvalue;
        float cf = 0.334f*factor * std::expf(-std::cbrt(factor*factor)*inv*qval2);

        waterASF = functions_WETSAXS::asf(99, qvalue) - cf;

        // populate qr values
        unsigned int ll_index = (lmax+1) * (lmax+1) * q_index;// + l*l + l + m;

        for(int l=0; l<=lmax; l++){

            bessel_index = lmax_plus_1_totalBeads*q_index + totalBeads*l;
            base_index = (l*l)*totalBeads; // always positive

            // calculate spherical bessels for given q, l over all atoms
            // populateSphBesselsArray(l, qr, bessel_at_ql);
            /*
             * for m = 0
             */
            // for given (l,m) calculate Spherical Harmonic
            // excluded volume term based on atomic coordinates
            partialSumREx = 0.0f;
            partialSumIEx = 0.0f;

            ylm_index = base_index + l * totalBeads; // for m = 0

            for (unsigned int n = 0; n < totalBeads; n++) { // over all atoms
                pBessel = &pSBJ[bessel_index + n];
                bessel_product_Ex = *pBessel; // for estimating excluded volume using water ASF

                // real term
                partialSumREx += bessel_product_Ex * pReal[ylm_index];  // excluded volume, real

                //imaginary term
                partialSumIEx += bessel_product_Ex * pImag[ylm_index]; // excluded volume, imaginary

                ++ylm_index;
            }

            int partials_index = ll_index + l*l + l;

            // excluded volume scattering for m=0
            ptrREx[partials_index] = waterASF * partialSumREx; // uses volume fraction with assigned waters
            ptrIEx[partials_index] = waterASF * partialSumIEx; // final should be negative of the amplitude (Babinet's)

            // biggest time saving will be here, need to exploit the symmetry of spherical hamormonics (for m = -m)
            for(int m=1; m<=l; m++) {
                // for given (l,m) calculate Spherical Harmonic
                // excluded volume term based on atomic coordinates
                partialSumREx = 0.0f;
                partialSumIEx = 0.0f;

                ylm_index = base_index + (l + m) * totalBeads; //incorrect

                // for each (l,m) calculate partial sum at given q value
                for (unsigned int n = 0; n < totalBeads; n++) { // over all atoms

                    pBessel = &pSBJ[bessel_index + n];
                    bessel_product_Ex = *pBessel; // for estimating excluded volume using water ASF

                    // real term
                    pYLM = &pReal[ylm_index];
                    partialSumREx += bessel_product_Ex * *pYLM;  // excluded volume, real

                    //imaginary term
                    pYLM = &pImag[ylm_index];
                    partialSumIEx += bessel_product_Ex * *pYLM; // excluded volume, imaginary

                    ++ylm_index;
                }

                // excluded volume
                partialSumREx *= waterASF;
                partialSumIEx *= waterASF;

                int update_index = partials_index + m;
                ptrREx[update_index] = partialSumREx; // uses volume fraction with assigned waters
                ptrIEx[update_index] = partialSumIEx; // final should be negative of the amplitude (Babinet's)

                // update negative m values
                update_index = partials_index - m;

                if (m & 1) {// if m is odd, flip sign of the real part
                    // std::cout << count << " ( " << (partials_index -m) <<  " ) "  << ptrRAlm[partials_index -m] << " " << (-1*partialSumRAp) << std::endl;
                    ptrREx[update_index] = -partialSumREx; // uses volume fraction with assigned waters
                    ptrIEx[update_index] = partialSumIEx;  // final should be negative of the amplitude (Babinet's)

                } else { // if m is even, flip sign of the Imaginary part
                    // std::cout << count << " ( " << (partials_index -m) <<  " ) "  << ptrIAlm[partials_index -m] << " " << (-1*partialSumIAp) << std::endl;
                    ptrREx[update_index] = partialSumREx;  // uses volume fraction with assigned waters
                    ptrIEx[update_index] = -partialSumIEx; // final should be negative of the amplitude (Babinet's)
                }
            }
        }
        ++q_index;
    }
}

/*
 * Model includes hydrogens (using unified form factor) but ignores hydration waters
 * average electron density of a protein should be around 0.4 electrons per cubic Angstrom
 *
 */
void AtomisticModel::setMatchPointScale(int * pAtomicNumbers){

    forwardScatteringParticle = 0.0f;
    //float * gaussian_radii = pdbModel.getAtomicGaussianRadii();
    auto totalParticleAtoms = pdbModel.getTotalCoordinates();
    float diff;

    // this should include hydrogens since using the unified form factors
    for(unsigned int n=0; n < totalParticleAtoms; n++){ // over all atoms
        diff = functions_WETSAXS::asf_at_q_zero(pAtomicNumbers[n]);
        forwardScatteringParticle += diff*diff;
    }

    for(unsigned int n=0; n < totalParticleAtoms; n++){ // over all atoms
        unsigned int next_j = n+1;
        float sumIt = 0;
        for(unsigned int j=next_j; j < totalParticleAtoms; j++){ // over all atoms
            sumIt += functions_WETSAXS::asf_at_q_zero(pAtomicNumbers[j]);
        }
        forwardScatteringParticle += 2*functions_WETSAXS::asf_at_q_zero(pAtomicNumbers[n])*sumIt;
    }


    /*
     * calculate forward scattering of excluded volume
     *
     */
    // water_asf at zero will be 10 electrons
    const float water_asf = functions_WETSAXS::asf_at_q_zero(99);

    // add the waters from the hydration.
    auto totalHydratingWaters = totalAtoms - totalParticleAtoms;
    forwardScatteringParticle += totalHydratingWaters*water_asf*water_asf;

    // const float pw_vol = 29.9*0.334; // 10 e_n
    // float factor;
    // for(unsigned int i=0; i<totalParticleAtoms; i++){ // should be adjusted to include the waters that form the hydration
    //     diff = gaussian_radii[i]/1.75; // 1.75 is gaussian sphere radius of water
    //
    //     factor = pw_vol*(1-diff*diff*diff);
    //     diff = water_asf;// - factor; // if radius is same as water, no correction
    //     forwardScatteringExVolume += diff*diff; // normalized to volume of water, looking at fractions of water volume
    // }
    //
    // for(unsigned int n=0; n < totalParticleAtoms; n++){ // over all atoms
    //     unsigned int next_j = n+1;
    //     diff = gaussian_radii[n]/1.75;
    //     factor = pw_vol*(1-diff*diff*diff);
    //     diff = (water_asf - factor);
    //     diff = water_asf;
    //     float sumIt = 0;
    //
    //     for(unsigned int j=next_j; j < totalParticleAtoms; j++){ // over all atoms
    //         float vdiff = gaussian_radii[j]/1.75;
    //         factor = pw_vol*(1-vdiff*vdiff*vdiff);
    //         sumIt += (water_asf - factor);
    //     }
    //     forwardScatteringExVolume += 2*diff*sumIt;
    // }

forwardScatteringExVolume = totalAtoms*totalAtoms*water_asf*water_asf;
    // forwardScatteringExVolume = 0.0f;
    // for (auto & weight : neighbor_weights) {
    //     forwardScatteringExVolume += weight;
    // }
    // forwardScatteringExVolume *= 2;
    // forwardScatteringExVolume += totalAtoms;
    // forwardScatteringExVolume *= water_asf*water_asf;

    //forwardScatteringExVolume += 377*377*pw_vol*pw_vol;
    //forwardScatteringExVolume = forwardScatteringH2Os*0.334*0.334 + 1067*1067*pw_vol*pw_vol;
    SASTOOLS_UTILS_H::logger("Forward Scattering Particle", formatNumber(forwardScatteringParticle, 1));
    SASTOOLS_UTILS_H::logger("Forward Scattering ExVol", formatNumber(forwardScatteringExVolume, 1));

    forwardScatteringParticle = std::sqrtf(forwardScatteringParticle);
    forwardScatteringExVolume = std::sqrtf(forwardScatteringExVolume);

    if (totalAtoms < 1 || forwardScatteringExVolume == 0)
        throw std::overflow_error("Divide by zero exception");

    SASTOOLS_UTILS_H::logger("density (en-/Angstrom^3)", formatNumber(forwardScatteringParticle/pdbModel.getVolume(), 3));
    //SASTOOLS_UTILS_H::logger("Total H2Os in Ex Volume", formatNumber(totalH2Os, 0));

    // proteins have a match point of 0.42 e_n per Angstrom3
    // nucleics have a match point of 0.55 e_n per Angstrom3
    matchPoint = forwardScatteringParticle/forwardScatteringExVolume;

}


/*
 * Model includes hydrogens implicitly through the unified form factors
 */
void AtomisticModel::setMatchPointScale(unsigned int totalHydrationWaters){

    //float water_asf = functions_WETSAXS::asf_at_q_zero(99)*totalHydrationWaters;

    if (totalAtoms < 1 || forwardScatteringExVolume == 0)
        throw std::overflow_error("Divide by zero exception");

//    matchPoint = (forwardScatteringParticle + water_asf)/( (forwardScatteringExVolume + water_asf)); // divide by zero possible
//    std::cout << matchPoint << " " << (forwardScatteringParticle)/( (forwardScatteringExVolume)) << std::endl;
    matchPoint = (forwardScatteringParticle)/( (forwardScatteringExVolume)); // divide by zero possible
}


void AtomisticModel::writeInVacuoModel(unsigned int lmax, unsigned int totalqvalues, std::vector < float > & qvalues){

    std::string filename = pdbModel.getFileStemName();  //includes extension

    std::string nameOf = filename + "_vac.int";
    FILE * pFile = fopen(nameOf.c_str(), "w");

    this->calculatePartialAmplitudesNoExcludedVolume(lmax, totalqvalues, qvalues);

    float * const ptrR = (totalAtoms != 0) ? almRealSum.data() : nullptr;
    float * const ptrI = (totalAtoms != 0) ? almImagSum.data() : nullptr;

    unsigned int count = 0;
    int q_index = 1;

    float sumR, sumI, pi_constant=16*M_PI*M_PI, i_calc;

    for(auto & qvalue : qvalues){

        i_calc = 0;

        for(int l=0; l<=lmax; l++){

            for(int m=-l; m<=l; m++){
                sumR = ptrR[count];
                sumI = ptrI[count];

                i_calc += (sumR*sumR + sumI*sumI); // sum of the squared amplitudes
                count++;
            }
        }

        i_calc *= pi_constant;
        fprintf(pFile, "%5i %.7f %E \n", q_index, qvalue, i_calc);

        q_index++;
    }

    fclose(pFile);
}


void AtomisticModel::createSphericalHarmonics(unsigned int lmax) {

    this->lmax = lmax;
    she = SphericalHarmonics(lmax, totalAtoms);
    she.populateSHETable(thetas.data(), phis.data());
}


void AtomisticModel::createSphericalHarmonicsEx() {

    she_x = SphericalHarmonics(lmax, beads.size());
    she_x.populateSHETable(bead_thetas.data(), bead_phis.data());

}

/*
 * must be run after SphericalHarmonics
 */
void AtomisticModel::createSphericalBessels(unsigned int totalqvalues, std::vector < float > & qvalues) {

    sbj = SphericalBessels(lmax, totalqvalues, totalAtoms, qvalues, &rvalues);
}


/*
 * must be run after SphericalHarmonics
 */
void AtomisticModel::createSphericalBesselsEx(unsigned int totalqvalues, std::vector < float > & qvalues) {

    sbj_x = SphericalBessels(lmax, totalqvalues, beads.size(), qvalues, &bead_rvalues);
}


const SphericalBessels &AtomisticModel::getSbj() const {
    return sbj;
}


/*
 * Total partials per q value is (lmax+1)^2
 * Total partials im array is totalqvalues*(lmax+1)*(lmax+1)
 * A_lm(q_i) index is given by (l*l + (l-m)*q_index
 *
 * Equivalent to in vacuo scattering
 *
 */
void AtomisticModel::calculatePartialAmplitudesNoExcludedVolume(unsigned int lmax, unsigned int totalqvalues, std::vector < float > & qvalues){

    this->createSphericalHarmonics(lmax);
    this->createSphericalBessels(totalqvalues, qvalues);

    // for each q-value calculate atomic scattering form factor
    // need to assemble into partial sums at each q-value?

    //std::vector<float> asf; // for each q-value, this is calculated over the number of unique atoms
    std::map<int, float> asf;

    for(auto & atom : atomList){
        asf.insert(std::pair<int, float>(atom, 0.0f));
    }

    // total partial amplitudes is (lmax+1)^2 per q-value
    almRealSum.resize((lmax+1)*(lmax+1)*totalqvalues);
    almImagSum.resize((lmax+1)*(lmax+1)*totalqvalues);

    float * const ptrRAlm = (totalAtoms != 0) ? almRealSum.data() : nullptr;
    float * const ptrIAlm = (totalAtoms != 0) ? almImagSum.data() : nullptr;

    float bessel_product;

    auto pSBJ = sbj.getpSphericalBessels();
    //auto pAtomicNumbers = pdbModel.getAtomicNumberVec();
    auto pAtomicNumbers = pdbModel.getAdjustedASFNumbers();

    auto pRealYlm = she.getpDataYLMReal();
    auto pImagYlm = she.getpDataYLMImag();

    float preSumR, preSumI;

    unsigned int bessel_index, count = 0, q_index=0;
    unsigned int long ylm_index, base_index;

    std::clock_t startTime = std::clock();

    for(auto & qvalue : qvalues){

        for(auto & ff : asf){ // calculate atomic scattering form factors for given qvalue
            ff.second = functions_WETSAXS::asf(ff.first, qvalue);
        }

        for(int l=0; l<=lmax; l++){

            bessel_index = (lmax+1)*totalAtoms*q_index + totalAtoms*l;
            base_index = (l*l)*totalAtoms;

            // calculate spherical bessel for given q and l for all atoms

            for(int m=-l; m<=l; m++){ // for given (l,m) calculate Spherical Harmonic
                // over all atoms
                preSumI = 0.0f;
                preSumR = 0.0f;

                ylm_index = base_index + (l+m)*totalAtoms; // incorrect

                // for each (l,m) calculate partial sum at given q value
                for(unsigned int n=0; n<totalAtoms; n++){

                    bessel_product = asf[ pAtomicNumbers[n] ] * pSBJ[bessel_index + n];

                    preSumR += bessel_product * pRealYlm[ylm_index];// real term
                    preSumI += bessel_product * pImagYlm[ylm_index];//imaginary term

                    ylm_index++;
                }

                ptrRAlm[count] = preSumR;
                ptrIAlm[count] = preSumI;

                count++;
            }
        }
        q_index++;
    }

    double    runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;
    logger("Partials TIME", formatNumber((float)runtime,8));
}

void AtomisticModel::addWatersToAtomicModel(Waters &waterModel ) {

    // add waters to the atomic model
    auto totalWaters = waterModel.getTotalWaters();

    auto wrvalues = waterModel.getRvalues();
    auto wphiValues = waterModel.getPhiValues();
    auto wthetaValues = waterModel.getThetaValues();

    unsigned int startHere = totalAtoms;

    if ( totalWaters < 1 ) {
        throw std::invalid_argument( "Water model not initialized" );
    }

    // resize populated spherical coordinates of PDB model to add the hydration

    totalAtoms += totalWaters;
    rvalues.resize(totalAtoms);
    thetas.resize(totalAtoms);
    phis.resize(totalAtoms);

    for (unsigned int i = 0; i<totalWaters; i++) {
        rvalues[startHere] = wrvalues[i];
        thetas[startHere] = wthetaValues[i];
        phis[startHere] = wphiValues[i];
        startHere += 1;
    }

    // atoms that make up the particle model end at => pdbModel.getTotalCoordinates()
    // converting spherical coordinates to harmonics is simple call to SHE and Bessel functions

    // could add the waters that form the hydration to the coordinate model?
    totalDistances = totalAtoms*(totalAtoms-1)/2;
    ex_volume_weights.resize(totalDistances, 0);

    this->populateExVolumeWeights(waterModel);
}

/*
 * Use this function to determine weights between atoms based on distances for excluded volume calculation
 *
 */
void AtomisticModel::populateExVolumeWeights(Waters &waterModel){

    auto totalWaters = waterModel.getTotalWaters();
    auto totalPDBAtoms = pdbModel.getTotalCoordinates();

    const float * pX = pdbModel.getCenteredXVec();
    const float * pY = pdbModel.getCenteredYVec();
    const float * pZ = pdbModel.getCenteredZVec();

    //float inv_delta = (float)lmax/ pdbModel.getDmax();
    float dis;
    float a_1 = 2.4 ;
    float r_1 = 2.8 ;
    float s_1 = 0.5/(0.28*0.28);
    float a_2 = 0.5 ;
    float r_2 = 4.5 ;
    float s_2 = 0.5/(0.55*0.55);
    float a_3 = 0.2 ;
    float r_3 = 7.1 ;
    float s_3 = 0.5/(0.9*0.9) ;
    float k=8;
    float normalization = 1.0/17.5;

    unsigned int counter = 0;
   // unsigned int num;
    float minLimit = 2.49*2.49;
    float maxLimit = 10.1*10.1;

    // using particle to particle atoms
    // then water to particle atoms
    // finally water to water atoms
    for(unsigned int i=0; i<totalPDBAtoms; i++){
        vector3 vec1(pX[i], pY[i], pZ[i]);

        for(unsigned int j=(i+1); j<totalPDBAtoms; j++){
            dis = (vector3(pX[j], pY[j], pZ[j])-vec1).sqlength();
            if (dis > minLimit && dis < maxLimit) {
                dis = sqrt(dis);
                // P(r) distribution of water => r^2*( g(r)-1 )
                // maximum value at P(2.8) is 1, normalized to peak value
                ex_volume_weights[counter] = normalization*dis*dis*((1 + a_1*expf(-std::powf(dis-r_1,2)*s_1) + a_2*expf(-std::powf(dis-r_2,2)*s_2) + a_3*expf(-std::powf(dis-r_3,2)*s_3))*(1-expf(-k*std::powf(dis-2.3,2))) - 1);
            }
            counter++;
        }
    }

    // now do the water-to-particle model cross-terms
    auto waters = waterModel.getHydration();
    Coords * pC2;

    for(unsigned int i=0; i < totalPDBAtoms; i++){

        vector3 vec1(pX[i], pY[i], pZ[i]);

        for(unsigned int j=0; j<totalWaters; j++){

            pC2 = &waters[j];
            dis = (vector3(pC2->x, pC2->y, pC2->z)-vec1).sqlength();

            if (dis > minLimit && dis < maxLimit) {
                dis = sqrt(dis);
                ex_volume_weights[counter] = normalization*dis*dis*((1 + a_1*expf(-std::powf(dis-r_1,2)*s_1) + a_2*expf(-std::powf(dis-r_2,2)*s_2) + a_3*expf(-std::powf(dis-r_3,2)*s_3))*(1-expf(-k*std::powf(dis-2.3,2))) - 1);
            }
            counter++;
        }
    }


    for(unsigned int i=0; i < totalWaters; i++){

        const Coords *pC = &waters[i];
        vector3 vec1(pC->x, pC->y, pC->z);

        for(unsigned int j=i+1; j<totalWaters; j++){
            pC2 = &waters[j];
            dis = (vector3(pC2->x, pC2->y, pC2->z)-vec1).sqlength();

            if (dis > minLimit && dis < maxLimit) {
                dis = sqrt(dis);
                ex_volume_weights[counter] = normalization*dis*dis*((1 + a_1*expf(-std::powf(dis-r_1,2)*s_1) + a_2*expf(-std::powf(dis-r_2,2)*s_2) + a_3*expf(-std::powf(dis-r_3,2)*s_3))*(1-expf(-k*std::powf(dis-2.3,2))) - 1);
            }
            counter++;
        }
    }
}


/*
 * bessel_alq aligned memory
 */
void AtomisticModel::populateSphBesselsArray(int l, std::vector<float> & qr, float * bessel_alq) {

    float * const pQR = qr.data();

    for (int i=0 ; i < totalAtoms; i++) {
        bessel_alq[ i ] =  boost::math::sph_bessel(l, pQR[i]);
    }

}

/*
 * cutoff is the squared distance limit
 * calculate unique index value if neighbors using formula
 * row*N - row*(row+1)/2 -row + col -1
 * row is always less than col
 * nearestNeighbors is a set so look up should be fast
 *
 * correlation length of water at ambient is around 1.3 Angstroms (Hura and Huang 2009 and 2010)
 *
 */
void AtomisticModel::createNeighborhoods(float cutoff){

    const float * pRv = rvalues.data();//pdbModel.getCenteredXVec();
    const float * pThetav = thetas.data();//pdbModel.getCenteredYVec();
    const float * pPhiv = phis.data();//pdbModel.getCenteredZVec();

    neighbors_offset.resize(totalAtoms+1);
    neighbors_offset[0] = 0;

    paired_r_values.clear();

    float minLimit = 2.57*2.57;
    float r1, r2, r1_2, r2_2, r12, theta1, theta2, diff, cosphi, phi1;

    // calculate spherical coordinates of atoms - must be centered
    unsigned int totalweights = 0;
    for(unsigned int row=0; row < totalAtoms; row++){

        r1 = pRv[row];
        theta1 = pThetav[row];
        phi1 = pPhiv[row];

        r1_2 = r1*r1;
        unsigned int totalNeighbors = 0;
        unsigned int col = 0;
        float weightsum = 0;
        for(; col < row; col++){

            r2 = pRv[col];
            r2_2 = r1_2 + r2*r2;

            r12 = 2*r1*r2;

            if (r2_2 - r12 < cutoff) { // calculate difference in spherical coordinates
                theta2 = pThetav[col];
                cosphi = cos(phi1 - pPhiv[col]);
                diff = r2_2 - r12*0.5*(cos(theta1 - theta2)*(1 + cosphi) + cos(theta1 + theta2)*(1 - cosphi));

                if (diff > minLimit && diff < cutoff) {
                    neighbor_indices.push_back(col);   // all neighbors concatenated
                    paired_r_values.push_back(sqrt(diff));
                    neighbor_weights.push_back(calculateWaterWeight(paired_r_values.back()));
                    // if (neighbor_weights.back() < 0) {
                    //     std::cout << "neighbor weight " << sqrtf(diff) << " " << neighbor_weights.back() << std::endl;
                    //     exit(0);
                    // }
                    weightsum += neighbor_weights.back();
                    totalNeighbors++;
                    totalweights++;
                }
            }
        }
        col++;

        for(; col < totalAtoms; col++){

            r2 = pRv[col];
            r2_2 = r1_2 + r2*r2;

            r12 = 2*r1*r2;

            if (r2_2 - r12 < cutoff) { // calculate difference in spherical coordinates
                theta2 = pThetav[col];
                cosphi = cos(phi1 - pPhiv[col]);
                diff = r2_2 - r12*0.5*(cos(theta1 - theta2)*(1 + cosphi) + cos(theta1 + theta2)*(1 - cosphi));

                if (diff > minLimit && diff < cutoff) {
                    neighbor_indices.push_back(col);   // all neighbors concatenated
                    paired_r_values.push_back(sqrt(diff));
                    neighbor_weights.push_back(calculateWaterWeight(paired_r_values.back()));
                    weightsum += neighbor_weights.back();
                    totalNeighbors++;
                    totalweights++;
                }
            }
        }
        neighbors_offset[row+1] = totalNeighbors + neighbors_offset[row]; // record total neighbors per atom
        //std::cout << "Total neighbors: " << row << " " << totalNeighbors << " " << weightsum << std::endl;
    }
}

void AtomisticModel::createHCPGridModel(float bin_width){

    const float * pXv = pdbModel.getCenteredXVec();
    const float * pYv = pdbModel.getCenteredYVec();
    const float * pZv = pdbModel.getCenteredZVec();

    pdbModel.writeCenteredCoordinatesToFile("centered");

    float radius = 2.8/2; // spacing between points in lattice, not radius
    float sqrt3 = sqrtf(3);
    float inv_rsqrt3= 1.0f/(radius * sqrt3);

    // create HCP grid model
    float inv_radius = 1.0/radius;
    const vector3 a1 = vector3(0.5*radius, -0.5*sqrt3*radius, 0);
    const vector3 a2 = vector3(0.5*radius,  0.5*sqrt3*radius, 0);
    const vector3 a3 = vector3(0, 0, (float)(2.0f*sqrt(6)/3.0f*radius)); // seems fine, spreads out

    const vector3 basis2 = a1*(2.0f/3.0f) + a2*(1.0f/3.0f) + a3/2.0f;
    logger("=> CREATING UNIVERSE","SPHERICAL SEARCH SPACE");

    float z_point = 1.5f/sqrtf(6.0)*inv_radius;
    float inv2sqrt6r = 0.5/sqrtf(6.0)*inv_radius;
    float invsqrt3r = 1.0f/sqrtf(3.0)*inv_radius;
    float sqrt3inv3r = sqrtf(3.0)/3.0f*radius;
    float sqrt6r = sqrtf(6.0)*radius;
    int n1_index, n2_index, n3_index, n1_index_b2, n2_index_b2, n3_index_b2;

    std::set<std::string> bead_indices;
    std::string index;

    // for a given atomic position, need all indices within 2.8 Angstrom?

    for(unsigned int i=0; i<totalAtoms; i++){

        const float * pX = &pXv[i];
        const float * pY = &pYv[i];

        n3_index = (int)std::round(z_point*pZv[i]);
        n2_index = (int)std::round( (sqrt3* *pX + *pY)*inv_rsqrt3);
        n1_index = (int)std::round( (sqrt3* *pX - *pY)*inv_rsqrt3);

        index = std::to_string(n1_index).append("_").append(std::to_string(n2_index)).append("_").append(std::to_string(n3_index));
        bead_indices.insert(index);

        // for a given index triple, nearest points within r given by
        // 6 closest points in the plane
        // (0,1,0)
        // (1,1,0)
        // (1,0,0)
        // (0,-1,0)
        // (-1,-1,0)
        // (-1,0,0)

        // Above and below the plane
        // (0,0,0) - basis
        // (1,1,0) - basis
        // (1,0,0) - basis

        // (0,0,0) + basis
        // (-1,0,0) + basis
        // (-1,-1,0) + basis

        // avoid repeats

        // (0,1,0)
        index = std::to_string(n1_index).append("_").append(std::to_string(n2_index + 1)).append("_").append(std::to_string(n3_index));
        bead_indices.insert(index);

        // (1,1,0)
        index = std::to_string(n1_index + 1).append("_").append(std::to_string(n2_index + 1)).append("_").append(std::to_string(n3_index));
        bead_indices.insert(index);

        // (1,0,0)
        index = std::to_string(n1_index + 1).append("_").append(std::to_string(n2_index)).append("_").append(std::to_string(n3_index));
        bead_indices.insert(index);

        // (0,-1,0)
        index = std::to_string(n1_index).append("_").append(std::to_string(n2_index-1)).append("_").append(std::to_string(n3_index));
        bead_indices.insert(index);

        // (-1,-1,0)
        index = std::to_string(n1_index-1).append("_").append(std::to_string(n2_index-1)).append("_").append(std::to_string(n3_index));
        bead_indices.insert(index);

        // (-1,0,0)
        index = std::to_string(n1_index-1).append("_").append(std::to_string(n2_index)).append("_").append(std::to_string(n3_index));
        bead_indices.insert(index);

        // (0,0,0) - basis
        vector3 temp = (a1 * n1_index + a2 * n2_index + a3 * n3_index) - basis2;
        n3_index_b2 = (int)std::round( (3*temp.z - sqrt6r)*inv2sqrt6r );
        n2_index_b2 = (int)std::round( invsqrt3r * (sqrt3*temp.x + temp.y - sqrt3inv3r) );
        n1_index_b2 = (int)std::round( invsqrt3r * (sqrt3*temp.x - temp.y - 2*sqrt3inv3r) );

        // convert this to the basis indices
        //index = std::to_string(n1_index).append("_").append(std::to_string(n2_index)).append("_").append(std::to_string(n3_index)).append("_-basis");
        index = std::to_string(n1_index_b2).append("_").append(std::to_string(n2_index_b2)).append("_").append(std::to_string(n3_index_b2)).append("_basis");
        bead_indices.insert(index);


        // (1,1,0) - basis
        temp = (a1 * (n1_index+1) + a2 * (n2_index+1) + a3 * n3_index) - basis2;
        n3_index_b2 = (int)std::round( (3*temp.z - sqrt6r)*inv2sqrt6r );
        n2_index_b2 = (int)std::round( invsqrt3r * (sqrt3*temp.x + temp.y - sqrt3inv3r) );
        n1_index_b2 = (int)std::round( invsqrt3r * (sqrt3*temp.x - temp.y - 2*sqrt3inv3r) );

        //index = std::to_string(n1_index+1).append("_").append(std::to_string(n2_index+1)).append("_").append(std::to_string(n3_index)).append("_-basis");
        index = std::to_string(n1_index_b2).append("_").append(std::to_string(n2_index_b2)).append("_").append(std::to_string(n3_index)).append("_basis");
        bead_indices.insert(index);

        // (1,0,0) - basis
        temp = (a1 * (n1_index+1) + a2 * (n2_index) + a3 * n3_index) - basis2;
        n3_index_b2 = (int)std::round( (3*temp.z - sqrt6r)*inv2sqrt6r );
        n2_index_b2 = (int)std::round( invsqrt3r * (sqrt3*temp.x + temp.y - sqrt3inv3r) );
        n1_index_b2 = (int)std::round( invsqrt3r * (sqrt3*temp.x - temp.y - 2*sqrt3inv3r) );

        index = std::to_string(n1_index_b2).append("_").append(std::to_string(n2_index_b2)).append("_").append(std::to_string(n3_index_b2)).append("_basis");
        bead_indices.insert(index);

        // (0,0,0) + basis
        temp = (a1 * (n1_index) + a2 * (n2_index) + a3 * n3_index) + basis2;
        n3_index_b2 = (int)std::round( (3*temp.z - sqrt6r)*inv2sqrt6r );
        n2_index_b2 = (int)std::round( invsqrt3r * (sqrt3*temp.x + temp.y - sqrt3inv3r) );
        n1_index_b2 = (int)std::round( invsqrt3r * (sqrt3*temp.x - temp.y - 2*sqrt3inv3r) );

        index = std::to_string(n1_index_b2).append("_").append(std::to_string(n2_index_b2)).append("_").append(std::to_string(n3_index_b2)).append("_basis");
        bead_indices.insert(index);

        // (-1,0,0) + basis
        temp = (a1 * (n1_index-1) + a2 * (n2_index) + a3 * n3_index) + basis2;
        n3_index_b2 = (int)std::round( (3*temp.z - sqrt6r)*inv2sqrt6r );
        n2_index_b2 = (int)std::round( invsqrt3r * (sqrt3*temp.x + temp.y - sqrt3inv3r) );
        n1_index_b2 = (int)std::round( invsqrt3r * (sqrt3*temp.x - temp.y - 2*sqrt3inv3r) );

        index = std::to_string(n1_index_b2).append("_").append(std::to_string(n2_index_b2)).append("_").append(std::to_string(n3_index_b2)).append("_basis");
        bead_indices.insert(index);

        // (-1,-1,0) + basis
        temp = (a1 * (n1_index-1) + a2 * (n2_index-1) + a3 * n3_index) + basis2;
        n3_index_b2 = (int)std::round( (3*temp.z - sqrt6r)*inv2sqrt6r );
        n2_index_b2 = (int)std::round( invsqrt3r * (sqrt3*temp.x + temp.y - sqrt3inv3r) );
        n1_index_b2 = (int)std::round( invsqrt3r * (sqrt3*temp.x - temp.y - 2*sqrt3inv3r) );

        index = std::to_string(n1_index_b2).append("_").append(std::to_string(n2_index_b2)).append("_").append(std::to_string(n3_index_b2)).append("_basis");
        bead_indices.insert(index);
    }

    std::vector<std::string> tempLine;
    boost::regex ifPlusB("\\basis"); // match any character

    for (auto & ind : bead_indices){
        tempLine.clear();

        boost::split(tempLine, ind, boost::is_any_of("_"), boost::token_compress_on);

        n1_index = (int)std::stoi(tempLine[0].c_str());
        n2_index = (int)std::stoi(tempLine[1].c_str());
        n3_index = (int)std::stoi(tempLine[2].c_str());

        vector3 temp = (a1 * n1_index + a2 * n2_index + a3 * n3_index);

        if (tempLine.size() > 3){ // determine if need to add or subtract the basis
            temp = temp + basis2;
        }

        beads.emplace_back(Bead(temp.x, temp.y, temp.z, 1));
    }

//    float smallest = FLT_MAX;
//    for(unsigned int i=0; i<beads.size(); i++){
//        const vector3 temp = beads[i].getVec();
//
//        for(unsigned int j=i+1; j<beads.size(); j++){
//            float len = (temp -  beads[j].getVec()).length();
//            if (len < smallest){
//                smallest = len;
//            }
//        }
//    }
//    std::cout << " Smallest " << smallest << std::endl;

    bead_phis.resize(beads.size());
    bead_thetas.resize(beads.size());
    bead_rvalues.resize(beads.size());
    float rval;

    unsigned int ind = 0;

    for (auto & bead : beads){
        rval = bead.getVec().length();
        if (rval > 0){
            bead_rvalues[ind] = rval;
            bead_thetas[ind] = acosf(bead.getZ()/rval);
            bead_phis[ind] = (bead.getX() == 0 || bead.getX() == 0.0f) ? 0.0f : atan2f(bead.getY(), bead.getX()); // need to check if this is a good way to test for 0
        } else {
            bead_rvalues[ind] = 0;
            bead_thetas[ind] = 0;
            bead_phis[ind] = 0;
        }

        ind++;
    }

    // write out bead model
    std::string nameOf = "beads.pdb";
    //const char * outputFileName = nameOf.c_str();
    FILE * pFile;
    pFile = fopen(nameOf.c_str(), "w");

    unsigned int residue_index=1;
    std::vector<std::string> chains(25);
    chains[0] = "A";
    chains[1] = "B";
    chains[2] = "C";
    chains[3] = "D";
    chains[4] = "E";
    chains[5] = "F";
    chains[6] = "G";
    chains[7] = "H";
    chains[8] = "I";
    chains[9] = "K";
    chains[10] = "L";
    chains[11] = "M";
    chains[12] = "N";
    chains[13] = "O";
    chains[14] = "P";
    chains[15] = "Q";
    chains[16] = "R";
    chains[17] = "S";
    chains[18] = "T";
    chains[19] = "U";
    chains[20] = "V";
    chains[21] = "W";
    chains[22] = "X";
    chains[23] = "Y";
    chains[24] = "J";

    std::string chainID = chains[0];

    //fprintf(pFile, "REMARK BEAD RADIUS %.3f\nREMARK CONTACTS PER BEAD %i\nREMARK VOLUME UPPER: %i LOWER: %i\n", this->bead_radius );
    //  fprintf(pFile, "REMARK TEMP RANGE %.3f => %.4E\nREMARK TOTAL TEMPERATURE STEPS %i\nREMARK TOTAL BEADS %i\n", "ATOM", i+1, "CA", "ALA", "A", residue_index.c_str(), currentBead->getX(), currentBead->getY(), currentBead->getZ() );
    unsigned int residueCnt = 1;
    unsigned locale = 0;

    for (auto & b : beads){

        if (residueCnt%9999 == 0){ // reinitialize counter
            residueCnt = 1;
            residue_index = 1;
            locale += 1;
            if (locale > 0){
                chainID = chains[locale];
            }
        }

        if (b.getContrast() > 0){
            printAtomLine(pFile, residueCnt, chainID, residue_index, b.getX(), b.getY(), b.getZ() );
        } else {
            printAtomLine(pFile, residueCnt, chains[locale+1], residue_index, b.getX(), b.getY(), b.getZ() );
        }

        residue_index += 1;
        residueCnt++;
    }

    fclose(pFile);

    // calculate scattering from the beads as my Solvent mask
}

/*
 * Use SIMDe to performed vector sum
 */
float AtomisticModel::addUp(const unsigned int total_width, float * vector_to_sum) {

    unsigned int window = 8; // window is 8 for _m128 and 16 for _m256

    const unsigned int alignedN = total_width - total_width % window;
    simde__m128 mmSum1, mmSum2, v0, v1;

    auto * sum_sse = (float *)_aligned_malloc(sizeof(float)*4, 16);

    for (int i = 0; i < alignedN; i+=window) {
        v0 = simde_mm_load_ps(&vector_to_sum[i]);
        mmSum1 = simde_mm_add_ps(mmSum1, v0);

        v1 = simde_mm_load_ps(&vector_to_sum[i+4]);
        mmSum1 = simde_mm_add_ps(mmSum1, v1);
    }

    simde_mm_store_ps(sum_sse, mmSum1);

    float sum = sum_sse[0] + sum_sse[1] + sum_sse[2] + sum_sse[3];

    for (int i = alignedN; i < total_width; ++i) {
        sum += vector_to_sum[i];
    }

    return sum;
}

/*
 * this can be distributed, treating each q-value separately
 *
 */
float AtomisticModel::excluded_volume_scattering(float probe_radius, float q_value) {

    const float * pXv = pdbModel.getX();
    const float * pYv = pdbModel.getY();
    const float * pZv = pdbModel.getZ();
    std::vector<std::vector<int>> cells;

    std::vector<vector3> pdb_vectors;

    // create envelope model
    vector3 low_vec = pdbModel.getLowerBoundVector() - probe_radius;
    vector3 high_vec = pdbModel.getUpperBoundVector() + probe_radius;

    std::cout << pdbModel.getLowerBoundVector().length() << std::endl;

    vector3 diff_vec = (high_vec - low_vec)/probe_radius;

    int nx = static_cast<int>(std::ceil(diff_vec.x)) + 1;
    int ny = static_cast<int>(std::ceil(diff_vec.y)) + 1;
    int nz = static_cast<int>(std::ceil(diff_vec.z)) + 1;

    cells.assign(static_cast<std::size_t>(nx)*ny*nz, {});

    for(unsigned int i=0; i<totalAtoms; i++) {
            pdb_vectors.emplace_back(pXv[i], pYv[i], pZv[i]);
            cells[flat(nx, ny, nz, probe_radius, &pXv[i], &pYv[i], &pZv[i], low_vec)].push_back(static_cast<int>(i));
    }

    float limit = probe_radius*probe_radius;
    unsigned int total_iterations = 100000;
    float integrated_radius = 14.0f;
    float integrated_volume = 4.0f/3.0f*M_PI*integrated_radius*integrated_radius*integrated_radius;
    float mc_normalization = integrated_volume/(float)total_iterations;

    // integrate from 0 to 14
    float radial_dis;

    std::mt19937 rng(static_cast<unsigned>(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
    // Uniform [0,1] for the radius transform
    std::uniform_real_distribution<float> uni(0.0, 1.0);
    // Standard normal for direction
    std::normal_distribution<float> normal(0.0, 1.0);

    float integrand = 0.0;
    // this can be distributed as each atom can be treated separately and summed
    for(unsigned int i=0; i<totalAtoms; i++) {
         vector3 vec(pXv[i], pYv[i], pZv[i]);

        std::vector<vector3> random_waters;                                     // one per thread, declared outside the loop
        water_shell::randomWaterShell(vec, rng, random_waters, integrated_radius);
        float sum = 0.0f;

        for (unsigned int j=0; j < random_waters.size(); j++) {

            vector3 * point = &random_waters[j];
            // check if point is within envelope
            bool found = false;
            vector3 diff = *point - vec;
            vector3 pp = (*point-low_vec)/probe_radius;

            // does the point overlap with the structure?
            forNear(pp, nx, ny, nz, cells, [&](int idx) {
                if (found) return;
                if ((*point - pdb_vectors[idx]).sqlength() <= limit) found = true;
            });


            if (found) { //evaluate structure factor
                radial_dis = diff.length();
                if (radial_dis >= 2.5) {
                    sum += radial_dis*radial_dis*( g_OO_SQ_matched(radial_dis) - 1.0f ) * sinc(radial_dis*q_value);
                } else {
                    sum += radial_dis*radial_dis*( - 1.0f ) * sinc(radial_dis*q_value);
                }
                //std::printf("%-6s%5i %4s %3s %1s%4s    %8.3f%8.3f%8.3f  1.00100.00\n", "ATOM", j," CA ", "ALA", "A", "1", point->x, point->y, point->z );
            }

            //std::printf("%-6s%5i %4s %3s %1s%4s    %8.3f%8.3f%8.3f  1.00100.00\n", "ATOM", j," CA ", "ALA", "A", "1", point->x, point->y, point->z );
        }

        // do the monte carlo
        // perform the integration for each point in excluded volume
        /*
         * peak 1 2.75 to 2.8
         * peak 2 3.65 to 4.5
         * peak 3 6.7 to 6.87
         * peak 4 8.85
         */

        // for (unsigned int j=0; j < total_iterations; j++) {
        //     float radial_dis = 0;
        //     vector3 point = randomPointInSphere(vec, integrated_radius, rng, uni, normal, radial_dis);
        //     // check if point is within envelope
        //     bool found = false;
        //     vector3 pp = (point-low_vec)/probe_radius;
        //
        //     forNear(pp, nx, ny, nz, cells, [&](int idx) {
        //         if (found) return;
        //         if ((point - pdb_vectors[idx]).sqlength() <= limit) found = true;
        //     });
        //
        //     if (found) { //evaluate structure factor
        //        if (radial_dis >= 2.5) {
        //            sum += radial_dis*radial_dis*( g_OO_SQ_matched(radial_dis) - 1.0f ) * sinc(radial_dis*q_value);
        //        } else {
        //            sum += radial_dis*radial_dis*( - 1.0f ) * sinc(radial_dis*q_value);
        //        }
        //         //std::printf("%-6s%5i %4s %3s %1s%4s    %8.3f%8.3f%8.3f  1.00100.00\n", "ATOM", j," CA ", "ALA", "A", "1", point.x, point.y, point.z );
        //     }
        // }

        integrand += sum;
        std::cout << " integrand: " << integrand << std::endl;
        exit(0);
    }


    return mc_normalization*integrand;

}
