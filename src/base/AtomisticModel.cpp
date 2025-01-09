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

    totalDistances = totalAtoms*(totalAtoms-1)/2;
    binned_distances.resize(totalDistances);

    const float * pXv = pdbModel.getCenteredXVec();
    const float * pYv = pdbModel.getCenteredYVec();
    const float * pZv = pdbModel.getCenteredZVec();

    const int * pAtomicNumbers = pdbModel.getAtomicNumberVec();

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

    atomList.insert(99);  // 99 is unique for water

    hydrogens = {
            {"ALA", 7},
            {"ASN", 8},
            {"ASP", 7},
            {"CYS", 7},
            {"DA", 14},
            {"DG", 11},
            {"DC", 14},
            {"DU", 13},
            {"GLN", 10},
            {"GLU", 9},
            {"GLY", 5},
            {"HIS", 11},
            {"ILE", 24},
            {"LEU", 13},
            {"LYS", 14},
            {"MET", 11},
            {"PHE", 11},
            {"PRO", 9},
            {"SER", 7},
            {"THR", 9},
            {"TRP", 12},
            {"TYR", 11},
            {"VAL", 11},
            {"rA", 14},
            {"rG", 11},
            {"rC", 14},
            {"rU", 13},
            {"HOH", 2},
            {"AMP", 14},
            {"GMP", 11}
    };
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

    this->lmax = lmax;
    unsigned int bessel_index, q_index=0, lmax_plus_1 = lmax + 1;


    if (recalculate){

        sbj.recalculate(totalqvalues, qvalues, &rvalues);
        //sbj_x.recalculate(totalqvalues, qvalues, &bead_rvalues);

    } else {

        createHCPGridModel(this->getDmax()/(float)lmax_plus_1);
        // SHE depend on atomic coordinates (r, theta, phi)
        this->createSphericalHarmonics(lmax);
        //this->createSphericalHarmonicsEx();

        // SBJ depends on qvalues and r
        this->createSphericalBessels(totalqvalues, qvalues);
        //this->createSphericalBesselsEx(totalqvalues, qvalues);
    }


    if (neighbors.empty()){
        this->createNeighborhoods(1.5*1.5);
    }

    // for each q-value calculate atomic scattering form factor
    // need to assemble into partial sums at each q-value?

    std::map<int, float> asf; // for each q-value, this is calculated over the number of unique atoms

    for(auto & atom : atomList){
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

    float *pYLM, bessel_product, bessel_product_Ex, waterASF, * pBessel;

    float partialSumRAp, partialSumIAp, partialSumREx, partialSumIEx;

    auto pSBJ = sbj.getpSphericalBessels();
    auto pAtomicNumbers = pdbModel.getAtomicNumberVec();

    this->setFractionalVolume(29.9); // set fractional volume of the waters in excluded volume
    this->setGaussianRadii();

    auto pReal = she.getpDataYLMReal();
    auto pImag = she.getpDataYLMImag();

    std::clock_t startTime = std::clock();

    unsigned int lmax_plus_1_totalAtoms = lmax_plus_1*totalAtoms;
    unsigned int long ylm_index, base_index;

    for(auto & qvalue : qvalues){

        for(auto & ff : asf){ // calculate atomic scattering form factors for given qvalue
            ff.second = functions_WETSAXS::asf(ff.first, qvalue);
        }

        waterASF = functions_WETSAXS::asf(99, qvalue);
        this->setFractionalVolumeAtQ(qvalue);

        // populate qr values
        unsigned int ll_index = (lmax+1)*(lmax+1)*q_index;// + l*l + l + m;

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
            partialSumIAp = 0.0f;

            // excluded volume term based on atomic coordinates
            partialSumREx = 0.0f;
            partialSumIEx = 0.0f;

            ylm_index = base_index + l * totalAtoms; // for m = 0

            for (unsigned int n = 0; n < totalAtoms; n++) { // over all atoms
                pBessel = &pSBJ[bessel_index + n];

                bessel_product = asf[pAtomicNumbers[n]] * *pBessel; // for each m value, using same bessel values at l
                bessel_product_Ex = (waterASF - fractionOccAtoms[n]) * *pBessel; // for estimating excluded volume using water ASF

                // real term
                pYLM = &pReal[ylm_index];
                partialSumRAp += bessel_product * *pYLM;     // particle, real
                partialSumREx += bessel_product_Ex * *pYLM;  // excluded volume, real

                //imaginary term
                pYLM = &pImag[ylm_index];
                partialSumIAp += bessel_product * *pYLM;    // particle, imaginary
                partialSumIEx += bessel_product_Ex * *pYLM; // excluded volume, imaginary

                ++ylm_index;
            }

            int partials_index = ll_index + l*l + l;

            // particle scattering for m=0
            ptrRAlm[partials_index] = partialSumRAp;
            ptrIAlm[partials_index] = partialSumIAp;

            // excluded volume scattering for m=0
            ptrREx[partials_index] = partialSumREx; // uses volume fraction with assigned waters
            ptrIEx[partials_index] = partialSumIEx; // final should be negative of the amplitude (Babinet's)

            // biggest time saving will be here, need to exploit the symmetry of spherical hamormonics (for m = -m)
            for(int m=1; m<=l; m++) {

                // for given (l,m) calculate Spherical Harmonic
                partialSumRAp = 0.0f;
                partialSumIAp = 0.0f;

                // excluded volume term based on atomic coordinates
                partialSumREx = 0.0f;
                partialSumIEx = 0.0f;

                ylm_index = base_index + (l + m) * totalAtoms; // incorrect

                // for each (l,m) calculate partial sum at given q value
                for (unsigned int n = 0; n < totalAtoms; n++) { // over all atoms

                    pBessel = &pSBJ[bessel_index + n];
                    bessel_product = asf[pAtomicNumbers[n]] * *pBessel; // for each m value, using same bessel values at l
                    bessel_product_Ex = (waterASF - fractionOccAtoms[n]) * *pBessel; // for estimating excluded volume using water ASF

                    // real term
                    pYLM = &pReal[ylm_index];
                    partialSumRAp += bessel_product * *pYLM;     // particle, real
                    partialSumREx += bessel_product_Ex * *pYLM;  // excluded volume, real

                    //imaginary term
                    pYLM = &pImag[ylm_index];
                    partialSumIAp += bessel_product * *pYLM;    // particle, imaginary
                    partialSumIEx += bessel_product_Ex * *pYLM; // excluded volume, imaginary

                    ++ylm_index;
                }

                int update_index = partials_index + m;

                // particle scattering
                ptrRAlm[update_index] = partialSumRAp;
                ptrIAlm[update_index] = partialSumIAp;

                // excluded volume
                ptrREx[update_index] = partialSumREx; // uses volume fraction with assigned waters
                ptrIEx[update_index] = partialSumIEx; // final should be negative of the amplitude (Babinet's)

                // partials_index = ll_index + l*l + l;
                update_index = partials_index - m;
                if (m & 1) {// if m is odd, flip sign of the real part
                    // std::cout << count << " ( " << (partials_index -m) <<  " ) "  << ptrRAlm[partials_index -m] << " " << (-1*partialSumRAp) << std::endl;
                    ptrRAlm[update_index] = -partialSumRAp;
                    ptrIAlm[update_index] = partialSumIAp;

                    ptrREx[update_index] = -partialSumREx; // uses volume fraction with assigned waters
                    ptrIEx[update_index] = partialSumIEx;  // final should be negative of the amplitude (Babinet's)

                } else { // if m is even, flip sign of the Imaginary part
                    // std::cout << count << " ( " << (partials_index -m) <<  " ) "  << ptrIAlm[partials_index -m] << " " << (-1*partialSumIAp) << std::endl;
                    ptrRAlm[update_index] = partialSumRAp;
                    ptrIAlm[update_index] = -partialSumIAp;

                    ptrREx[update_index] = partialSumREx;  // uses volume fraction with assigned waters
                    ptrIEx[update_index] = -partialSumIEx; // final should be negative of the amplitude (Babinet's)
                }
            }
        }
        ++q_index;
    }

    //calculateExVolPartialAmplitudes(lmax, totalqvalues, qvalues, recalculate);
    std::cout << " Finished Partials - starting cross " << std::endl;

//    float delta_r = pdbModel.getDmax()/(15*2);
//    std::vector<float> sinqr_qr(15*2);
//    std::vector<float> r_values(15*2);
//
//    for(int i=0; i < sinqr_qr.size(); i++){
//        r_values[i] = delta_r*i + 0.5*delta_r;
//    }
//
//    pdbModel.setVectorModelToCenter();
//    auto * pCoords = pdbModel.getModel();
//
//    std::vector<float> new_qvalues(800);
//    float deltaq = (0.55 - 0.008)/800;
//    float incr = 0.008;
//    for(auto & qq : new_qvalues){
//        qq = incr;
//        incr += deltaq;
//    }
//    for(auto & qvalue : new_qvalues){
//
//        waterASF = functions_WETSAXS::asf(99, qvalue);
//        waterASF *= waterASF;
//
//        for(int i=0; i < sinqr_qr.size(); i++){
//            float qr_val = r_values[i]*qvalue;
//            sinqr_qr[i] = sin(qr_val)/qr_val;
//        }
//
//        float f_i_f_i_product = 0;
//        float neighbors = 0;
//
//        for (unsigned int n = 0; n < totalAtoms; n++) { // over all atoms
//
//            int neighborhood = n*totalAtoms - n*(n+1)/2 - n; // + col -1;
//            const vector3 * n_vec = &pCoords[n];
//
//            float * val = &fractionOccAtoms[n];
//
//            f_i_f_i_product += 1;// * *val * *val; // sinqr_qr for r = 0;
//
//            for (unsigned int m = n+1; m < totalAtoms; m++) { // over all atoms
//                // check if atoms are neighbors
//                if (nearestNeighbors.find((neighborhood + m) - 1) == nearestNeighbors.end()){ // non-neighbors
//                    //f_i_f_i_product += *val*fractionOccAtoms[m]*2*sinqr_qr[(int)std::floor((*n_vec - pCoords[m]).length()/delta_r)];
//                    f_i_f_i_product += 2*sinqr_qr[(int)std::floor((*n_vec - pCoords[m]).length()/delta_r)];
//                } else {
//                    //neighbors += *val*fractionOccAtoms[m]*2*sinqr_qr[(int)std::floor((*n_vec - pCoords[m]).length()/delta_r)];
//                    neighbors += sinqr_qr[(int)std::floor((*n_vec - pCoords[m]).length()/delta_r)];
//                }
//            }
//        }
//
//        float b_factor = expf(-0.5*qvalue*qvalue*1.2);
//        float i_of_q = waterASF*(f_i_f_i_product + b_factor*2*neighbors);
//        std::cout << q_index << " " << qvalue << " " << i_of_q << std::endl;
//    }

//    clmNeighbors.resize(qvalues.size());
//    clmNonNeighbors.resize(qvalues.size());
//
//    float * bessel_n, * bessel_k;
//    float * n_plm_real, * n_plm_imag;
//    unsigned int long ylm_index_k;
//    float bessel_frac_n, bessel_frac_k;
//    float * const ptrExNeighbors = (totalAtoms != 0) ? clmNeighbors.data() : nullptr;
//    float * const ptrExNonNeighbors = (totalAtoms != 0) ? clmNonNeighbors.data() : nullptr;
//
//    q_index = 0;
//    for(auto & qvalue : qvalues){
//
//        waterASF = functions_WETSAXS::asf(99, qvalue);
//        waterASF *= waterASF;
//
//        float i_q_non = 0, i_q_nn = 0;
//
//        for(int l=0; l<=lmax; l++){
//
//            bessel_index = lmax_plus_1_totalAtoms*q_index + totalAtoms*l;
//            base_index = (l*l)*totalAtoms; // always positive
//
//            // calculate spherical bessels for given q, l over all atoms
//
//            /*
//             * for m = 0
//             */
//            ylm_index = base_index + l * totalAtoms; // for m = 0
//            ylm_index_k = base_index + l * totalAtoms;
//
//            float f_i_f_i_product = 0, cross_term_real, cross_term_imag;
//
//            for(auto const & non_n : non_neighbors){ // iterates over all atoms
//
//                //bessel_n = &pSBJ[bessel_index + non_n.first]; // every atom has a bessel value
//                //bessel_frac_n = pSBJ[bessel_index + non_n.first] * fractionOccAtoms[non_n.first];
//                bessel_frac_n = pSBJ[bessel_index + non_n.first];
//                // real and imag terms
//                n_plm_real = &pReal[ylm_index]; // j_l * P_lm * cos(phi)
//                n_plm_imag = &pImag[ylm_index]; // j_l * P_lm * sin(phi)
//
//                f_i_f_i_product += bessel_frac_n * bessel_frac_n * (*n_plm_real * *n_plm_real + *n_plm_imag * *n_plm_imag);
//                cross_term_real = 0;
//                cross_term_imag = 0;
//
//                for (auto const & k : non_n.second){ // cross-terms, 2x because of i,j and j,i are identical
//                    //f_i_f_i_product += 2*bessel_n*pSBJ[bessel_index + k]*(n_plm_real*pReal[ylm_index_k + k] + n_plm_imag*pImag[ylm_index_k + k]);
//                    //bessel_frac_k = pSBJ[bessel_index + k] * fractionOccAtoms[k];
//                    bessel_frac_k = pSBJ[bessel_index + k];
//                    cross_term_real += bessel_frac_k*pReal[ylm_index_k + k];
//                    cross_term_imag += bessel_frac_k*pImag[ylm_index_k + k];
//                }
//
//                f_i_f_i_product += 2 * bessel_frac_n *(*n_plm_real*cross_term_real + *n_plm_imag*cross_term_imag);
//                ++ylm_index;
//            }
//
//            i_q_non += f_i_f_i_product; //add the m=0 term
//            // nearest neighbors
//            float nneighbors = 0;
//            ylm_index = base_index + l * totalAtoms;
//            ylm_index_k = base_index + l * totalAtoms;
//
//            for(auto const & n_n : neighbors){ // iterates over all atoms
//
//                // real and imag terms
//                cross_term_real = 0;
//                cross_term_imag = 0;
//
//                for (auto const & k : n_n.second){ // cross-terms
//                   // bessel_k = &pSBJ[bessel_index + k];
//                    //bessel_frac_k = pSBJ[bessel_index + k] * fractionOccAtoms[k];
//                    bessel_frac_k = pSBJ[bessel_index + k];
//                    cross_term_real += bessel_frac_k*pReal[ylm_index_k + k];
//                    cross_term_imag += bessel_frac_k*pImag[ylm_index_k + k];
//                }
//                nneighbors += pSBJ[bessel_index + n_n.first] * (pReal[ylm_index] * cross_term_real + pImag[ylm_index] * cross_term_imag);
//                //nneighbors += pSBJ[bessel_index + n_n.first] * fractionOccAtoms[n_n.first] *(pReal[ylm_index] * cross_term_real + pImag[ylm_index] * cross_term_imag);
//
//                ++ylm_index;
//            }
//
//            i_q_nn += 2*nneighbors; // add the m=0 term
//
//            //int partials_index = ll_index + l*l + l;
//            // excluded volume scattering for m=0 for given l
//            //ptrExNeighbors[partials_index] = waterASF * nneighbors; // uses volume fraction with assigned waters
//            //ptrExNonNeighbors[partials_index] = waterASF * f_i_f_i_product; // final should be negative of the amplitude (Babinet's)
//
//            // biggest time saving will be here, need to exploit the symmetry of spherical hamormonics (for m = -m)
//            for(int m=1; m<=l; m++) {
//                // for given (l,m) calculate Spherical Harmonic
//                // excluded volume term based on atomic coordinates
//                ylm_index = base_index + (l + m) * totalAtoms; //
//                ylm_index_k = base_index + (l + m) * totalAtoms;
//
//                f_i_f_i_product = 0;
//
//                // for each (l,m) calculate partial sum at given q value
//                for(auto const & non_n : non_neighbors){ // over all atoms via map
//
//                    //bessel_n = &pSBJ[bessel_index + non_n.first]; // for estimating excluded volume using water ASF
//                    //bessel_frac_n = pSBJ[bessel_index + non_n.first] * fractionOccAtoms[non_n.first];
//                    bessel_frac_n = pSBJ[bessel_index + non_n.first];
//                    // real and imag terms
//                    n_plm_real = &pReal[ylm_index]; //  P_lm * cos(phi)
//                    n_plm_imag = &pImag[ylm_index]; //  P_lm * sin(phi)
//
//                    f_i_f_i_product += bessel_frac_n * bessel_frac_n * (*n_plm_real * *n_plm_real + *n_plm_imag * *n_plm_imag);
//                    cross_term_real = 0;
//                    cross_term_imag = 0;
//
//                    for (auto const & k : non_n.second){ // cross-terms
//                        //bessel_k = &pSBJ[bessel_index + k];
//                        //bessel_frac_k = pSBJ[bessel_index + k] * fractionOccAtoms[k];
//                        bessel_frac_k = pSBJ[bessel_index + k];
//                        cross_term_real += bessel_frac_k*pReal[ylm_index_k + k];
//                        cross_term_imag += bessel_frac_k*pImag[ylm_index_k + k];
//                    }
//
//                    f_i_f_i_product += 2 * bessel_frac_n*(*n_plm_real * cross_term_real + *n_plm_imag * cross_term_imag);
//                    ++ylm_index;
//                }
//
//                // nearest neighbors (nneighbors)
//                nneighbors = 0;
//                ylm_index = base_index + (l + m) * totalAtoms;
//
//                for(auto const & n_n : neighbors){ // iterates over all atoms
//                    // real and imag terms
//                    cross_term_real = 0;
//                    cross_term_imag = 0;
//                    for (auto const & k : n_n.second){ // cross-terms
//                        //bessel_k = &pSBJ[bessel_index + k];
//                        //bessel_frac_k = pSBJ[bessel_index + k] * fractionOccAtoms[k];
//                        bessel_frac_k = pSBJ[bessel_index + k];
//                        cross_term_real += bessel_frac_k*pReal[ylm_index_k + k];
//                        cross_term_imag += bessel_frac_k*pImag[ylm_index_k + k];
//                    }
//                    nneighbors += pSBJ[bessel_index + n_n.first] * (pReal[ylm_index] * cross_term_real + pImag[ylm_index] * cross_term_imag);
//                    //nneighbors += pSBJ[bessel_index + n_n.first] * fractionOccAtoms[n_n.first] * (pReal[ylm_index] * cross_term_real + pImag[ylm_index] * cross_term_imag);
//                    ++ylm_index;
//                }
//                nneighbors *= 2; // multiply by 2 for (j,k) and (k,j)
//
//                i_q_nn += 2*nneighbors; // multiply by 2x for m = -1 and m = 1;
//                i_q_non += 2*f_i_f_i_product; // multiply by 2x for m = -1 and m = 1;
//            }
//        }
//
//        ptrExNeighbors[q_index] = waterASF*i_q_nn;
//        ptrExNonNeighbors[q_index] = waterASF*i_q_non;
//
//        // std::cout << q_index << " " << qvalue << " " << ptrExNeighbors[q_index] + ptrExNonNeighbors[q_index] << std::endl;
//        ++q_index;
//    }

    // calculate I(0) for protein and excluded volume

    try {
        setMatchPointScale(pAtomicNumbers, fractionOccAtoms.data());
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
 * Unless the model includes hydrogens, it will only consider non-hydrogen atoms and will under-estimate forward scattering
 * average electron density of a protein should be around 0.4 electrons per cubic Angstrom
 *
 */
void AtomisticModel::setMatchPointScale(int * pAtomicNumbers, float * fractionalOccupancies){

    forwardScatteringParticle = 0.0f;
    forwardScatteringExVolume = 0.0f;

    // water_asf at zero will be 10 electrons
    const float water_asf = functions_WETSAXS::asf_at_q_zero(99);
    float totalH2Os= 0;
    float diff;

    for(unsigned int n=0; n < totalAtoms; n++){ // over all atoms
        diff = functions_WETSAXS::asf_at_q_zero(pAtomicNumbers[n]);
        forwardScatteringParticle += diff*diff;
        totalH2Os += fractionalOccupancies[n];
    }

    float pw_vol = 29.9*0.334; // 10 e_n
    float factor;

    for(unsigned int i=0; i<totalAtoms; i++){
        diff = gaussian_radii[i]/1.75;
        factor = pw_vol*(1-diff*diff*diff);
        diff = (water_asf - factor);
        forwardScatteringExVolume += diff*diff; // normalized to volume of water, looking at fractions of water volume
    }

//    forwardScatteringExVolume =  beads.size()*(water_asf - pw_vol*(1 - 0.6366*0.6366*0.6366/(1.75*1.75*1.75)));
//    forwardScatteringExVolume *= forwardScatteringExVolume;

    forwardScatteringParticle = std::sqrtf(forwardScatteringParticle);
    forwardScatteringExVolume = std::sqrtf(forwardScatteringExVolume);

    totalH2Os = std::ceil(totalH2Os);

    if (totalAtoms < 1 || forwardScatteringExVolume == 0)
        throw std::overflow_error("Divide by zero exception");

    SASTOOLS_UTILS_H::logger("density", formatNumber(forwardScatteringParticle/pdbModel.getVolume(), 3));
    SASTOOLS_UTILS_H::logger("Total H2Os in Ex Volume", formatNumber(totalH2Os, 0));

    // proteins have a match point of 0.42 e_n per Angstrom3
    // nucleics have a match point of 0.55 e_n per Angstrom3
    matchPoint = forwardScatteringParticle/forwardScatteringExVolume;
}


/*
 * Unless the model includes hydrogens, it will only consider non-hydrogen atoms and will under-estimate forward scattering
 */
void AtomisticModel::setMatchPointScale(unsigned int totalHydrationWaters){

    float water_asf = functions_WETSAXS::asf_at_q_zero(99)*totalHydrationWaters;

    if (totalAtoms < 1 || forwardScatteringExVolume == 0)
        throw std::overflow_error("Divide by zero exception");

    matchPoint = std::sqrtf((forwardScatteringParticle + water_asf)*(forwardScatteringParticle + water_asf)/( (forwardScatteringExVolume+water_asf)*(forwardScatteringExVolume+water_asf))); // divide by zero possible
//    matchPoint = std::sqrtf((forwardScatteringParticle + water_asf)*(forwardScatteringParticle + water_asf)/(forwardScatteringExVolume*forwardScatteringExVolume)); // divide by zero possible
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
    auto pAtomicNumbers = pdbModel.getAtomicNumberVec();

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

void AtomisticModel::populateBinnedDistances(unsigned int lmax){

    const float * pX = pdbModel.getX();
    const float * pY = pdbModel.getY();
    const float * pZ = pdbModel.getZ();

    float inv_delta = (float)lmax/ pdbModel.getDmax();

    unsigned int counter = 0;
    unsigned int num;

    for(unsigned int i=0; i<totalAtoms; i++){
        vector3 vec1(pX[i], pY[i], pZ[i]);

        for(unsigned int j=(i+1); j<totalAtoms; j++){

            num = (unsigned int) std::floor ((vector3(pX[j], pY[j], pZ[j])-vec1).length()*inv_delta);

            if (num >= lmax){
                binned_distances[counter] = num - 1;
            } else {
                binned_distances[counter] = num;
            }

//            if (std::floor((vector3(pX[j], pY[j], pZ[j])-vec1).length()*inv_delta) >= lmax){
//                std::cout << "greater than :: " << std::floor((vector3(pX[j], pY[j], pZ[j])-vec1).length()*inv_delta) << " " << lmax << std::endl;
//            }
//            binned_distances[counter] = std::floor((vector3(pX[j], pY[j], pZ[j])-vec1).length()*inv_delta);

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

    const float * pXv = pdbModel.getCenteredXVec();
    const float * pYv = pdbModel.getCenteredYVec();
    const float * pZv = pdbModel.getCenteredZVec();

    if (!neighbors.empty()){
        neighbors.clear();
        non_neighbors.clear();
    }

    // calculate spherical coordinates of atoms - must be centered
    for(unsigned int row=0; row < totalAtoms; row++){

        const vector3 vec1 = vector3(pXv[row], pYv[row], pZv[row]);

        //const unsigned int row_index = row*totalAtoms - row*(row+1)/2 - row;
        //const unsigned int row_index = row*(totalAtoms - (row+1)/2 - 1);

        auto mit = neighbors.emplace(std::pair<unsigned int, std::vector<unsigned int>>(row, std::vector<unsigned int>()));
        auto non_it = non_neighbors.emplace(std::pair<unsigned int, std::vector<unsigned int>>(row, std::vector<unsigned int>()));

//        std::vector<unsigned int> neighborIndices;
//        std::vector<unsigned int> nonNeighborIndices;

        for(unsigned int col=(row+1); col < totalAtoms; col++){

            if ((vector3(pXv[col], pYv[col], pZv[col]) - vec1).sqlength() < cutoff){ // they are neighbors
                // calculate unique index value if neighbors using formula
                // nearestNeighbors.insert((row_index + col) - 1);
                //neighborIndices.push_back(col);
                mit.first->second.push_back(col);
            } else { // non neighbors
                non_it.first->second.push_back(col);
                //nonNeighborIndices.push_back(col);
            }
        }

//        neighbors.insert(std::pair<unsigned int, std::vector<unsigned int>>(row, std::vector<unsigned int>(neighborIndices)));
//        non_neighbors.insert(std::pair<unsigned int, std::vector<unsigned int>>(row, std::vector<unsigned int>(nonNeighborIndices)));
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

        if (tempLine.size() > 3){ // determine if need ot add or subtract the basis
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
