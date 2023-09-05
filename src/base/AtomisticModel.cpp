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

        rval = sqrtf(xval*xval + yval*yval +zval*zval);

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

    if (recalculate){

        sbj.recalculate(totalqvalues, qvalues, &rvalues);

    } else {
        // SHE depend on atomic coordinates (r, theta, phi)
        this->createSphericalHarmonics(lmax);

        // SBJ depends on qvalues and r
        this->createSphericalBessels(totalqvalues, qvalues);
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
    almRealExcludedVolSum.resize((lmax+1)*(lmax+1)*totalqvalues);; // partial sums
    almImagExcludedVolSum.resize((lmax+1)*(lmax+1)*totalqvalues);

    float * const ptrRAlm = (totalAtoms != 0) ? almRealSum.data() : nullptr;
    float * const ptrIAlm = (totalAtoms != 0) ? almImagSum.data() : nullptr;

    float * const ptrREx = (totalAtoms != 0) ? almRealExcludedVolSum.data() : nullptr;
    float * const ptrIEx = (totalAtoms != 0) ? almImagExcludedVolSum.data() : nullptr;

    float *pYLM, bessel_product, bessel_product_Ex, waterASF, vol, * pBessel;

    float partialSumRAp, partialSumIAp, partialSumREx, partialSumIEx;

    auto pSBJ = sbj.getpSphericalBessels();
    auto pAtomicNumbers = pdbModel.getAtomicNumberVec();
    auto pAtomicVolumes = pdbModel.getAtomicVolumeVec();

    this->setFractionalVolume(29.9); // set fractional volume of the waters in excluded volume

    auto pReal = she.getpDataYLMReal();
    auto pImag = she.getpDataYLMImag();

    std::clock_t startTime = std::clock();

    unsigned int bessel_index, q_index=0, lmax_plus_1 = lmax + 1;
    unsigned int lmax_plus_1_totalAtoms = lmax_plus_1*totalAtoms;
    unsigned int long ylm_index, base_index;

    for(auto & qvalue : qvalues){

        for(auto & ff : asf){ // calculate atomic scattering form factors for given qvalue
            ff.second = functions_WETSAXS::asf(ff.first, qvalue);
        }

        waterASF = functions_WETSAXS::asf(99, qvalue);

        // populate qr values
        int ll_index = (lmax+1)*(lmax+1)*q_index;// + l*l + l + m;

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
                bessel_product_Ex = fractionOccAtoms[n] * *pBessel; // for estimating excluded volume using water ASF

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
            ptrREx[partials_index] = waterASF * partialSumREx; // uses volume fraction with assigned waters
            ptrIEx[partials_index] = waterASF * partialSumIEx; // final should be negative of the amplitude (Babinet's)


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
                    bessel_product_Ex = fractionOccAtoms[n] * *pBessel; // for estimating excluded volume using water ASF

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
                partialSumREx *= waterASF;
                partialSumIEx *= waterASF;

                ptrREx[update_index] = partialSumREx; // uses volume fraction with assigned waters
                ptrIEx[update_index] = partialSumIEx; // final should be negative of the amplitude (Babinet's)

                update_index = partials_index - m;
                if (m & 1) {// if m is odd, flip sign of the real part
                    // std::cout << count << " ( " << (partials_index -m) <<  " ) "  << ptrRAlm[partials_index -m] << " " << (-1*partialSumRAp) << std::endl;
                    ptrRAlm[update_index] = -partialSumRAp;
                    ptrIAlm[update_index] = partialSumIAp;

                    ptrREx[update_index] = -partialSumREx; // uses volume fraction with assigned waters
                    ptrIEx[update_index] = partialSumIEx; // final should be negative of the amplitude (Babinet's)

                } else { // if m is even, flip sign of the Imaginary part
                    // std::cout << count << " ( " << (partials_index -m) <<  " ) "  << ptrIAlm[partials_index -m] << " " << (-1*partialSumIAp) << std::endl;
                    ptrRAlm[update_index] = partialSumRAp;
                    ptrIAlm[update_index] = -partialSumIAp;

                    ptrREx[update_index] = partialSumREx; // uses volume fraction with assigned waters
                    ptrIEx[update_index] = -partialSumIEx; // final should be negative of the amplitude (Babinet's)
                }

            }
        }
        ++q_index;
    }

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



/*
 * Unless the model includes hydrogens, it will only consider non-hydrogen atoms and will under-estimate forward scattering
 * average electron density of a protein should be around 0.4 electrons per cubic Angstrom
 *
 */
void AtomisticModel::setMatchPointScale(int * pAtomicNumbers, float * fractionalOccupancies){

    forwardScatteringParticle = 0.0f;
    forwardScatteringExVolume = 0.0f;

    // water_asf at zero will be 10 electrons
    float water_asf = functions_WETSAXS::asf_at_q_zero(99);

    for(unsigned int n=0;n<totalAtoms;n++){ // over all atoms
        forwardScatteringParticle += functions_WETSAXS::asf_at_q_zero(pAtomicNumbers[n]);
        forwardScatteringExVolume += water_asf*fractionalOccupancies[n]; // izero of water model mask
    }

    if (totalAtoms < 1 || forwardScatteringExVolume == 0)
        throw std::overflow_error("Divide by zero exception");


    SASTOOLS_UTILS_H::logger("density", formatNumber(forwardScatteringParticle/pdbModel.getVolume(), 3));


    // proteins have a match point of 0.42 e_n per Angstrom3
    // nucleics have a match point of 0.55 e_n per Angstrom3
    matchPoint = forwardScatteringParticle/forwardScatteringExVolume; // divide by zero possible
}


/*
 * Unless the model includes hydrogens, it will only consider non-hydrogen atoms and will under-estimate forward scattering
 */
void AtomisticModel::setMatchPointScale(unsigned int totalHydrationWaters){

    float water_asf = functions_WETSAXS::asf_at_q_zero(99);

    if (totalAtoms < 1 || forwardScatteringExVolume == 0)
        throw std::overflow_error("Divide by zero exception");

    matchPoint = (forwardScatteringParticle + totalHydrationWaters*water_asf)/(forwardScatteringExVolume + totalHydrationWaters*water_asf); // divide by zero possible

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

/*
 * must be run after SphericalHarmonics
 */
void AtomisticModel::createSphericalBessels(unsigned int totalqvalues, std::vector < float > & qvalues) {

    sbj = SphericalBessels(lmax, totalqvalues, totalAtoms, qvalues, &rvalues);

}

const SphericalBessels &AtomisticModel::getSbj() const {
    return sbj;
}


/*
 * Total partials per q value is (lmax+1)^2
 * Total partials im array is totalqvalues*(lmax+1)*(lmax+1)
 * A_lm(q_i) index is given by (l*l + (l-m)*q_index
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

//    auto pAtomicVolumes = pdbModel.getAtomicVolumeVec();
//    std::vector<float> fractionOccAtoms(totalAtoms);
//    for(unsigned int i=0; i<totalAtoms; i++){
//        fractionOccAtoms[i] =  pAtomicVolumes[i]/29.9;
//    }

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
