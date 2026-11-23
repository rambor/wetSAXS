
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

#include <algorithm>
#include <chrono>

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
        thetas[n] = (rval > 0.0f) ? acosf(zval/rval) : 0.0f;
        phis[n] = (xval == 0 || xval == 0.0f) ? 0.0f : atan2f(yval, xval); // need to check if this is a good way to test for 0

        atomList.insert(pAtomicNumbers[n]); // unique atomic numbers
    }

    // need to reassign atomic numbers to include hydrogens

    atomList.insert(99);  // 99 is unique for water
}

/*
 * Partial amplitudes for the particle (A_p) and its excluded volume (A_x).
 *
 * Total partials per q value is (lmax+1)^2
 * Total partials in array is totalqvalues*(lmax+1)*(lmax+1)
 * A_lm(q_i) index is (lmax+1)^2*q_index + l*l + l + m
 *
 * A_p : real atomic form factors (unified atoms)
 * A_x : displaced bulk solvent at the same coordinates. Each atom carries a smooth (Gaussian) blob of
 *       bulk water with the atom's volume V_k, written as a Gaussian water minus the part of the water
 *       that does not fit:
 *           g_k(q) = G(q; 29.9) - G(q; 29.9 - V_k),   G(q; V) = 0.334 V exp(-q^2 V^(2/3) / 4pi)
 *       At q = 0, g_k = 0.334 V_k (electrons of solvent displaced).
 *       A_x is stored with a POSITIVE sign; Babinet's minus sign is applied in Fit (-c_x A_xw).
 *
 * Only the particle's own atoms (n < pdbModel.getTotalCoordinates()) are summed. Hydration waters are
 * handled by the Waters class (A_w and their displaced bulk water), never here.
 *
 * The solvent-fluctuation term (g(r) of water) is an intensity and is NOT part of these amplitudes.
 */
void AtomisticModel::calculatePartialAmplitudes(unsigned int lmax,
                                                unsigned int totalqvalues,
                                                std::vector < float > & qvalues,
                                                bool recalculate){

    const unsigned int totalPDBAtoms = pdbModel.getTotalCoordinates();
    this->lmax = lmax;
    const unsigned int lmax_plus_1 = lmax + 1;

    if (recalculate){
        sbj.recalculate(totalqvalues, qvalues, &rvalues);
    } else {
        this->createSphericalHarmonics(lmax);               // depends on (theta, phi)
        this->createSphericalBessels(totalqvalues, qvalues); // depends on q and r
    }

    const unsigned int partialsPerQ = lmax_plus_1*lmax_plus_1;
    almRealSum.resize(partialsPerQ*totalqvalues);
    almImagSum.resize(partialsPerQ*totalqvalues);
    almRealExcludedVolSum.resize(partialsPerQ*totalqvalues);
    almImagExcludedVolSum.resize(partialsPerQ*totalqvalues);

    float * const ptrRAlm = almRealSum.data();
    float * const ptrIAlm = almImagSum.data();
    float * const ptrREx = almRealExcludedVolSum.data();
    float * const ptrIEx = almImagExcludedVolSum.data();

    const float * pSBJ = sbj.getpSphericalBessels();
    const int * pASF = pdbModel.getAdjustedASFNumbers();
    const float * pVol = pdbModel.getAtomicVolumeVec();
    const float * pReal = she.getpDataYLMReal();
    const float * pImag = she.getpDataYLMImag();

    // unique ASF codes -> slot, so each thread evaluates each form factor once per q
    const std::vector<int> uniqueASF(atomList.begin(), atomList.end());
    std::vector<int> slotOfAtom(totalPDBAtoms);
    for(unsigned int n = 0; n < totalPDBAtoms; n++){
        slotOfAtom[n] = (int)(std::lower_bound(uniqueASF.begin(), uniqueASF.end(), pASF[n]) - uniqueASF.begin());
    }

    const float rho0 = 0.334f;                 // bulk water, electrons per A^3
    const float vWater = 29.9f;                // volume of one water, A^3
    const float inv4pi = 1.0f/(4.0f*(float)M_PI);
    const float waterWidth = std::cbrt(vWater*vWater)*inv4pi;

    // q-independent parts of the leftover Gaussian, G(q; dV) = leftAmp * exp(-q^2 * leftWidth),
    // with dV = 29.9 - V_k (may be negative for atoms larger than a water)
    std::vector<float> leftAmp(totalPDBAtoms), leftWidth(totalPDBAtoms);
    for(unsigned int n = 0; n < totalPDBAtoms; n++){
        const float dV = vWater - pVol[n];
        leftAmp[n]   = rho0*dV;
        leftWidth[n] = std::cbrt(dV*dV)*inv4pi;
    }

    // stride of the bessel and ylm tables is totalAtoms (may include appended waters)
    const unsigned long lmax_plus_1_totalAtoms = (unsigned long)lmax_plus_1*totalAtoms;
    const int totalQ = (int)qvalues.size();

    auto startTime = std::chrono::steady_clock::now();

    // each q writes only its own block of partials: q values are independent
    #pragma omp parallel
    {
        std::vector<float> ffq(uniqueASF.size()), fp(totalPDBAtoms), fx(totalPDBAtoms);

        #pragma omp for schedule(dynamic)
        for(int q_index = 0; q_index < totalQ; q_index++){

            const float qvalue = qvalues[q_index];
            for(size_t k = 0; k < uniqueASF.size(); k++){
                ffq[k] = functions_WETSAXS::asf(uniqueASF[k], qvalue);
            }

            const float q2 = qvalue*qvalue;
            const float gWater = rho0*vWater*std::exp(-q2*waterWidth);

            for(unsigned int n = 0; n < totalPDBAtoms; n++){
                fp[n] = ffq[slotOfAtom[n]];
                fx[n] = gWater - leftAmp[n]*std::exp(-q2*leftWidth[n]);
            }

            const unsigned int ll_index = partialsPerQ*(unsigned int)q_index;

            for(int l=0; l<=(int)lmax; l++){

                const unsigned long bessel_index = lmax_plus_1_totalAtoms*q_index + (unsigned long)totalAtoms*l;
                const unsigned long base_index = (unsigned long)(l*l)*totalAtoms;
                const unsigned int partials_index = ll_index + l*l + l;

                // m = 0 : imaginary part of Y_l0 is zero
                {
                    double sumRAp = 0.0, sumREx = 0.0;
                    unsigned long ylm_index = base_index + (unsigned long)l*totalAtoms;

                    for(unsigned int n = 0; n < totalPDBAtoms; n++, ylm_index++){
                        const float by = pSBJ[bessel_index + n]*pReal[ylm_index];
                        sumRAp += fp[n]*by;
                        sumREx += fx[n]*by;
                    }

                    ptrRAlm[partials_index] = (float)sumRAp;
                    ptrIAlm[partials_index] = 0.0f;
                    ptrREx[partials_index] = (float)sumREx;
                    ptrIEx[partials_index] = 0.0f;
                }

                // m > 0, and m < 0 by symmetry of a real density: A_{l,-m} = (-1)^m conj(A_{l,m})
                for(int m=1; m<=l; m++) {

                    double sumRAp = 0.0, sumIAp = 0.0, sumREx = 0.0, sumIEx = 0.0;
                    unsigned long ylm_index = base_index + (unsigned long)(l + m)*totalAtoms;

                    for(unsigned int n = 0; n < totalPDBAtoms; n++, ylm_index++){
                        const float b = pSBJ[bessel_index + n];
                        const float yR = b*pReal[ylm_index];
                        const float yI = b*pImag[ylm_index];
                        sumRAp += fp[n]*yR;
                        sumIAp += fp[n]*yI;
                        sumREx += fx[n]*yR;
                        sumIEx += fx[n]*yI;
                    }

                    const float rAp = (float)sumRAp, iAp = (float)sumIAp, rEx = (float)sumREx, iEx = (float)sumIEx;

                    unsigned int idx = partials_index + m;
                    ptrRAlm[idx] = rAp;  ptrIAlm[idx] = iAp;
                    ptrREx[idx]  = rEx;  ptrIEx[idx]  = iEx;

                    idx = partials_index - m;
                    if (m & 1) { // odd m: flip sign of the real part
                        ptrRAlm[idx] = -rAp; ptrIAlm[idx] = iAp;
                        ptrREx[idx]  = -rEx; ptrIEx[idx]  = iEx;
                    } else {     // even m: flip sign of the imaginary part
                        ptrRAlm[idx] = rAp;  ptrIAlm[idx] = -iAp;
                        ptrREx[idx]  = rEx;  ptrIEx[idx]  = -iEx;
                    }
                }
            }
        }
    }

    try {
        setMatchPointScale(pdbModel.getAdjustedASFNumbers());
    } catch (std::overflow_error &e ){
        std::cout << e.what() << " :: ";
    }

    const double runtime = std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime).count();
    logger("Total TIME (Model)", formatNumber((float)runtime,8));
}

/*
 * Forward-scattering amplitudes (electrons) of the particle and of its excluded volume.
 *   particle        : sum of f_k(0) over the particle's atoms (unified atoms include H)
 *   excluded volume : 0.334 * sum of V_k  (= A_x at q = 0)
 * matchPoint = their ratio, i.e. mean particle density / bulk water density (proteins ~1.28).
 * It is reported for information. With A_x in electrons, c_x is the buffer density relative to
 * 0.334 e/A^3 directly and must NOT be multiplied by matchPoint in Fit.
 * Hydration waters are not included here.
 */
void AtomisticModel::setMatchPointScale(int * pAtomicNumbers){

    const unsigned int totalParticleAtoms = pdbModel.getTotalCoordinates();
    const float * pVol = pdbModel.getAtomicVolumeVec();

    double sumF = 0.0, sumV = 0.0;
    for(unsigned int n=0; n < totalParticleAtoms; n++){
        sumF += functions_WETSAXS::asf_at_q_zero(pAtomicNumbers[n]);
        sumV += pVol[n];
    }

    forwardScatteringParticle = (float)sumF;            // amplitude, electrons
    forwardScatteringExVolume = (float)(0.334*sumV);     // amplitude, electrons

    if (totalParticleAtoms < 1 || forwardScatteringExVolume == 0)
        throw std::overflow_error("Divide by zero exception");

    SASTOOLS_UTILS_H::logger("Particle electrons", formatNumber(forwardScatteringParticle, 1));
    SASTOOLS_UTILS_H::logger("Excluded volume electrons", formatNumber(forwardScatteringExVolume, 1));
    SASTOOLS_UTILS_H::logger("density (en-/Angstrom^3)", formatNumber(forwardScatteringParticle/(float)sumV, 3));

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

    if ( totalWaters < 1 ) {
        throw std::invalid_argument( "Water model not initialized" );
    }

    // waters are appended after the particle's atoms; resetting (not accumulating) makes repeated
    // calls safe, e.g. once per round in Fit::chiFreeSearch
    unsigned int startHere = pdbModel.getTotalCoordinates();
    totalAtoms = startHere + totalWaters;
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
// float AtomisticModel::addUp(const unsigned int total_width, float * vector_to_sum) {
//
//     unsigned int window = 8; // window is 8 for _m128 and 16 for _m256
//
//     const unsigned int alignedN = total_width - total_width % window;
//     simde__m128 mmSum1, mmSum2, v0, v1;
//
//     auto * sum_sse = (float *)_aligned_malloc(sizeof(float)*4, 16);
//
//     for (int i = 0; i < alignedN; i+=window) {
//         v0 = simde_mm_load_ps(&vector_to_sum[i]);
//         mmSum1 = simde_mm_add_ps(mmSum1, v0);
//
//         v1 = simde_mm_load_ps(&vector_to_sum[i+4]);
//         mmSum1 = simde_mm_add_ps(mmSum1, v1);
//     }
//
//     simde_mm_store_ps(sum_sse, mmSum1);
//
//     float sum = sum_sse[0] + sum_sse[1] + sum_sse[2] + sum_sse[3];
//
//     for (int i = alignedN; i < total_width; ++i) {
//         sum += vector_to_sum[i];
//     }
//
//     return sum;
// }

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
