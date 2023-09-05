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

#ifndef WETSAXS_SPHERICALHARMONICS_H
#define WETSAXS_SPHERICALHARMONICS_H

#include <math.h>
#include <vector>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <sastools/Coords.h>
#include <sastools/utils.h>

class SphericalHarmonics {

private:

    unsigned int lmax, numAtoms;
    unsigned int ylm_size;

    std::vector<float> y_lm_real;
    std::vector<float> y_lm_imag;


public:
    SphericalHarmonics()=default;
    SphericalHarmonics(unsigned int lmax, unsigned int numAtoms);

    //    SphericalHarmonics(int lmax, int numAtoms, std::vector<Coords> & coords);

    // copy constructor - prevent copying
    SphericalHarmonics(const SphericalHarmonics & model)= delete;
    // copy assignment - prevent copying
    SphericalHarmonics & operator=(const SphericalHarmonics & model)= delete;

    // move assignment operator
    SphericalHarmonics & operator=(SphericalHarmonics && model) noexcept {

        if (&model == this){
            return *this;
        }

        lmax = std::move(model.lmax);
        numAtoms = std::move(model.numAtoms);
        ylm_size = std::move(model.ylm_size);

        y_lm_real = std::move(model.y_lm_real);
        y_lm_imag = std::move(model.y_lm_imag);

        return *this;
    }

    SphericalHarmonics (SphericalHarmonics && model) noexcept {
        *this = std::move(model);
    }

    float * getpDataYLMReal() { return y_lm_real.data();}
    float * getpDataYLMImag() { return y_lm_imag.data();}


    int getYLMSize() { return ylm_size; }

    void populateSHETable( const float * pThetas,  const float * pPhis);

//    void calculateSHE(int lvalue, int mvalue, const float *pThetas, const float *pPhis);
};


#endif //WETSAXS_SPHERICALHARMONICS_H
