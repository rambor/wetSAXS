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
// Created by Robert Rambo on 20/05/2022.
//

#ifndef WETSAXS_SPHERICALBESSELS_H
#define WETSAXS_SPHERICALBESSELS_H


#include <vector>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <boost/math/special_functions/bessel.hpp>
#include <sastools/Coords.h>
#include <sastools/simde-no-tests-master/simde-arch.h>
#include <sastools/simde-no-tests-master/simde-common.h>
#include <sastools/vector3.h>
#include <sastools/utils.h>

class SphericalBessels {

private:

    int bessel_size;
    int lmax;
    int qvaluesSize;
    int numAtoms;

    //float * besselalq; // aligned memory

    std::vector<float> besselalq;
    void calculate(std::vector<float> &qr);

public:
    SphericalBessels();

    SphericalBessels(int bs, int qs, int numatoms, std::vector < float > & qvalues, std::vector < float > * rvalues);
    SphericalBessels(int bs, int qs, int na, std::vector < float > & qvalues, std::vector<Coords> & coords);

    ~SphericalBessels(){
        //_aligned_free(besselalq);
    }

    // copy constructor - prevent copying
    SphericalBessels(const SphericalBessels & model)= delete;
    // copy assignment - prevent copying
    SphericalBessels & operator=(const SphericalBessels & model)= delete;

    // move assignment operator
    SphericalBessels & operator=(SphericalBessels && model) noexcept {

        if (&model == this){
            return *this;
        }

        bessel_size = std::move(model.bessel_size);
        lmax = std::move(model.lmax);
        qvaluesSize = std::move(model.qvaluesSize);
        numAtoms = std::move(model.numAtoms);
        //std::swap(besselalq, model.besselalq);
        besselalq = std::move(model.besselalq);
        return *this;
    }

    SphericalBessels (SphericalBessels && model) noexcept {
        *this = std::move(model);
        //std::swap(*this, model);
    }


    float * getpSphericalBessels() { return besselalq.data();}

    int getSize(){return bessel_size;}

    // qs total values in qvalues

    void recalculate(int totalQvalues, std::vector<float> & qvalues, std::vector < float > * rvalues);
};


#endif //WETSAXS_SPHERICALBESSELS_H
