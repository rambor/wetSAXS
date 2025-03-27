// Copyright (c) 2025.
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
// Created by Robert Rambo on 26/02/2025.
//

#ifndef WETSAXS_ALIGNER_H
#define WETSAXS_ALIGNER_H

#include <string>
#include <iostream>
#include <boost/filesystem.hpp>
#include <utility>
#include <fstream>
#include <boost/regex.hpp>
#include <random>

// read in alignment file
// align all models in file
// create density map for ensembles
class Aligner {

private:
    std::vector<std::string> pdbs;
    unsigned int startResID, endResID, totalPDBs, totalRounds = 5000;
    std::string chainID;
    std::string directory;

public:
    Aligner(std::string filename, std::string directory);

    unsigned int getStartResID(){ return startResID;}
    unsigned int getEndResID(){ return endResID;}
    unsigned int getTotalPDBs() {return totalPDBs; }

    void align();

    void median_SVD_alignment(std::mt19937 & gen,
                              std::vector<unsigned int> & indices,
                              std::vector<vector3> & ref,
                              std::vector<vector3> & tar,
                              boost::numeric::ublas::matrix<float> & rotation,
                              vector3 & centering_vec_ref,
                              vector3 & centering_vec_tar);

    void getRotationMatrix(boost::numeric::ublas::matrix<float> &a_matrix,
                      boost::numeric::ublas::matrix<float> &i_m,
                      boost::numeric::ublas::matrix<float> &rotation);
};


#endif //WETSAXS_ALIGNER_H
