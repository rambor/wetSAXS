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


#include <sastools/utils.h>
#include <sastools/FileClass.h>
#include <sastools/PDBModel.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <cfloat>
#include "Aligner.h"
#include "SVD.h"
#include "HCPGridModel.h"


Aligner::Aligner(std::string filename, std::string directory) {

    // open file and read in names, check if PDB
    std::cout << "Reading COM file" << std::endl;

    if  (!boost::filesystem::exists( filename )){
        throw std::invalid_argument("  => CANNOT READ COMMAND FILE : " + filename);
    }

    /*
     * Read resid range and filenames
     * Load/Associate additional information with chains
     * Assign ss keys
     */
    std::vector<std::string> com_lines;
    std::ifstream data (filename.c_str());

    if (data.is_open()) {
        boost::regex hashtag("#");
        std::string line;

        while(!data.eof()) //
        {
            getline(data, line); //this function grabs a line and moves to next line in file

            if (std::isspace(line[0])){ // remove any leading white space
                line.erase(line.begin(), std::find_if(line.begin(), line.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
            }

            boost::trim(line);
            if (line[0] != '#' && line.size() > 2){ // parse
                com_lines.push_back(line);
            }
        }
    }
    data.close();

    logger("Total Lines Read", std::to_string(com_lines.size()));
    boost::trim(com_lines[0]);
    std::vector<std::string> elements;
    boost::split(elements, com_lines[0], boost::is_any_of("\t  "), boost::token_compress_on);

    boost::regex toFormat("to", boost::regex::icase);

    if (boost::regex_match( elements[0], boost::regex("ALIGNMENT_RESIDS", boost::regex::icase)) &&
            boost::regex_search( elements[1], boost::regex("to", boost::regex::icase))) {
            //10to101A
        std::vector<std::string> resis;
        boost::split(resis, elements[1], boost::is_any_of("to"), boost::token_compress_on);
        startResID = std::stoi(resis[0]);// << " -- " << resis[1];

        if (isalpha(resis[1][ resis[1].size() -1 ]) == 0){ //check that chain is specified
            throw std::invalid_argument("  => Missing Chain ID : " + resis[1]);
        } else {
            chainID = resis[1][ resis[1].size() -1 ];
            resis[1].pop_back();
            endResID = std::stoi(resis[1]);
        }

    } else {
        throw std::invalid_argument("  => MISSING : ALIGNMENT_RESIDS");
    }

    // check if directory exists and is writeable
    if (!boost::filesystem::exists(directory)){
        throw std::invalid_argument("** ERROR FILE => DIRECTORY NOT FOUND : " + directory);
    } else {
        this->directory = directory;
    }

    com_lines.erase(com_lines.begin());

    for(auto & ll : com_lines){
        try{
            FileClass basefile = FileClass(directory+"/"+ll);

            if (basefile.isPDB()){
                pdbs.push_back(directory+"/"+ll);
            }

        } catch (const std::invalid_argument& fl){
            std::cerr << " " << fl.what() << std::endl;
        }
    }

    totalPDBs = pdbs.size();
    logger("Total files to align", std::to_string(totalPDBs));
}

void Aligner::align(){
    // only align on the calphas
    std::random_device rd;
    std::mt19937 gen(rd());
    PDBModel reference(pdbs[0], false, false);
    unsigned int totalAtoms = reference.getTotalCoordinates();

    // need the CA within specified range

    auto pAtomTupe = reference.getPointerToAtomTypes();
    auto pResIDs = reference.getResIDIterator();

    std::vector<vector3> ref_vec;
    std::vector<vector3> ref_all_ca_vec;
    vector3 row1, row2, row3;
    std::vector<unsigned int> ref_indices;
    std::vector<unsigned int> indices;

    unsigned int counter=0;
    for(unsigned int i=0; i<totalAtoms; i++){
        std::string att = pAtomTupe[i];
        boost::trim(att);
        if (att== "CA"){
          if (pResIDs[i] >= startResID && pResIDs[i] <= endResID){
              ref_indices.push_back(pResIDs[i]);
              ref_vec.emplace_back(vector3(reference.getX()[i],reference.getY()[i], reference.getZ()[i]));
              indices.push_back(counter);
              counter++;
          }
        }
    } // all coordinates should be centered in the output ?

    // for each remaining pdb in pdbs, align to reference
    boost::numeric::ublas::matrix<float> rotation;
    vector3 centering_vec_tar, centering_vec_ref;

    HCPGridModel hcpGridModel(5.5/2);

    for(unsigned int i=1; i<totalPDBs;i++){
        PDBModel target(pdbs[i], false, false);
        // need the CA within specified range
        pAtomTupe = target.getPointerToAtomTypes();
        pResIDs = target.getResIDIterator();

        vector3 tar_center(0,0,0);
        std::vector<vector3> tar_vec;

        // should be in register with ref_indices, for each refID, find the one in target
        unsigned startAt = 0;
        for (auto & refID : ref_indices){
            for(unsigned int t=startAt; t<totalAtoms; t++){
                std::string att = pAtomTupe[t];
                boost::trim(att);
                if (att== "CA" && pResIDs[t] == refID) {
                      tar_vec.emplace_back(vector3(target.getX()[t], target.getY()[t], target.getZ()[t]));
                      startAt+=1;
                      break;
                }
            }
        }

        // perform the median based alignment
        this->median_SVD_alignment(gen,
                            indices,
                            ref_vec,
                            tar_vec,
                            rotation,
                            centering_vec_ref,
                            centering_vec_tar);

        // rotate the tar_all_ca_vec
        std::vector<vector3> tar_all_ca_vec;
        for(unsigned int t=0; t<totalAtoms; t++){
            std::string att = pAtomTupe[t];
            boost::trim(att);
            if (att == "CA" || att == "N" || att == "C" || att == "O") {
                tar_all_ca_vec.emplace_back(vector3(target.getX()[t], target.getY()[t], target.getZ()[t]));
            }
        }

        // write to file?
        row1 = vector3(rotation(0,0), rotation(0,1), rotation(0,2));
        row2 = vector3(rotation(1,0), rotation(1,1), rotation(1,2));
        row3 = vector3(rotation(2,0), rotation(2,1), rotation(2,2));

        // apply rotation to entire coordinates and calculate median RMSD
        unsigned int count = 1;
        for (auto & vec : tar_all_ca_vec) {
            vec -= centering_vec_tar;
            // rotate and translate to reference
            vec = vector3(row1.dot(vec), row2.dot(vec), row3.dot(vec)) + centering_vec_ref - *reference.getCenteringVector();
            // convert and update HCP grid model
            hcpGridModel.addToHCPGrid(vec);

//            char buffer [50];
//
//            sprintf(buffer, "%-6s%5i %4s %3s %1s%4s    %8.3f%8.3f%8.3f  1.00100.00\n", "ATOM", count, " O  ", "CA", chainID.c_str(),
//                    "A",
//                    vec.x,
//                    vec.y,
//                    vec.z);
//            std::cout << buffer;
//            count +=1;
        }
    }

    logger("FNISHED","ALIGNMENT");
    // create density map from HCP grid model
    logger("CREATING","HCP Model");
    hcpGridModel.convertToDensityMap("temp");
}


void Aligner::median_SVD_alignment(std::mt19937 & gen,
                                   std::vector<unsigned int> & indices,
                                   std::vector<vector3> &ref,
                                   std::vector<vector3> &tar,
                                   boost::numeric::ublas::matrix<float> & rotation,
                                   vector3 & centering_vec_ref,
                                   vector3 & centering_vec_tar) {

    // randPercent = rand(4..(totalCommon*0.51).to_i)
    unsigned int totalAtomsInUse = ref.size();
    unsigned int totalToSelectFrom = std::ceil(totalAtomsInUse * 0.51);

    unsigned int median_pos = (totalAtomsInUse % 2 == 0) ? (totalAtomsInUse/2 - 1) : (totalAtomsInUse-1)/2;

    std::vector<float> residuals(totalAtomsInUse);
    float d_x, d_y, d_z;
    std::uniform_int_distribution<unsigned int> distribution(4, totalToSelectFrom);

    vector3 ref_vec, tar_vec;
    vector3 row1, row2, row3;

    float residual, bestResidual = FLT_MAX;

    boost::numeric::ublas::matrix<float> bestRotation;
    boost::numeric::ublas::matrix<float> i_m(3,3);
    i_m(0,0) = 1.0;
    i_m(1,1) = 1.0;
    i_m(2,2) = 1.0;

    vector3 best_center;

    // select random subset of coordinates
    for(unsigned int r=0; r<totalRounds; r++){
        unsigned int useThis = distribution(gen);

        std::shuffle(indices.begin(), indices.end(), gen);
        std::sort(indices.begin(), indices.begin() + useThis);

        //assemble matrix
        boost::numeric::ublas::matrix<float> p_ref(useThis, 3);
        boost::numeric::ublas::matrix<float> q_tar(3, useThis);

        vector3 ref_center(0,0,0);
        vector3 tar_center(0,0,0);

        for(unsigned int j=0; j<useThis; j++){
            ref_center += ref[indices[j]];
            tar_center += tar[indices[j]];
        }

        ref_center /= (float)useThis;
        tar_center /= (float)useThis;

        // fill matrix after centering coordinates
        for(unsigned int j=0; j<useThis; j++){
            ref_vec = ref[indices[j]] - ref_center;
            p_ref(j,0) = ref_vec.x;
            p_ref(j,1) = ref_vec.y;
            p_ref(j,2) = ref_vec.z;

            tar_vec = tar[indices[j]] - tar_center;
            q_tar(0,j) = tar_vec.x;
            q_tar(1,j) = tar_vec.y;
            q_tar(2,j) = tar_vec.z;
        }

        /*
        * # A = U S V^T
        * # p reference
        * # from target to reference
        * # covariance matrix
        */
        //a = q*p.transpose;
        boost::numeric::ublas::matrix<float> a_matrix = boost::numeric::ublas::prod(q_tar, p_ref);
        this->getRotationMatrix(a_matrix, i_m, rotation);

        row1 = vector3(rotation(0,0), rotation(0,1), rotation(0,2));
        row2 = vector3(rotation(1,0), rotation(1,1), rotation(1,2));
        row3 = vector3(rotation(2,0), rotation(2,1), rotation(2,2));

        // apply rotation to entire coordiantes and calculate median RMSD
        for(unsigned int j=0; j<totalAtomsInUse; j++){

            ref_vec = ref[j] - ref_center;
            tar_vec = tar[j] - tar_center;

            d_x = ref_vec.x - row1.dot(tar_vec);
            d_y = ref_vec.y - row2.dot(tar_vec);
            d_z = ref_vec.z - row3.dot(tar_vec);

            residuals[j] = d_x*d_x + d_y*d_y + d_z*d_z;
        }

        std::sort(residuals.begin(), residuals.end());
        residual = ( median_pos & 1) ? residuals[median_pos] : (residuals[median_pos] + residuals[median_pos+1])*0.5;

        // get median residual
        if (residual < bestResidual){
            bestResidual = residual;
            bestRotation = boost::numeric::ublas::matrix(rotation);
            centering_vec_tar = tar_center;
            centering_vec_ref = ref_center;
            std::cout << r << " " << "Best Residual " << bestResidual << std::endl;
        }
    }

    // only c-alphas positions
    std::cout << "Best Residual " << bestResidual;
    // rotate the coordinates
    // write to file
    rotation = bestRotation;
}

void Aligner::getRotationMatrix(boost::numeric::ublas::matrix<float> & a_matrix,
                                boost::numeric::ublas::matrix<float> & i_m,
                                boost::numeric::ublas::matrix<float> & rotation) {

    boost::numeric::ublas::matrix<float> um;
    boost::numeric::ublas::matrix<float> sm;
    boost::numeric::ublas::matrix<float> vm;
    boost::numeric::ublas::matrix<float> v_u;
    boost::numeric::ublas::matrix<float> astar;
    boost::numeric::ublas::matrix<float> vm_i_m;
    boost::numeric::ublas::matrix<float> u_transpose;

    SVD svd;
    svd.svd(a_matrix, um, sm, vm);
    svd.pseudo_inverse(3, astar, um, sm, vm);

    // u, v, s = a.SV_decomp
    // sgn = (v*u.transpose).det <=> 0
    // m[3-1, 3-1] = sgn
    // rot = (v*m)*u.transpose
    u_transpose = boost::numeric::ublas::trans(um);
    v_u = boost::numeric::ublas::prod(vm, u_transpose);
    i_m(2,2) = svd.sgn(svd.determinant(v_u));
    vm_i_m = boost::numeric::ublas::prod(vm, i_m);
    rotation = boost::numeric::ublas::prod(vm_i_m, u_transpose);
}
