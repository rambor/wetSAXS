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
// Created by Robert Rambo on 28/02/2025.
//
#include "HCPGridModel.h"

/*
 * Each point in HCP lattice is used to sample from via Gaussian kernel
 * underlying density map will be grid points
 *
 * n?_index_max and _min must be set
 */
void HCPGridModel::convertToDensityMap(std::string name){

    int na, nb, nc;

    float grid_spacing = radius; // should match HCP grid

//    float inv_bandwidth = 1.0f/(cutoff/2.355f); // should be same as fwhm_sigma?
//    float inv_bandwidth = 1.0f/fwhm_sigma;
//    float inv_bandwidth_squared = 0.5f*inv_bandwidth*inv_bandwidth;
//    float inv_bandwidth_2PI = 1.0f/sqrtf(2.0f*M_PI)*inv_bandwidth;

    float mapSum = 0.0f, mapSumSquared = 0.0f;
    float mapCount = 0.0f;

    std::string tempHeader = setMapAndHeaderParametersForXPLOR(na, nb, nc, grid_spacing);
    convertIncidentsToVectorAmplitudes();

    std::cout << kmax << std::endl;
    std::cout << na << std::endl;
    std::cout << nb << std::endl;
    std::cout << nc << std::endl;
    std::cout << index_max << std::endl;
    std::cout << index_min << std::endl;
    std::cout << "Grid Spacing " << upper_limit << std::endl;
    std::cout << "Incidents " << incidents.size() << std::endl;

    char buffer[80];

//    std::vector<float> averages;
//    std::vector<vector3> coords_averages;
    std::map<std::string, float > kernel_mapping;
    std::map<std::string, int > kernel_mapping_counter;
    float val;

    std::vector<std::string> tempLine;
    boost::regex ifPlusB("\\basis"); // match any character

    auto * const pPosition = (positions.size() != 0) ? positions.data() : nullptr;
    float * const pAmp = (amplitudes.size() != 0) ? amplitudes.data() : nullptr;

    int n_x, n_y, n_z;
    std::vector<int> indices = {-1,0,1};

    for(unsigned int m=0; m<positions.size(); m++){
        const vector3 & pPos = pPosition[m];
        vector3 pVec = pPos/grid_spacing;

        n_x = std::round(pVec.x);
        n_y = std::round(pVec.y);
        n_z = std::round(pVec.z);

        for(int i=0; i<3; i++){
            int v1 = indices[i] + n_x;
            for(int j=0; j<3; j++){
                int v2 = indices[j] + n_y;
                for(int k=0; k<3; k++){
                    int v3 = indices[k] + n_z;
                    std::string key = std::to_string(v1)+"_"+std::to_string(v2)+"_"+std::to_string(v3);

                    vector3 grid_point = vector3(v1, v2, v3)*grid_spacing;
                    float length = (pPos - grid_point).length();
                    val = convolutionFunction(length);

                    auto it = kernel_mapping.find(key);
                    auto cit = kernel_mapping_counter.find(key);
                    if (it == kernel_mapping.end()){
                        kernel_mapping.emplace(key, pAmp[m]*val);
                        kernel_mapping_counter.emplace(key, 1);
                    } else {
                        it->second += pAmp[m]*val;
                        cit->second += 1;
                    }
                }
            }
        }
    }

    // find largest amplitude to normalize the map
    float maxAmp = -FLT_MIN;
    for(auto & mapping : kernel_mapping){
        float value = mapping.second/(float)kernel_mapping_counter.find(mapping.first)->second;
        if (value > maxAmp){
            maxAmp = value;
        }
    }



    // kmax is the longest radial distance from center
    for(int z=0; z<nc; z++){

        float zsection = -kmax + z*grid_spacing;
        n_z = std::round(zsection/grid_spacing);

        std::snprintf(buffer, 80, "%8i\n", z);
        tempHeader.append(buffer);
        unsigned int stringIndex = 0;

        for(int y=0; y<nb; y++){

            float ysection = -kmax + y*grid_spacing;
            n_y = std::round(ysection/grid_spacing);

            for(int x=0; x<na; x++){

                n_x = std::round((-kmax+x*grid_spacing)/grid_spacing);

                std::string key = std::to_string(n_x)+"_"+std::to_string(n_y)+"_"+std::to_string(n_z);

                auto it = kernel_mapping.find(key);
                auto cit = kernel_mapping_counter.find(key);

                float kernelSum = 0.0f;

                if (it != kernel_mapping.end()){
                    kernelSum = it->second/(float)cit->second/maxAmp;
                }

                std::snprintf(buffer, 80, "%12.5E", kernelSum);
                tempHeader.append(buffer);
//                kernel_mapping.emplace(key, kernelSum);

                stringIndex += 1;

                if (stringIndex%6 == 0 ){
                    tempHeader += "\n";
                }

                mapSum += kernelSum;
                mapSumSquared += kernelSum*kernelSum;
                mapCount += 1.0;
            }

        }

        if (stringIndex % 6 != 0){
            tempHeader += "\n";
        }
    }

    std::snprintf(buffer, 80, "%8i\n",-9999);
    tempHeader.append(buffer);

    float ave = mapSum/mapCount;
    std::snprintf(buffer, 80, "%12.4E %12.4E\n", ave, std::sqrt(mapSumSquared/mapCount - ave*ave));
    tempHeader.append(buffer);

    // write to file
    std::string map = name + "_map.xplor";
    const char * outputFileName = map.c_str();

    FILE * pFile = fopen(outputFileName, "w");
    fprintf(pFile, tempHeader.c_str());
    fclose(pFile);
}

float HCPGridModel::convolutionFunction(float length){
    //return (cutoff - length)/cutoff; // linear model
    if (length <= upper_limit){ // if bead overlaps with atom, keep it -> cutoff = PI/qmax
        /*
         * inv_sigma = 1.0f/(2.0f*fwhm_sigma*fwhm_sigma);
         * fwhm_sigma = 1.5f*cutoff/2.355f;
         */
        return 1.0f/std::pow(length, 1.2);
        //return (cutoff - length)/cutoff;
    } else { // sets neighbors to include in calculation
        return 0;
    }

    // return expf(-0.5f*length*length/(fwhm_sigma*fwhm_sigma));
}


std::string HCPGridModel::setMapAndHeaderParametersForXPLOR(int & na, int & nb, int & nc, float grid_spacing){

//    std::cout << " BOUNDING BOX "<< std::endl;
//    std::cout << "   AXIS   MIN      MAX     LENGTH"  << std::endl;
    float cminx = -kmax;
    float cmaxx = kmax;
    float cminy = -kmax;
    float cmaxy = kmax;
    float cminz = -kmax;
    float cmaxz = kmax;
    float a_side = cmaxx - cminx;
    float b_side = cmaxy - cminy;
    float c_side = cmaxz - cminz;

//    printf("  => X %8.4f %8.4f %8.4f \n", cminx, cmaxx, a_side);
//    printf("  => Y %8.4f %8.4f %8.4f \n", cminy, cmaxy, b_side);
//    printf("  => Z %8.4f %8.4f %8.4f \n", cminz, cmaxz, c_side);

    auto startingNA = (int)std::round(cminx/grid_spacing);
    auto startingNB = (int)std::round(cminy/grid_spacing);
    auto startingNC = (int)std::round(cminz/grid_spacing);

    auto stoppingNA = (int)std::round(cmaxx/grid_spacing);
    auto stoppingNB = (int)std::round(cmaxy/grid_spacing);
    auto stoppingNC = (int)std::round(cmaxz/grid_spacing);

    na = std::abs(startingNA) + std::abs(stoppingNA) + 1;
    nb = std::abs(startingNB) + std::abs(stoppingNB) + 1;
    nc = std::abs(startingNC) + std::abs(stoppingNC) + 1;

    std::string tempHeader = "\n";
    tempHeader += "        4 !NTITLE\n";
    tempHeader += "REMARK 265 EXPERIMENTAL DETAILS\n";
    tempHeader += "REMARK 265\n";
    tempHeader += "REMARK 265 EXPERIMENT TYPE : X-RAY SOLUTION SCATTERING\n";
    tempHeader += "REMARK 265 DATA ACQUISITION\n";

    char buffer[80];
    std::snprintf(buffer, 80, "%8i%8i%8i%8i%8i%8i%8i%8i%8i\n", na, startingNA, stoppingNA, nb, startingNB, stoppingNB, nc, startingNC, stoppingNC);
    tempHeader.append(buffer);

    std::snprintf(buffer, 80, "%12.5E%12.5E%12.5E%12.5E%12.5E%12.5E \n",a_side, b_side, c_side, 90.0, 90.0, 90.0);

    tempHeader.append(buffer);
    // std::cout << tempHeader << std::endl;
    // write the matrix
    tempHeader += "ZYX\n";

    return tempHeader;
}

void HCPGridModel::convertIncidentsToVectorAmplitudes(){
    positions.clear();
    amplitudes.clear();

    std::vector<std::string> tempLine;
    boost::regex ifPlusB("\\basis"); // match any character
    int n1_index, n2_index, n3_index;

    for (auto & hcp_grid_point : incidents){
        tempLine.clear();
        boost::split(tempLine, hcp_grid_point.first, boost::is_any_of("_"), boost::token_compress_on);

        n1_index = (int)std::stoi(tempLine[0].c_str());
        n2_index = (int)std::stoi(tempLine[1].c_str());
        n3_index = (int)std::stoi(tempLine[2].c_str());

        //vector3 bead_vec = (a1 * n1_index + a2 * n2_index + a3 * n3_index);
        positions.emplace_back(vector3(a1 * n1_index + a2 * n2_index + a3 * n3_index));
        amplitudes.emplace_back(hcp_grid_point.second);

//        if (tempLine.size() > 3){ // determine if need to add or subtract the basis
//            positions.back() += basis2;
//        }
    }
}