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

#ifndef WETSAXS_HCPGRIDMODEL_H
#define WETSAXS_HCPGRIDMODEL_H


#include <cmath>
#include <sastools/vector3.h>
#include <string>
#include <map>
#include <vector>
#include <cfloat>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <iostream>

class HCPGridModel {

    const float radius; // spacing between points in lattice, not radius
    float kmax;
    const float sqrt3 = sqrtf(3);
    const float inv_rsqrt3, inv_radius, z_point;
    const float inv2sqrt6r, invsqrt3r, sqrt3inv3r, sqrt6r;
    float upper_limit;

    const vector3 a1;
    const vector3 a2;
    const vector3 a3;
    const vector3 basis2;

    std::map<std::string, float> incidents;
    std::vector< vector3> positions;
    std::vector< float> amplitudes;
    int index_max, index_min;
    std::vector< vector3> bases;
public:
    /*
     * radius is the distance between lattice points
     * c-alpha to c-alpha distance on average 5.5 Angstrom
     * radius => 1.8?
     */
    HCPGridModel(float r=1) : radius(r),
                            inv_radius(1/r),
                            inv_rsqrt3(1.0f/(radius * sqrt3)),
                            z_point(1.5f/sqrtf(6.0)*inv_radius),
                            a1(vector3(0.5*r, -0.5*sqrt3*r, 0)),
                            a2(vector3(0.5*r,  0.5*sqrt3*r, 0)),
                            a3(vector3(0, 0, (float)(2.0f*sqrt(6)/3.0f*r))),
                            basis2(a1*(2.0f/3.0f) + a2*(1.0f/3.0f) + a3/2.0f),
                            inv2sqrt6r(0.5/sqrtf(6.0)*inv_radius),
                            invsqrt3r(1.0f/sqrtf(3.0)*inv_radius),
                            sqrt3inv3r(sqrtf(3.0)/3.0f*radius),
                            sqrt6r(sqrtf(6.0)*radius){

        index_max = -INT16_MAX;
        index_min = INT16_MAX;
        kmax = FLT_MIN;
        upper_limit = 4*radius;
        bases.emplace_back(vector3(0,1,0));
        bases.emplace_back(vector3(1,1,0));
        bases.emplace_back(vector3(-1,1,0));
        bases.emplace_back(vector3(0,-1,0));
        bases.emplace_back(vector3(1,-1,0));
        bases.emplace_back(vector3(-1,-1,0));
        bases.emplace_back(vector3(-1,0,0));
        bases.emplace_back(vector3(1,0,0));

        bases.emplace_back(vector3(0,1,1));
        bases.emplace_back(vector3(1,1,1));
        bases.emplace_back(vector3(-1,1,1));
        bases.emplace_back(vector3(0,-1,1));
        bases.emplace_back(vector3(1,-1,1));
        bases.emplace_back(vector3(-1,-1,1));
        bases.emplace_back(vector3(-1,0,1));
        bases.emplace_back(vector3(1,0,1));
        bases.emplace_back(vector3(0,0,1));

        bases.emplace_back(vector3(0,1,-1));
        bases.emplace_back(vector3(1,1,-1));
        bases.emplace_back(vector3(-1,1,-1));
        bases.emplace_back(vector3(0,-1,-1));
        bases.emplace_back(vector3(1,-1,-1));
        bases.emplace_back(vector3(-1,-1,-1));
        bases.emplace_back(vector3(-1,0,-1));
        bases.emplace_back(vector3(1,0,-1));
        bases.emplace_back(vector3(0,0,-1));


    }

    float getRadius(){return radius;}
    float getInvRadius(){return inv_radius;}
    float getInvRSqrt3(){ return inv_rsqrt3;}

    void clearIncidents(){ incidents.clear();}


    /*
     * this lattice is sampling the positions of the atoms, so the distance between the atom and lattice point needs to be calculated
     * the model should be centered
     *
     */
    void addToHCPGrid(vector3 & point){

        int n1_index, n2_index, n3_index, n1_index_b2, n2_index_b2, n3_index_b2;

        const float * pX = &point.x;
        const float * pY = &point.y;
        const float * pZ = &point.z;

        float vecLen = point.length();
        if (vecLen > kmax){
            kmax = vecLen;
        }

        n3_index = (int)std::round(z_point* *pZ);
        n2_index = (int)std::round( (sqrt3* *pX + *pY)*inv_rsqrt3);
        n1_index = (int)std::round( (sqrt3* *pX - *pY)*inv_rsqrt3);

        setMinMax(n1_index, n2_index, n3_index);

        std::string index = std::to_string(n1_index).append("_").append(std::to_string(n2_index)).append("_").append(std::to_string(n3_index));
        auto mit = incidents.find(index);
        if (mit == incidents.end()){
            incidents.insert(std::make_pair(index,1.0));
        } else {
            mit->second += 1.0;
        }

        float partial_amplitude = 1/6.0;

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
//        index = std::to_string(n1_index).append("_").append(std::to_string(n2_index + 1)).append("_").append(std::to_string(n3_index));
//        mit = incidents.find(index);
//        if (mit == incidents.end()){
//            incidents.insert(std::make_pair(index,partial_amplitude));
//        } else {
//            mit->second += partial_amplitude;
//        }
//
//        // (1,1,0)
//        index = std::to_string(n1_index + 1).append("_").append(std::to_string(n2_index + 1)).append("_").append(std::to_string(n3_index));
//        mit = incidents.find(index);
//        if (mit == incidents.end()){
//            incidents.insert(std::make_pair(index,partial_amplitude));
//        } else {
//            mit->second += partial_amplitude;
//        }
//
//        // (1,0,0)
//        index = std::to_string(n1_index + 1).append("_").append(std::to_string(n2_index)).append("_").append(std::to_string(n3_index));
//        mit = incidents.find(index);
//        if (mit == incidents.end()){
//            incidents.insert(std::make_pair(index,partial_amplitude));
//        } else {
//            mit->second += partial_amplitude;
//        }
//
//        // (0,-1,0)
//        index = std::to_string(n1_index).append("_").append(std::to_string(n2_index-1)).append("_").append(std::to_string(n3_index));
//        mit = incidents.find(index);
//        if (mit == incidents.end()){
//            incidents.insert(std::make_pair(index,partial_amplitude));
//        } else {
//            mit->second += partial_amplitude;
//        }
//
//        // (-1,-1,0)
//        index = std::to_string(n1_index-1).append("_").append(std::to_string(n2_index-1)).append("_").append(std::to_string(n3_index));
//        mit = incidents.find(index);
//        if (mit == incidents.end()){
//            incidents.insert(std::make_pair(index,partial_amplitude));
//        } else {
//            mit->second += partial_amplitude;
//        }
//
//        // (-1,0,0)
//        index = std::to_string(n1_index-1).append("_").append(std::to_string(n2_index)).append("_").append(std::to_string(n3_index));
//        mit = incidents.find(index);
//        if (mit == incidents.end()){
//            incidents.insert(std::make_pair(index,partial_amplitude));
//        } else {
//            mit->second += partial_amplitude;
//        }
//
//        // (0,0,0) - basis
//        vector3 temp = (a1 * n1_index + a2 * n2_index + a3 * n3_index) - basis2;
//        n3_index_b2 = (int)std::round( (3*temp.z - sqrt6r) * inv2sqrt6r );
//        n2_index_b2 = (int)std::round( invsqrt3r * (sqrt3*temp.x + temp.y - sqrt3inv3r) );
//        n1_index_b2 = (int)std::round( invsqrt3r * (sqrt3*temp.x - temp.y - 2*sqrt3inv3r) );
//
//        // convert this to the basis indices
//        //index = std::to_string(n1_index).append("_").append(std::to_string(n2_index)).append("_").append(std::to_string(n3_index)).append("_-basis");
//        index = std::to_string(n1_index_b2).append("_").append(std::to_string(n2_index_b2)).append("_").append(std::to_string(n3_index_b2)).append("_basis");
//        mit = incidents.find(index);
//        if (mit == incidents.end()){
//            incidents.insert(std::make_pair(index,partial_amplitude));
//        } else {
//            mit->second += partial_amplitude;
//        }
//
//        // (1,1,0) - basis
//        temp = (a1 * (n1_index+1) + a2 * (n2_index+1) + a3 * n3_index) - basis2;
//        n3_index_b2 = (int)std::round( (3*temp.z - sqrt6r)*inv2sqrt6r );
//        n2_index_b2 = (int)std::round( invsqrt3r * (sqrt3*temp.x + temp.y - sqrt3inv3r) );
//        n1_index_b2 = (int)std::round( invsqrt3r * (sqrt3*temp.x - temp.y - 2*sqrt3inv3r) );
//
//        //index = std::to_string(n1_index+1).append("_").append(std::to_string(n2_index+1)).append("_").append(std::to_string(n3_index)).append("_-basis");
//        index = std::to_string(n1_index_b2).append("_").append(std::to_string(n2_index_b2)).append("_").append(std::to_string(n3_index)).append("_basis");
//        mit = incidents.find(index);
//        if (mit == incidents.end()){
//            incidents.insert(std::make_pair(index,partial_amplitude));
//        } else {
//            mit->second += partial_amplitude;
//        }
//
//
//        // (1,0,0) - basis
//        temp = (a1 * (n1_index+1) + a2 * (n2_index) + a3 * n3_index) - basis2;
//        n3_index_b2 = (int)std::round( (3*temp.z - sqrt6r)*inv2sqrt6r );
//        n2_index_b2 = (int)std::round( invsqrt3r * (sqrt3*temp.x + temp.y - sqrt3inv3r) );
//        n1_index_b2 = (int)std::round( invsqrt3r * (sqrt3*temp.x - temp.y - 2*sqrt3inv3r) );
//
//        index = std::to_string(n1_index_b2).append("_").append(std::to_string(n2_index_b2)).append("_").append(std::to_string(n3_index_b2)).append("_basis");
//        mit = incidents.find(index);
//        if (mit == incidents.end()){
//            incidents.insert(std::make_pair(index,partial_amplitude));
//        } else {
//            mit->second += partial_amplitude;
//        }
//
//        // (0,0,0) + basis
//        temp = (a1 * (n1_index) + a2 * (n2_index) + a3 * n3_index) + basis2;
//        n3_index_b2 = (int)std::round( (3*temp.z - sqrt6r)*inv2sqrt6r );
//        n2_index_b2 = (int)std::round( invsqrt3r * (sqrt3*temp.x + temp.y - sqrt3inv3r) );
//        n1_index_b2 = (int)std::round( invsqrt3r * (sqrt3*temp.x - temp.y - 2*sqrt3inv3r) );
//
//        index = std::to_string(n1_index_b2).append("_").append(std::to_string(n2_index_b2)).append("_").append(std::to_string(n3_index_b2)).append("_basis");
//        mit = incidents.find(index);
//        if (mit == incidents.end()){
//            incidents.insert(std::make_pair(index,partial_amplitude));
//        } else {
//            mit->second += partial_amplitude;
//        }
//
//        // (-1,0,0) + basis
//        temp = (a1 * (n1_index-1) + a2 * (n2_index) + a3 * n3_index) + basis2;
//        n3_index_b2 = (int)std::round( (3*temp.z - sqrt6r)*inv2sqrt6r );
//        n2_index_b2 = (int)std::round( invsqrt3r * (sqrt3*temp.x + temp.y - sqrt3inv3r) );
//        n1_index_b2 = (int)std::round( invsqrt3r * (sqrt3*temp.x - temp.y - 2*sqrt3inv3r) );
//
//        index = std::to_string(n1_index_b2).append("_").append(std::to_string(n2_index_b2)).append("_").append(std::to_string(n3_index_b2)).append("_basis");
//        mit = incidents.find(index);
//        if (mit == incidents.end()){
//            incidents.insert(std::make_pair(index,partial_amplitude));
//        } else {
//            mit->second += partial_amplitude;
//        }
//
//        // (-1,-1,0) + basis
//        temp = (a1 * (n1_index-1) + a2 * (n2_index-1) + a3 * n3_index) + basis2;
//        n3_index_b2 = (int)std::round( (3*temp.z - sqrt6r)*inv2sqrt6r );
//        n2_index_b2 = (int)std::round( invsqrt3r * (sqrt3*temp.x + temp.y - sqrt3inv3r) );
//        n1_index_b2 = (int)std::round( invsqrt3r * (sqrt3*temp.x - temp.y - 2*sqrt3inv3r) );
//
//        index = std::to_string(n1_index_b2).append("_").append(std::to_string(n2_index_b2)).append("_").append(std::to_string(n3_index_b2)).append("_basis");
//        mit = incidents.find(index);
//        if (mit == incidents.end()){
//            incidents.insert(std::make_pair(index,partial_amplitude));
//        } else {
//            mit->second += partial_amplitude;
//        }
    }

    // convert incidents to density map
    void convertToDensityMap(std::string outputfilename);

    void setMinMax(int n1_index, int n2_index, int n3_index){

        if (n3_index < index_min){
            index_min = n3_index;
        }

        if (n2_index < index_min){
            index_min = n2_index;
        }

        if (n1_index < index_min){
            index_min = n1_index;
        }

        if (n3_index > index_max){
            index_max = n3_index;
        }

        if (n2_index > index_max){
            index_max = n2_index;
        }

        if (n1_index > index_max){
            index_max = n1_index;
        }
    }

    float convolutionFunction(float length);

    std::string setMapAndHeaderParametersForXPLOR(int & na, int & nb, int & nc, float grid_spacing);
    void convertIncidentsToVectorAmplitudes();
};


#endif //WETSAXS_HCPGRIDMODEL_H
