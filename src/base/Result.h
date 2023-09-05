// Copyright (c) 2023.
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
// Created by Robert Rambo on 11/01/2023.
//

#ifndef WETSAXS_RESULT_H
#define WETSAXS_RESULT_H

#include <vector>

class Result {

public:
    float score;

private:

    float scale;
    float cx;
    float bfactor;
    float dw;         // Durbin-Watson
    float chi;        // chi
    std::vector<float> i_calc; // model
    std::vector<unsigned int> selectedIndices; // indices used by Fit that correspond

public:

    Result(){};

    Result(float score, float scale, float cx, float b, float dw, float chi, std::vector<float> & icalc, const std::vector<unsigned int> & indices) : score(score), scale(scale), cx(cx), bfactor(b), dw(dw), chi(chi) {

        this->i_calc.swap(icalc); // if i_calc is empty, now guarantee icalc is empty

        unsigned int total = i_calc.size();
        selectedIndices.resize(total);
        const unsigned int * index = indices.data();
        unsigned int * sindex = selectedIndices.data();

        for(unsigned int i; i<total; i++){
            sindex[i] = index[i];
        }
    };

    float getScore(){ return score;}
    float getDurbinWatson(){ return dw;}
    float getChi(){ return chi;}
    float getBfactor(){ return bfactor;}
    float getCX(){ return cx; }

    std::vector<float> & getICalc(){ return i_calc;}
    std::vector<unsigned int> & getSelectedIndices(){ return selectedIndices;}

};

#endif //WETSAXS_RESULT_H
