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
    float co;
    float dw;         // Durbin-Watson
    float chi;        // chi
    std::vector<float> i_calc; // model
    std::vector<unsigned int> selectedIndices; // indices used by Fit that correspond

public:

    Result() {
        std::cout << i_calc.size() << std::endl;
    }

    Result(float score, float scale, float cx, float co, float dw, float chi, std::vector<float> & icalc, const std::vector<unsigned int> & indices) :
    score(score),
    scale(scale),
    cx(cx),
    co(co),
    dw(dw),
    chi(chi) {

        if (!icalc.empty()) {
            this->i_calc.swap(icalc); // if i_calc is empty, now guarantee icalc is empty
        }

        unsigned int total = i_calc.size();
        selectedIndices = std::vector<unsigned int>(total);
        const unsigned int * index = indices.data();
        unsigned int * sindex = selectedIndices.data();

        for(unsigned int i; i<total; i++){
            sindex[i] = index[i];
        }
    };


    // copy constructor
    Result(const Result &src) = delete;

    // move constructor
    Result(Result &&src) {
        score = src.score;
        scale = src.scale;
        cx = src.cx;
        co = src.co;
        dw = src.dw;         // Durbin-Watson
        chi = src.chi;        // chi

        if (!src.i_calc.empty()) {
            i_calc = std::vector<float>();
            std::swap(src.i_calc, i_calc);
        }

        if (!src.selectedIndices.empty()) {
            selectedIndices = std::vector<unsigned int>();
            std::swap(src.selectedIndices, selectedIndices);
        }

    }

    ~Result()=default;

    // copy assignment operator
    Result& operator=(const Result &rhs) = delete;

    // move assignment operator // transfer ownership
    Result& operator=(Result &&rhs) noexcept {

        if (&rhs == this){
            return *this;
        }

        score = rhs.score;
        scale = rhs.scale;
        cx = rhs.cx;
        co = rhs.co;
        dw = rhs.dw;         // Durbin-Watson
        chi = rhs.chi;        // chi

        if (!rhs.i_calc.empty()) {
            i_calc = std::vector<float>();
            std::swap(rhs.i_calc, i_calc);
        }

        if (!rhs.selectedIndices.empty()) {
            selectedIndices = std::vector<unsigned int>();
            std::swap(rhs.selectedIndices, selectedIndices);
        }
        return *this;
    }

    float getScore(){ return score;}
    float getDurbinWatson(){ return dw;}
    float getChi(){ return chi;}
    float getCX(){ return cx; }
    float getCO(){ return co; }
    float getScale() {return scale; }


    std::vector<float> & getICalc(){ return i_calc;}
    std::vector<unsigned int> & getSelectedIndices(){ return selectedIndices;}

};

#endif //WETSAXS_RESULT_H
