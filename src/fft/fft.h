#pragma once

#include <vector>

#include "Log.hpp"
#include "kiss_fft.h"

using namespace std;

class fft_wrapper
{
public:
    fft_wrapper();
    fft_wrapper(uint32_t sampleRate, uint32_t windowSize, bool isInverse);
    ~fft_wrapper();

    void init(uint32_t sampleRate, uint32_t windowSize, bool isInverse);
    void setRealInput(const double * realData);
    std::vector<std::pair<double, double>> computeFrequencies(bool sort = false, int bandLimit = 0);

    uint32_t getWindowSize() {
        return windowSize_;
    }

private:
    std::vector<kiss_fft_cpx>   in_;
    std::vector<kiss_fft_cpx>   out_;
    std::vector<double>         mag_;

    kiss_fft_cfg    cfg_;
    uint32_t        windowSize_;
    uint32_t        sampleRate_;
};
