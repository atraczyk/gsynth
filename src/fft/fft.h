#pragma once

#include <vector>

#include "Log.hpp"
#include "kiss_fft.h"

struct fftData {
    double frequency;
    double amplitude;
    double phase;
};
using fftDataBlob = std::vector<fftData>;

class fft_wrapper
{
public:
    fft_wrapper();
    fft_wrapper(uint32_t sampleRate, uint32_t windowSize);
    ~fft_wrapper();

    void init(uint32_t sampleRate, uint32_t windowSize);
    void setRealInput(const double * realData);
    fftDataBlob computeStft(bool useThreshold = false);
    std::vector<double> computeInverseStft();

    uint32_t getWindowSize() {
        return windowSize_;
    }

private:
    std::vector<kiss_fft_cpx>   in_;
    std::vector<kiss_fft_cpx>   out_;
    std::vector<kiss_fft_cpx>   inverse_;

    kiss_fft_cfg    cfg_;
    kiss_fft_cfg    inverse_cfg_;
    uint32_t        windowSize_;
    uint32_t        sampleRate_;
};
