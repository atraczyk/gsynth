#pragma once

#include <vector>

#include "Log.hpp"
#include "kiss_fft.h"

#ifndef M_PI
 #define M_PI   3.14159265358979323846264338327950288
#endif
#ifndef M_PI_2
#define M_PI_2  1.57079632679489661923132169163975144
#endif
#ifndef M_2PI
#define M_2PI   6.28318530717958647692528676655900576
#endif


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
