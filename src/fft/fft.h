#pragma once

#include <vector>

#include "Log.hpp"
#include "kiss_fft.h"

#define FFT_R(a, i) *(((kiss_fft_scalar*)&a) + (i << 1))
#define FFT_I(a, i) *(((kiss_fft_scalar*)&a) + (i << 1) + 1)

struct fftData {
    float frequency;
    float amplitude;
    float phase;
};
using fftDataBlob = std::vector<fftData>;

enum FFTDirection {
    In = 1,
    Out = 0
};

class fft_wrapper
{
public:
    fft_wrapper();
    fft_wrapper(uint32_t windowSize);
    ~fft_wrapper();

    void init(uint32_t windowSize);

    void computeA(kiss_fft_cpx* data, size_t frameSize = 0);
    void computeInverseA(kiss_fft_cpx* data, size_t frameSize = 0);

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
