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

enum FFTDirection {
    In = 1,
    Out = 0
};

class fft_wrapper
{
public:
    fft_wrapper();
    fft_wrapper(uint32_t sampleRate, uint32_t windowSize);
    ~fft_wrapper();

    void init(uint32_t sampleRate, uint32_t windowSize);
    void setData(   FFTDirection direction,
                    const std::vector<double>& realData,
                    const std::vector<double>& imagData = std::vector<double>());
    void setData(   FFTDirection direction,
                    const std::vector<kiss_fft_cpx>& realData);
    fftDataBlob computeStft();
    std::vector<kiss_fft_cpx> getRawOutput();
    std::vector<kiss_fft_cpx> getRawInverse();
    std::vector<double> computeInverseStft();

    void compute(std::vector<kiss_fft_cpx>& data);
    void computeInverse(std::vector<kiss_fft_cpx>& data);

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
