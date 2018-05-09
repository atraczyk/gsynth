#pragma once

#include <vector>

#include "fft.h"

constexpr const uint32_t maxFrameLength = 4096;

class PitchShifter
{
public:
    PitchShifter();
    PitchShifter(const uint32_t& sampleRate, const uint32_t& grainSize, const uint32_t& overlap);
    ~PitchShifter();

    void PitchShifter::pitchShift(
        float pitchShiftRatio,
        uint32_t grainSize_,
        long overlappingSamples,
        long nSamples,
        float sampleRate,
        float* source,
        float* out);

    void PitchShifter::pitchShift(
        float pitchShiftRatio,
        long nSamples,
        float* source,
        float* out);

    std::array<float, maxFrameLength>           srcBuffer_{};
    std::array<float, maxFrameLength>           dstBuffer_{};

    uint32_t        sampleRate_;
    uint32_t        grainSize_;
    uint32_t        overlap_;

    std::array<kiss_fft_cpx, maxFrameLength>    fftData_{};
    std::array<float, maxFrameLength / 2 + 1>   lastPhase_{};
    std::array<float, maxFrameLength / 2 + 1>   sumPhase_{};
    std::array<float, 2 * maxFrameLength>       outputAccum_{};
    std::array<float, maxFrameLength>           analysisFrq_{};
    std::array<float, maxFrameLength>           analysisMag_{};
    std::array<float, maxFrameLength>           synthesisFrq_{};
    std::array<float, maxFrameLength>           synthesisMag_{};

    fft_wrapper     fft_;

private:
    

    
};
