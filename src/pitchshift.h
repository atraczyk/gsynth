#pragma once

#include <vector>

#include "ringbuffer.h"
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

private:
    std::array<float, maxFrameLength>           srcRingBuffer_{};
    std::array<float, maxFrameLength>           dstRingBuffer_{};
    std::array<kiss_fft_cpx, maxFrameLength>    fftData_{};
    std::array<float, maxFrameLength / 2 + 1>   lastPhase_{};
    std::array<float, maxFrameLength / 2 + 1>   sumPhase_{};
    std::array<float, 2 * maxFrameLength>       outputAccum_{};
    std::array<float, maxFrameLength>           analysisFrq_{};
    std::array<float, maxFrameLength>           analysisMag_{};
    std::array<float, maxFrameLength>           synthesisFrq_{};
    std::array<float, maxFrameLength>           synthesisMag_{};

    lock_free_queue<float>                      srcRingBuffer0_{};
    lock_free_queue<float>                      dstRingBuffer0_{};

    uint32_t        sampleRate_;
    uint32_t        grainSize_;
    uint32_t        overlap_;
    fft_wrapper     fft_;
};
