#include <map>
#include <algorithm>

#include "fft.h"
#include "common.h"

fft_wrapper::fft_wrapper()
{
}

fft_wrapper::fft_wrapper(uint32_t sampleRate, uint32_t windowSize)
{
    init(sampleRate, windowSize);
}

fft_wrapper::~fft_wrapper()
{
    if (cfg_) {
        free(cfg_);
    }
}

void fft_wrapper::init(uint32_t sampleRate, uint32_t windowSize)
{
    sampleRate_ = sampleRate;
    windowSize_ = windowSize;
    in_.resize(windowSize_);
    out_.resize(windowSize_);
    inverse_.resize(windowSize_);
    cfg_ = kiss_fft_alloc(windowSize_, false, NULL, NULL);
    inverse_cfg_ = kiss_fft_alloc(windowSize_, true, NULL, NULL);
}

void
fft_wrapper::setInput(const double * realData, const double * imagData)
{
    for (unsigned i = 0; i < windowSize_; i++) {
        in_[i].r = realData[i];
        in_[i].i = imagData != nullptr ? imagData[i] : 0;
    }
}

void
fft_wrapper::setOutput(const double * realData, const double * imagData)
{
    for (unsigned i = 0; i < windowSize_; i++) {
        out_[i].r = realData[i];
        out_[i].i = imagData != nullptr ? imagData[i] : 0;
    }
}

fftDataBlob
fft_wrapper::computeStft()
{
    fftDataBlob dataBlob;

    kiss_fft(cfg_, &in_[0], &out_[0]);

    auto magWindowSize = windowSize_ / 2;
    double freqPerBin = sampleRate_ / windowSize_;
    double frequency, amplitude, phase;
    double threshold = 0.001; // -60.0; -60dB

    for (uint32_t i = 0; i < windowSize_; i++) {
        amplitude = sqrt(
            static_cast<double>(out_[i].r * out_[i].r) +
            static_cast<double>(out_[i].i * out_[i].i)
        );
        frequency = i * freqPerBin;
        phase = atan2(static_cast<double>(out_[i].i), static_cast<double>(out_[i].r));
        dataBlob.emplace_back(fftData{ frequency, amplitude, phase });
    }

    if (SYNTHTYPE == SynthType::resynth) {
        // limit bands?
        auto bandLimit = 64;
        if (bandLimit) {
            // sort by amplitude first
            std::sort(dataBlob.begin(), dataBlob.end(), [](fftData a, fftData b) {
                return a.amplitude > b.amplitude;
            });
            dataBlob.resize(bandLimit);
        }
    } 

    return dataBlob;
}

std::vector<kiss_fft_cpx>
fft_wrapper::getRawOutput()
{
    return out_;
}

std::vector<double>
fft_wrapper::computeInverseStft()
{
    std::vector<double> dataBlob;
    kiss_fft(inverse_cfg_, &out_[0], &inverse_[0]);
    for (const auto& bin : inverse_) {
        dataBlob.emplace_back(bin.r / (1 * windowSize_));
    }
    return dataBlob;
}
