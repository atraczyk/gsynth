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
fft_wrapper::setRealInput(const double * realData)
{
    for (unsigned i = 0; i < windowSize_; i++) {
        in_[i].r = realData[i];
        in_[i].i = 0;
    }
}

fftDataBlob
fft_wrapper::computeStft(bool useThreshold)
{
    fftDataBlob dataBlob;

    kiss_fft(cfg_, &in_[0], &out_[0]);

    using magInd = std::pair<double, double>;
    std::vector<magInd> magInds;

    auto magWindowSize = windowSize_ / 2;
    double freqPerBin = sampleRate_ / windowSize_;
    double frequency, amplitude, phase;
    double threshold = 0.001; // -60.0; -60dB

    for (uint32_t i = 0; i < magWindowSize; i++) {
        amplitude = sqrt((out_[i].r * out_[i].r) + (out_[i].i * out_[i].i)) /
            static_cast<double>(magWindowSize); // dB: 20 * log10(amplitude);

        if (useThreshold && (amplitude < threshold)) {
            continue;
        }
        
        frequency = i * freqPerBin;

        if (out_[i].i == 0.0) {
            phase = 0.0;
        }
        else if (out_[i].r == 0.0) {
            phase = out_[i].i > 0.0 ? M_PI_2 : -M_PI_2;
        }
        else {
            phase = atan2(out_[i].i, out_[i].r);
        }

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

std::vector<double>
fft_wrapper::computeInverseStft(double pitchShift)
{
    auto shift = static_cast<int>(fabs(pitchShift / 10));
    if (shift && shift < 20) {
        auto mid = out_.end() - out_.size() / 2;
        auto rmid = out_.rend() - out_.size() / 2;

        std::rotate(out_.begin(), out_.begin() + shift, mid);
        std::rotate(out_.rbegin(), out_.rbegin() + shift, rmid);
        for (int i = 0; i < shift; i++) {
            *(mid - i) = { 0 ,0 };
            *(rmid - i) = { 0 ,0 };
        }
    }

    std::vector<double> dataBlob;
    kiss_fft(inverse_cfg_, &out_[0], &inverse_[0]);
    for (const auto& bin : inverse_) {
        dataBlob.emplace_back(bin.r / (1 * windowSize_));
    }
    return dataBlob;
}

