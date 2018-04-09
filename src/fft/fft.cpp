#include <map>
#include <algorithm>

#include "fft.h"

fft_wrapper::fft_wrapper()
{
}

fft_wrapper::fft_wrapper(uint32_t sampleRate, uint32_t windowSize, bool isInverse)
{
    init(sampleRate, windowSize, isInverse);
}

fft_wrapper::~fft_wrapper()
{
    if (cfg_) {
        free(cfg_);
    }
}

void fft_wrapper::init(uint32_t sampleRate, uint32_t windowSize, bool isInverse)
{
    sampleRate_ = sampleRate;
    windowSize_ = windowSize;
    in_.resize(windowSize_);
    out_.resize(windowSize_);
    mag_.resize(windowSize_);
    cfg_ = kiss_fft_alloc(windowSize_, isInverse, NULL, NULL);
}

void
fft_wrapper::setRealInput(const double * realData)
{
    for (int i = 0; i < windowSize_; i++) {
        in_[i].r = realData[i];
        in_[i].i = 0;
    }
}

std::vector<std::pair<double, double>>
fft_wrapper::computeFrequencies(bool sort, bool useThreshold, int bandLimit)
{
    std::vector<std::pair<double, double>> freqs;

    kiss_fft(cfg_, &in_[0], &out_[0]);

    using magInd = std::pair<double, double>;
    std::vector<magInd> magInds;

    auto magWindowSize = windowSize_ / 2;

    for (uint32_t i = 1; i < windowSize_ / 2; i++) {
        mag_[i] = sqrt((out_[i].r * out_[i].r) + (out_[i].i * out_[i].i)) / static_cast<double>(magWindowSize);
        magInds.emplace_back(std::make_pair(mag_[i], i));
    }

    // sort by magnitude?
    if (sort) {
        std::sort(magInds.begin(), magInds.end(), [](magInd a, magInd b) {
            return a.first > b.first;
        });
    }

    // limit bands?
    if (bandLimit) {
        magInds.resize(bandLimit);
    }

    double threshold = 0.001; // -60.0; -60dB
    auto applyThreshold = useThreshold && bandLimit;
    for (const auto& mi : magInds) {
        auto power = mi.first; // dB: 20 * log10(mi.first);
        if (applyThreshold && (power < threshold)) {
            continue;
        }
        auto weightedIndex = static_cast<double>(mi.second);
        auto currIndex = mi.second;
        if (currIndex - 1 >= 0)
            weightedIndex = weightedIndex - (mag_[currIndex - 1] / mag_[currIndex]);
        if (currIndex + 1 < magWindowSize)
            weightedIndex = weightedIndex + (mag_[currIndex + 1] / mag_[currIndex]);
        auto freq = static_cast<double>(weightedIndex) * sampleRate_ / windowSize_;
        freqs.emplace_back(std::make_pair(freq, power));
    }

    return freqs;
}
