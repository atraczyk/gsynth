#include <map>
#include <algorithm>

#include "fft.h"
#include "common.h"

fft_wrapper::fft_wrapper()
{
}

fft_wrapper::fft_wrapper(uint32_t windowSize)
    : cfg_(nullptr)
    , inverse_cfg_(nullptr)
{
    init(windowSize);
}

fft_wrapper::~fft_wrapper()
{
    if (cfg_) {
        free(cfg_);
    }
    if (inverse_cfg_) {
        free(inverse_cfg_);
    }
}

void
fft_wrapper::init(uint32_t windowSize)
{
    windowSize_ = windowSize;
    if (cfg_) {
        free(cfg_);
    }
    if (inverse_cfg_) {
        free(inverse_cfg_);
    }
    cfg_ = kiss_fft_alloc(windowSize_, false, nullptr, nullptr);
    inverse_cfg_ = kiss_fft_alloc(windowSize_, true, nullptr, nullptr);
}

void
fft_wrapper::computeA(kiss_fft_cpx* data, size_t frameSize)
{
    if (frameSize && cfg_ && frameSize != windowSize_) {
        free(cfg_);
        cfg_ = kiss_fft_alloc(windowSize_, false, nullptr, nullptr);
    }
    kiss_fft(cfg_, &data[0], &data[0]);
}

void
fft_wrapper::computeInverseA(kiss_fft_cpx* data, size_t frameSize)
{
    if (frameSize && inverse_cfg_ && frameSize != windowSize_) {
        free(inverse_cfg_);
        inverse_cfg_ = kiss_fft_alloc(windowSize_, true, nullptr, nullptr);
    }
    kiss_fft(inverse_cfg_, &data[0], &data[0]);
}