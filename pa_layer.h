#pragma once

#include "audiolayer.h"

#include <memory>
#include <array>

class PortAudioLayer: AudioLayer {

public:
    PortAudioLayer();
    ~PortAudioLayer();

    virtual void startStream();

    virtual void stopStream();

    uint64_t outFrame_;
    uint64_t inFrame_;

private:
    PortAudioLayer(const PortAudioLayer&) = delete;
    PortAudioLayer& operator=(const PortAudioLayer&) = delete;

    struct PortAudioLayerImpl;
    std::unique_ptr<PortAudioLayerImpl> pimpl_;
};
