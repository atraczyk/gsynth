#pragma once

#include <memory>
#include <array>

#include "audiolayer.h"

class PortAudioLayer: public AudioLayer {

public:
    PortAudioLayer();
    PortAudioLayer(AudioFormat out_format, AudioFormat in_format);
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
