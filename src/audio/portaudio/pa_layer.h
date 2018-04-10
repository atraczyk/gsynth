#pragma once

#include <memory>
#include <array>

#include "audiolayer.h"
#include "fft.h"

class PortAudioLayer: public AudioLayer {

public:
    PortAudioLayer();
    PortAudioLayer(AudioFormat out_format, AudioFormat in_format);
    ~PortAudioLayer();

    virtual void startStream();
    virtual void stopStream();

    fftDataBlob getFftData();

private:
    PortAudioLayer(const PortAudioLayer&) = delete;
    PortAudioLayer& operator=(const PortAudioLayer&) = delete;

    struct PortAudioLayerImpl;
    std::unique_ptr<PortAudioLayerImpl> pimpl_;
};
