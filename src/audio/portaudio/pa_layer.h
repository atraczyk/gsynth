#pragma once

#include <memory>
#include <array>

#include "audiolayer.h"

using freqData = std::vector<std::pair<double, double>>;

class PortAudioLayer: public AudioLayer {

public:
    PortAudioLayer();
    PortAudioLayer(AudioFormat out_format, AudioFormat in_format);
    ~PortAudioLayer();

    virtual void startStream();
    virtual void stopStream();

    freqData getFrequencyData();

private:
    PortAudioLayer(const PortAudioLayer&) = delete;
    PortAudioLayer& operator=(const PortAudioLayer&) = delete;

    struct PortAudioLayerImpl;
    std::unique_ptr<PortAudioLayerImpl> pimpl_;
};
