#include "audiolayer.h"

#include "Log.hpp"

#include <ctime>
#include <algorithm>

AudioLayer::AudioLayer()
    : isCaptureMuted_(false)
    , isPlaybackMuted_(false)
    , captureGain_(1.0)
    , playbackGain_(1.0)
    , audioFormat_(AudioFormat::MONO())
    , audioInputFormat_(AudioFormat::MONO())
{
}

AudioLayer::AudioLayer(AudioFormat out_format, AudioFormat in_format)
    : isCaptureMuted_(false)
    , isPlaybackMuted_(false)
    , captureGain_(1.0)
    , playbackGain_(1.0)
    , audioFormat_(out_format)
    , audioInputFormat_(in_format)
{
}

AudioLayer::~AudioLayer()
{
}

void AudioLayer::flushMain()
{
    std::lock_guard<std::mutex> lock(mutex_);
    // TODO: flush buffer
}
