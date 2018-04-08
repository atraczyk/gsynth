#pragma once

#include "ringbuffer.h"

#include <chrono>
#include <mutex>
#include <memory>
#include <vector>
#include <condition_variable>
#include <atomic>

using AudioSample = int16_t;
static const constexpr unsigned DEFAULT_SAMPLE_RATE = 44100;

struct AudioFormat {
    unsigned sample_rate;
    unsigned nb_channels;

    constexpr AudioFormat(unsigned sr, unsigned c) : sample_rate(sr), nb_channels(c) {}

    inline bool operator == (const AudioFormat &b) const
    {
        return ((b.sample_rate == sample_rate) && (b.nb_channels == nb_channels));
    }

    inline bool operator != (const AudioFormat &b) const
    {
        return !(*this == b);
    }

    static const AudioFormat NONE() { return AudioFormat{ 0, 0 }; }
    static const AudioFormat MONO() { return AudioFormat{ DEFAULT_SAMPLE_RATE, 1 }; }
    static const AudioFormat STEREO() { return AudioFormat{ DEFAULT_SAMPLE_RATE, 2 }; }
};

enum class DeviceType {
    PLAYBACK,
    CAPTURE
};

class AudioLayer {

private:
    AudioLayer(const AudioLayer&) = delete;
    AudioLayer& operator=(const AudioLayer&) = delete;

protected:
    enum class Status {
        Idle,
        Starting,
        Started
    };

public:
    AudioLayer();
    AudioLayer(AudioFormat out_format, AudioFormat in_format);
    virtual ~AudioLayer();

    virtual void startStream() = 0;
    virtual void stopStream() = 0;

    void flushMain();

    bool isCaptureMuted() const {
        return isCaptureMuted_;
    }

    void muteCapture(bool muted) {
        isCaptureMuted_ = muted;
    }

    bool isPlaybackMuted() const {
        return isPlaybackMuted_;
    }

    void mutePlayback(bool muted) {
        isPlaybackMuted_ = muted;
    }

    void setCaptureGain(double gain) {
        captureGain_ = gain;
    }

    double getCaptureGain() const {
        return captureGain_;
    }

    void setPlaybackGain(double gain) {
        playbackGain_ = gain;
    }

    double getPlaybackGain() const {
        return playbackGain_;
    }

    unsigned int getSampleRate() const {
        return audioFormat_.sample_rate;
    }

    AudioFormat getFormat() const {
        return audioFormat_;
    }

    AudioFormat getInputFormat() const {
        return audioInputFormat_;
    }

protected:
    bool isCaptureMuted_{ false };
    bool isPlaybackMuted_{ false };
    double captureGain_;
    double playbackGain_;

    std::atomic<Status> status_{ Status::Idle };

    AudioFormat audioFormat_;
    AudioFormat audioInputFormat_;

    mutable std::mutex mutex_ { };

};