#pragma once

#include "audiobuffer.h"
#include "ringbuffer.h"

#include <chrono>
#include <mutex>
#include <memory>
#include <vector>
#include <condition_variable>
#include <atomic>

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
    virtual ~AudioLayer();

    virtual void startStream() = 0;
    virtual void stopStream() = 0;

    inline bool isStarted() const {
        return status_ == Status::Started;
    }

    template< class Rep, class Period >
    bool waitForStart(const std::chrono::duration<Rep, Period>& rel_time) const {
        std::unique_lock<std::mutex> lk(mutex_);
        startedCv_.wait_for(lk, rel_time, [&] {return isStarted(); });
        return isStarted();
    }

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

    bool isRingtoneMuted() const {
        return isRingtoneMuted_;
    }

    void muteRingtone(bool muted) {
        isRingtoneMuted_ = muted;
    }

    /**
    * Range: [-1.0, 1.0]
    */
    void setCaptureGain(double gain) {
        captureGain_ = gain;
    }

    double getCaptureGain() const {
        return captureGain_;
    }

    /**
    * Range: [-1.0, 1.0]
    */
    void setPlaybackGain(double gain)
    {
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

protected:
    bool isCaptureMuted_;
    bool isPlaybackMuted_;
    bool isRingtoneMuted_{ false };
    double captureGain_;
    double playbackGain_;

    AudioBuffer playbackBuffer_;
    AudioBuffer playbackResampleBuffer_;

    std::atomic<Status> status_{ Status::Idle };
    mutable std::condition_variable startedCv_;

    AudioFormat audioFormat_;
    AudioFormat audioInputFormat_;

    mutable std::mutex mutex_ { };

    /**
    * Remove audio offset that can be introduced by certain cheap audio device
    */
    //DcBlocker dcblocker_{};

    /**
    * Manage sampling rate conversion
    */
    //std::unique_ptr<Resampler> resampler_;

private:

    /**
    * Time of the last incoming call notification
    */
    std::chrono::system_clock::time_point lastNotificationTime_;
};