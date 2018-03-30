#include "pa_layer.h"

#include "ringbuffer.h"
#include "Log.hpp"

#include <portaudio.h>

#include <algorithm>
#include <cmath>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265
#endif

enum Direction {
    Input = 0,
    Output = 1,
    End = 2
};

struct PortAudioLayer::PortAudioLayerImpl
{
    PortAudioLayerImpl(PortAudioLayer&);
    ~PortAudioLayerImpl();

    void init(PortAudioLayer&);
    void terminate() const;
    void initStream(PortAudioLayer&);

    std::vector<std::string> getDeviceByType(bool) const;

    PaDeviceIndex indexIn_;
    PaDeviceIndex indexOut_;

    AudioBuffer playbackBuff_;

    std::shared_ptr<RingBuffer> mainRingBuffer_;

    std::array<PaStream*, static_cast<int>(Direction::End)> streams_;

    int paOutputCallback(PortAudioLayer& userData,
                         const AudioSample* inputBuffer,
                         AudioSample* outputBuffer,
                         unsigned long framesPerBuffer,
                         const PaStreamCallbackTimeInfo* timeInfo,
                         PaStreamCallbackFlags statusFlags);

    int paInputCallback(PortAudioLayer& userData,
                        const AudioSample* inputBuffer,
                        AudioSample* outputBuffer,
                        unsigned long framesPerBuffer,
                        const PaStreamCallbackTimeInfo* timeInfo,
                        PaStreamCallbackFlags statusFlags);
};

PortAudioLayer::PortAudioLayer()
    : AudioLayer()
    , pimpl_{ new PortAudioLayerImpl(*this) }
    , outFrame_(0)
    , inFrame_(0)
{
}

PortAudioLayer::~PortAudioLayer()
{
}

void
PortAudioLayer::startStream()
{
    {
        std::lock_guard<std::mutex> lock(mutex_);
        if (status_ != Status::Idle)
            return;
        status_ = Status::Started;
    }
    pimpl_->initStream(*this);
}

void
PortAudioLayer::stopStream()
{
    if (status_ != Status::Started)
        return;

    DBGOUT("Stop PortAudio Streams");

    for (auto& st_ptr : pimpl_->streams_) {
        auto err = Pa_StopStream(st_ptr);
        if (err != paNoError)
            DBGOUT("Pa_StopStream error : %s", Pa_GetErrorText(err));

        err = Pa_CloseStream(st_ptr);
        if (err != paNoError)
            DBGOUT("Pa_StopStream error : %s", Pa_GetErrorText(err));
    }

    {
        std::lock_guard<std::mutex> lock{ mutex_ };
        status_ = Status::Idle;
    }
}

PortAudioLayer::PortAudioLayerImpl::PortAudioLayerImpl(PortAudioLayer& parent)
    : playbackBuff_{ 0, parent.audioFormat_ }
    , indexIn_(paNoDevice)
    , indexOut_(paNoDevice)
{
    init(parent);
}

PortAudioLayer::PortAudioLayerImpl::~PortAudioLayerImpl()
{
    terminate();
}

std::vector<std::string>
PortAudioLayer::PortAudioLayerImpl::getDeviceByType(bool playback) const
{
    std::vector<std::string> ret;
    int numDevices = 0;

    numDevices = Pa_GetDeviceCount();
    if (numDevices < 0)
        DBGOUT("PortAudioLayer error : %s", Pa_GetErrorText(numDevices));
    else {
        for (int i = 0; i < numDevices; i++) {
            const auto deviceInfo = Pa_GetDeviceInfo(i);
            if (playback) {
                if (deviceInfo->maxOutputChannels > 0)
                    ret.push_back(deviceInfo->name);
            }
            else {
                if (deviceInfo->maxInputChannels > 0)
                    ret.push_back(deviceInfo->name);
            }
        }
    }
    return ret;
}

float getElapsedTime()
{
    using namespace std::chrono;
    static high_resolution_clock::time_point start = high_resolution_clock::now();
    return duration_cast<duration<float>>(high_resolution_clock::now() - start).count();
}

double getFloatPrecision(const double& value, const double& precision)
{
    return (floor((value * pow(10, precision) + 0.5)) / pow(10, precision));
}

int
PortAudioLayer::PortAudioLayerImpl::paOutputCallback(   PortAudioLayer& parentLayer,
                                                        const AudioSample* inputBuffer,
                                                        AudioSample* outputBuffer,
                                                        unsigned long framesPerBuffer,
                                                        const PaStreamCallbackTimeInfo* timeInfo,
                                                        PaStreamCallbackFlags statusFlags)
{
    // unused arguments
    (void)inputBuffer;
    (void)timeInfo;
    (void)statusFlags;
    (void)inputBuffer;

    if (framesPerBuffer == 0) {
        DBGOUT("No frames for output.");
        return paContinue;
    }

    auto startOutFrame = parentLayer.outFrame_;
    auto endOutFrame = startOutFrame + framesPerBuffer;
    AudioSample *out = (AudioSample*)outputBuffer;

    // test
    static double gain = 0.2f;
    static double freq = 440.0f;
    auto elapsedTimePerFrame = 1.0 / static_cast<double>(parentLayer.audioFormat_.sample_rate);

    for (; parentLayer.outFrame_ < endOutFrame; parentLayer.outFrame_++) {
        auto t = parentLayer.outFrame_ * elapsedTimePerFrame;

        auto value = sin(2 * M_PI * t * freq );
        auto int16Data = static_cast<AudioSample>(value * gain * 32768.0f);
        *out++ = int16Data;
        *out++ = int16Data;
    }

    return paContinue;
}

int
PortAudioLayer::PortAudioLayerImpl::paInputCallback(    PortAudioLayer& parentLayer,
                                                        const AudioSample* inputBuffer,
                                                        AudioSample* outputBuffer,
                                                        unsigned long framesPerBuffer,
                                                        const PaStreamCallbackTimeInfo* timeInfo,
                                                        PaStreamCallbackFlags statusFlags)
{
    // unused arguments
    (void)outputBuffer;
    (void)timeInfo;
    (void)statusFlags;

    if (framesPerBuffer == 0) {
        DBGOUT("No frames for input.");
        return paContinue;
    }

    // TODO: get frames from input

    return paContinue;
}

void
PortAudioLayer::PortAudioLayerImpl::init(PortAudioLayer& parent)
{
    DBGOUT("Init PortAudioLayer");
    const auto err = Pa_Initialize();
    if (err != paNoError) {
        DBGOUT("PortAudioLayer error : %s", Pa_GetErrorText(err));
        terminate();
    }

    indexOut_ = indexOut_ == paNoDevice ? Pa_GetDefaultOutputDevice() : indexOut_;
    indexIn_ = indexIn_ == paNoDevice ? Pa_GetDefaultInputDevice() : indexIn_;

    if (indexOut_ != paNoDevice) {
        if (const auto outputDeviceInfo = Pa_GetDeviceInfo(indexOut_)) {
            parent.audioFormat_.nb_channels = outputDeviceInfo->maxOutputChannels;
            parent.audioFormat_.sample_rate = outputDeviceInfo->defaultSampleRate;
        }
        else {
            indexOut_ = paNoDevice;
        }
    }

    if (indexIn_ != paNoDevice) {
        if (const auto inputDeviceInfo = Pa_GetDeviceInfo(indexIn_)) {
            parent.audioInputFormat_.nb_channels = inputDeviceInfo->maxInputChannels;
            parent.audioInputFormat_.sample_rate = inputDeviceInfo->defaultSampleRate;
        }
        else {
            indexIn_ = paNoDevice;
        }
    }

    std::fill(std::begin(streams_), std::end(streams_), nullptr);
}

void
PortAudioLayer::PortAudioLayerImpl::terminate() const
{
    DBGOUT("PortAudioLayer terminate.");
    auto err = Pa_Terminate();
    if (err != paNoError)
        DBGOUT("PortAudioLayer error : %s", Pa_GetErrorText(err));
}

static AudioFormat
openStreamDevice(PaStream** stream,
    PaDeviceIndex device, Direction direction,
    PaStreamCallback* callback, void* user_data)
{
    auto is_out = direction == Direction::Output;
    auto device_info = Pa_GetDeviceInfo(device);

    PaStreamParameters params;
    params.device = device;
    params.channelCount = is_out ? 2 : device_info->maxInputChannels;
    params.sampleFormat = paInt16;
    params.suggestedLatency = is_out ? device_info->defaultLowOutputLatency : device_info->defaultLowInputLatency;
    params.hostApiSpecificStreamInfo = nullptr;

    auto outRate = 44100;

    auto err = Pa_OpenStream(
        stream,
        is_out ? nullptr : &params,
        is_out ? &params : nullptr,
        is_out ? outRate : device_info->defaultSampleRate,
        paFramesPerBufferUnspecified,
        paNoFlag,
        callback,
        user_data);

    if (err != paNoError)
        DBGOUT("PortAudioLayer error : %s", Pa_GetErrorText(err));

    DBGOUT("stream %d initialized with %d channels @ %.0fHz", is_out, params.channelCount, device_info->defaultSampleRate);
     return AudioFormat {   static_cast<unsigned>(is_out ? outRate : device_info->defaultSampleRate),
                            static_cast<unsigned>(params.channelCount) } ;
}

void
PortAudioLayer::PortAudioLayerImpl::initStream(PortAudioLayer& parent)
{
    //parent.dcblocker_.reset();

    DBGOUT("Open PortAudio Output Stream");
    if (indexOut_ != paNoDevice) {
        parent.audioFormat_ =  openStreamDevice(&streams_[Direction::Output],
            indexOut_,
            Direction::Output,
            [] (const void* inputBuffer,
            void* outputBuffer,
            unsigned long framesPerBuffer,
            const PaStreamCallbackTimeInfo* timeInfo,
            PaStreamCallbackFlags statusFlags,
            void* userData) -> int {
                auto layer = static_cast<PortAudioLayer*>(userData);
                return layer->pimpl_->paOutputCallback( *layer,
                                                        static_cast<const AudioSample*>(inputBuffer),
                                                        static_cast<AudioSample*>(outputBuffer),
                                                        framesPerBuffer,
                                                        timeInfo,
                                                        statusFlags);
            },
            &parent);
    }
    else {
        DBGOUT("Error: No valid output device. There will be no sound.");
    }

    DBGOUT("Open PortAudio Input Stream");
    if (indexIn_ != paNoDevice) {
        parent.audioInputFormat_ = openStreamDevice(
            &streams_[Direction::Input],
            indexIn_,
            Direction::Input,
            [] (const void* inputBuffer,
            void* outputBuffer,
            unsigned long framesPerBuffer,
            const PaStreamCallbackTimeInfo* timeInfo,
            PaStreamCallbackFlags statusFlags,
            void* userData) -> int {
                auto layer = static_cast<PortAudioLayer*>(userData);
                return layer->pimpl_->paInputCallback(  *layer,
                                                        static_cast<const AudioSample*>(inputBuffer),
                                                        static_cast<AudioSample*>(outputBuffer),
                                                        framesPerBuffer,
                                                        timeInfo,
                                                        statusFlags);
            },
            &parent);
    }
    else {
        DBGOUT("Error: No valid input device. There will be no mic.");
    }

    DBGOUT("Start PortAudio Streams");
    for (auto& st_ptr : streams_) {
        if (st_ptr) {
            auto err = Pa_StartStream(st_ptr);
            if (err != paNoError)
                DBGOUT("PortAudioLayer error : %s", Pa_GetErrorText(err));
        }
    }
}