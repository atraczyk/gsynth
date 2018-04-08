#include "pa_layer.h"

#include "ringbuffer.h"
#include "fft.h"
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

    bool running_ = true;

    std::vector<std::string> getDeviceByType(bool) const;

    PaDeviceIndex indexIn_;
    PaDeviceIndex indexOut_;

    AudioBuffer playbackBuff_;

    std::shared_ptr<RingBuffer> primaryRingBuffer_;
    std::shared_ptr<RingBuffer> secondaryRingBuffer_;

    bool canPlay_;

    fft_wrapper fft_;
    std::thread fftThread_;
    std::mutex fftMutex_;
    std::condition_variable fftCv_;

    std::array<PaStream*, static_cast<int>(Direction::End)> streams_;

    void processFFT();

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

PortAudioLayer::PortAudioLayer(AudioFormat out_format, AudioFormat in_format)
    : AudioLayer(out_format, in_format)
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
    , primaryRingBuffer_(std::make_shared<RingBuffer>())
    , secondaryRingBuffer_(std::make_shared<RingBuffer>())
    , fft_(parent.audioFormat_.sample_rate, parent.audioFormat_.sample_rate * 0.025, false)
    , canPlay_(false)
{
    init(parent);
}

PortAudioLayer::PortAudioLayerImpl::~PortAudioLayerImpl()
{
    terminate();
    running_ = false;
    fftCv_.notify_all();
    fftThread_.join();
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

    auto endOutFrame = parentLayer.outFrame_ + framesPerBuffer;
    AudioSample *out = (AudioSample*)outputBuffer;

    // test
    static double gain = 0.1f;
    static double freq = 440.0f;
    auto elapsedTimePerFrame = 1.0 / static_cast<double>(parentLayer.audioFormat_.sample_rate);
    auto periodFrames = static_cast<int>((1.0 / freq) / elapsedTimePerFrame);

    auto readyBufSize = secondaryRingBuffer_.get()->getBufSize();
    DBGOUT("readyBufSize: %d", readyBufSize);
    // avoid buffer underrun at start
    if (readyBufSize >= framesPerBuffer * 2) {
        canPlay_ = true;
    }
    if (readyBufSize >= framesPerBuffer && canPlay_) {
        DBGOUT("playing buffer");
        for (; parentLayer.outFrame_ < endOutFrame; parentLayer.outFrame_++) {
            AudioSample processed;
            //primaryRingBuffer_.get()->tryPop(&processed);
            secondaryRingBuffer_.get()->tryPop(&processed);
            *out++ = processed;
            *out++ = processed;
        }
    } else {
        DBGOUT("output buffer underrun");
        for (; parentLayer.outFrame_ < endOutFrame; parentLayer.outFrame_++) {
            *out++ = 0;
            *out++ = 0;
        }
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

    auto endInFrame = parentLayer.inFrame_ + framesPerBuffer;
    AudioSample *in = (AudioSample*)inputBuffer;

    static double gain = 0.5f;
    static double freq = 300.0f;
    auto elapsedTimePerFrame = 1.0 / static_cast<double>(parentLayer.audioFormat_.sample_rate);
    auto periodFrames = static_cast<int>((1.0 / freq) / elapsedTimePerFrame);

    for (; parentLayer.inFrame_ < endInFrame; parentLayer.inFrame_++) {
        auto t = (parentLayer.inFrame_) * elapsedTimePerFrame;
        //auto value = sin(2 * M_PI * t * freq);
        auto value = 0.75 * sin(2 * M_PI * t * 354) + 0.75 * cos(2 * M_PI * t * 649);
        auto int16Data = static_cast<AudioSample>(value * gain * 32768.0f);

        // generate sine test
        //primaryRingBuffer_.get()->tryPush(int16Data);

        // mic test
        primaryRingBuffer_.get()->tryPush(*in++);
    }
    fftCv_.notify_all();

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
            parent.audioFormat_.nb_channels = std::min( static_cast<unsigned>(outputDeviceInfo->maxOutputChannels),
                                                        parent.audioFormat_.nb_channels);
            parent.audioFormat_.sample_rate = std::min( static_cast<unsigned>(outputDeviceInfo->defaultSampleRate),
                                                        parent.audioFormat_.sample_rate);
        }
        else {
            indexOut_ = paNoDevice;
        }
    }

    if (indexIn_ != paNoDevice) {
        if (const auto inputDeviceInfo = Pa_GetDeviceInfo(indexIn_)) {
            parent.audioInputFormat_.nb_channels = std::min( static_cast<unsigned>(inputDeviceInfo->maxInputChannels),
                                                             parent.audioInputFormat_.nb_channels);
            parent.audioInputFormat_.sample_rate = std::min( static_cast<unsigned>(inputDeviceInfo->defaultSampleRate),
                                                             parent.audioInputFormat_.sample_rate);
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
openStreamDevice(PortAudioLayer& parent, PaStream** stream,
    PaDeviceIndex device, Direction direction,
    PaStreamCallback* callback, void* user_data)
{
    auto is_out = direction == Direction::Output;
    auto device_info = Pa_GetDeviceInfo(device);

    auto rate = is_out ? parent.getFormat().sample_rate : parent.getInputFormat().sample_rate;
    auto channels = is_out ? parent.getFormat().nb_channels : parent.getInputFormat().nb_channels;

    PaStreamParameters params;
    params.device = device;
    params.channelCount = channels;
    params.sampleFormat = paInt16;
    params.suggestedLatency = is_out ? device_info->defaultLowOutputLatency : device_info->defaultLowInputLatency;
    params.hostApiSpecificStreamInfo = nullptr;

    auto err = Pa_OpenStream(
        stream,
        is_out ? nullptr : &params,
        is_out ? &params : nullptr,
        rate,
        paFramesPerBufferUnspecified,
        paNoFlag,
        callback,
        user_data);

    if (err != paNoError)
        DBGOUT("PortAudioLayer error : %s", Pa_GetErrorText(err));

    DBGOUT("stream %d initialized with %d channels @ %dHz", is_out, channels, rate);
    return AudioFormat {   static_cast<unsigned>(rate),
                           static_cast<unsigned>(channels) } ;
}

void
PortAudioLayer::PortAudioLayerImpl::initStream(PortAudioLayer& parent)
{
    //parent.dcblocker_.reset();

    DBGOUT("Open PortAudio Output Stream");
    if (indexOut_ != paNoDevice) {
        parent.audioFormat_ =  openStreamDevice(
            parent,
            &streams_[Direction::Output],
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
            parent,
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

    DBGOUT("Start FFT Handler");
    fftThread_ = std::thread(&PortAudioLayerImpl::processFFT, this);
}

void
PortAudioLayer::PortAudioLayerImpl::processFFT()
{
    while (running_) {
        std::unique_lock<std::mutex> lk(fftMutex_);
        //DBGOUT("processFFT - waiting for buffer");
        fftCv_.wait(lk, [this] {
            //DBGOUT("getBufSize() %d, getWindowSize() %d", primaryRingBuffer_.get()->getBufSize(), fft_.getWindowSize());
            return (primaryRingBuffer_.get()->getBufSize() >= fft_.getWindowSize()) || !running_;
        });
        //DBGOUT("processFFT - processing");
        AudioSample processed;
        auto bufSize = primaryRingBuffer_.get()->getBufSize();
        std::vector<double> data;
        auto fftWindowSize = fft_.getWindowSize();
        //DBGOUT("bufSize: %d, framesPerBuffer %d, fftWindowSize %d", bufSize, framesPerBuffer, fftWindowSize);
        if (bufSize >= fftWindowSize) {
            while (data.size() <= fftWindowSize) {
                if (!primaryRingBuffer_.get()->tryPop(&processed))
                    continue;
                secondaryRingBuffer_.get()->tryPush(processed);
                auto floatData = static_cast<double>(processed / 32768.0f);
                data.emplace_back(floatData);
            }
            fft_.setRealInput(&data[0]);
            auto freqs = fft_.computeFrequencies();

            for (int i = 0; i < freqs.size(); i++) {
                if (i == 0) {
                    DBGOUT("FREQ: i: %d, f: %0.4f, p: %0.4f ", i, freqs.at(i).first, freqs.at(i).second);
                }
                else if (i < 2) {
                    DBGOUT("freq: i: %d, f: %0.4f, p: %0.4f ", i, freqs.at(i).first, freqs.at(i).second);
                }
            }
        }
    }
}