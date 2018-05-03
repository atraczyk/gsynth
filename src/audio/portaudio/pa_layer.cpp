#include "pa_layer.h"

#include "ringbuffer.h"
#include "Log.hpp"
#include "filter.h"
#include "soundfile.h"
#include "common.h"

#include <portaudio.h>

#include <algorithm>
#include <cmath>
#include <vector>

#ifndef _WIN32
#include <ncurses.h>
#endif // !_WIN32

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
    bool canPlay_;

    std::vector<std::string> getDeviceByType(bool) const;

    PaDeviceIndex indexIn_;
    PaDeviceIndex indexOut_;

    std::shared_ptr<RingBuffer> inputBuffer_;
    std::shared_ptr<RingBuffer> processedBuffer_;

    fft_wrapper fft_;
    std::thread fftThread_;
    std::mutex fftMutex_;
    std::condition_variable fftCv_;

    fftDataBlob currentFftData_;
    std::vector<double> currentFftInverse_;
    std::mutex fftDataMutex_;

    LowPassFilter lpFilter;

    uint64_t outFrame_;
    uint64_t inFrame_;
    
    double pitchShift_ = 0;

    std::vector<float> out_;

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
{
}

PortAudioLayer::PortAudioLayer(AudioFormat out_format, AudioFormat in_format)
    : AudioLayer(out_format, in_format)
    , pimpl_{ new PortAudioLayerImpl(*this) }
{
}

PortAudioLayer::~PortAudioLayer()
{
    WriteWaveFile("out.wav", pimpl_->out_, 1, audioFormat_.sample_rate, 2);
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

fftDataBlob
PortAudioLayer::getFftData()
{
    fftDataBlob ret;
    {
        std::lock_guard<std::mutex> lock(pimpl_->fftDataMutex_);
        ret = pimpl_->currentFftData_;
    }
    return ret;
}

void 
PortAudioLayer::setPitchShift(double pitchShift)
{
    pimpl_->pitchShift_ = pitchShift;
}

PortAudioLayer::PortAudioLayerImpl::PortAudioLayerImpl(PortAudioLayer& parent)
    : indexIn_(paNoDevice)
    , indexOut_(paNoDevice)
    , inputBuffer_(std::make_shared<RingBuffer>())
    , processedBuffer_(std::make_shared<RingBuffer>())
    , fft_(parent.audioFormat_.sample_rate, parent.audioFormat_.sample_rate * 0.05)
    , canPlay_(false)
    , outFrame_(0)
    , inFrame_(0)
    , lpFilter(parent.audioFormat_.sample_rate, 4000.0, 1.0, 8.0)
{
    init(parent);
}

PortAudioLayer::PortAudioLayerImpl::~PortAudioLayerImpl()
{
    running_ = false;
    fftCv_.notify_all();
    if (fftThread_.joinable()) {
        fftThread_.join();
    }
    terminate();
}

std::vector<std::string>
PortAudioLayer::PortAudioLayerImpl::getDeviceByType(bool playback) const
{
    std::vector<std::string> ret;
    int numDevices = 0;

    numDevices = Pa_GetDeviceCount();
    if (numDevices < 0) {
        DBGOUT("PortAudioLayer error : %s", Pa_GetErrorText(numDevices));
    } else {
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

    auto startOutFrame = outFrame_;
    auto endOutFrame = outFrame_ + framesPerBuffer;
    AudioSample *out = (AudioSample*)outputBuffer;

    auto readyBufSize = processedBuffer_.get()->getBufSize();
    //DBGOUT("readyBufSize: %d", readyBufSize);
    
    // avoid buffer underrun at start
    if (readyBufSize >= framesPerBuffer * 2) {
        canPlay_ = true;
    }

    if (readyBufSize >= framesPerBuffer && canPlay_) {
        auto _fftDataBlob = parentLayer.getFftData();
        auto elapsedTimePerFrame = 1.0 / static_cast<double>(parentLayer.audioFormat_.sample_rate);
        for (; outFrame_ < endOutFrame; outFrame_++) {
            auto t = (outFrame_)* elapsedTimePerFrame;
            AudioSample processed;
            processedBuffer_.get()->try_pop(processed);

            double value = 0.0;
            if (SYNTHTYPE == SynthType::resynth) {
                auto fMult = (1.0 + pitchShift_ / 200);
                for (int i = 0; i < _fftDataBlob.size() && i < 32; i++) {
                    value += 0.5 * _fftDataBlob.at(i).amplitude *
                        sin(M_2_PI * t * _fftDataBlob.at(i).frequency * fMult + _fftDataBlob.at(i).phase);
                }
            }
            if (SYNTHTYPE == SynthType::ifft) {
                auto i = outFrame_ - startOutFrame;
                value = currentFftInverse_.at(outFrame_ - startOutFrame);
                //0.5 * (1.0 + cos(2.0 * pi* x / (256 - 1.0))
                //auto mult = 0.5 - 0.5 * (cos(2.0 * M_PI * k / (thisGrainSize - 1.0)));
                value *= (1.0 + cos(2.0 * M_PI * i / (fft_.getWindowSize() - 1.0)));
            }
            value = value > 1.0 ? 1.0 : (value < -1.0 ? -1.0 : value);
            auto int16Data = static_cast<AudioSample>(value * 32768.0f);
            *out++ = int16Data;
            *out++ = int16Data;

            //out_.emplace_back(processed * 0.000030517578125f);//value);

            //*out++ = processed;
            //*out++ = processed;
        }
    } else {
        DBGOUT("output buffer underrun");
        for (; outFrame_ < endOutFrame; outFrame_++) {
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

    auto endInFrame = inFrame_ + framesPerBuffer;
    AudioSample *in = (AudioSample*)inputBuffer;

    double gain = 1.0f;
    double freq = 100.0f + pitchShift_;
    auto elapsedTimePerFrame = 1.0 / static_cast<double>(parentLayer.audioFormat_.sample_rate);
    auto periodFrames = static_cast<int>((1.0 / freq) / elapsedTimePerFrame);

    for (; inFrame_ < endInFrame; inFrame_++) {
        auto t = (inFrame_) * elapsedTimePerFrame;
        auto value = sin(M_2_PI * t * freq);
        //auto value = 0.5 * sin(M_2_PI * t * 354) + 0.5 * cos(M_2_PI * t * 649);
        auto int16Data = static_cast<AudioSample>(value * gain * 32768.0f);

        // generate sine test
        //inputBuffer_.get()->tryPush(int16Data);

        // mic test
        inputBuffer_.get()->try_push(*in++);
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
    PaDeviceIndex device, Direction direction, uint32_t bufferLength,
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
        bufferLength,//paFramesPerBufferUnspecified,
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
    DBGOUT("Open PortAudio Output Stream");
    if (indexOut_ != paNoDevice) {
        parent.audioFormat_ =  openStreamDevice(
            parent,
            &streams_[Direction::Output],
            indexOut_,
            Direction::Output,
            fft_.getWindowSize(),
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
            fft_.getWindowSize(),
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
        fftCv_.wait(lk, [this] {
            return (inputBuffer_.get()->getBufSize() >= fft_.getWindowSize()) || !running_;
        });
        if (!running_) {
            continue;
        }
        AudioSample sample;
        std::vector<double> data;
        double value;
        auto fftWindowSize = fft_.getWindowSize();
        unsigned i = 0;
        while (data.size() <= fftWindowSize) {
            if (!inputBuffer_.get()->try_pop(sample))
                continue;
            processedBuffer_.get()->try_push(sample);
            value = lpFilter.processSample(sample * 0.000030517578125f);
            double windowMult;
            windowMult = 0.5 * (1.0 - cos(2.0 * M_PI * i / (fftWindowSize - 1.0)));
            data.emplace_back(static_cast<double>(windowMult * value));
            ++i;
        }
        fft_.setData(FFTDirection::In, data);
        {
            std::lock_guard<std::mutex> lk(fftDataMutex_);
            currentFftData_ = fft_.computeStft();
            if (SYNTHTYPE == SynthType::ifft) {
                currentFftInverse_ = fft_.computeInverseStft();
            }
        }
    }
}
