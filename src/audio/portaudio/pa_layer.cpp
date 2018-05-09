#include "pa_layer.h"

#include "ringbuffer.h"
#include "Log.hpp"
#include "filter.h"
#include "soundfile.h"
#include "pitchshift.h"
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

using ringbuffer = lock_free_queue<AudioSample>;

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

    std::shared_ptr<ringbuffer> inputBuffer_;
    std::shared_ptr<ringbuffer> processedBuffer_;

    fft_wrapper fft_;
    std::thread fftThread_;
    std::mutex fftMutex_;
    std::condition_variable fftCv_;

    fftDataBlob currentFftData_;
    std::vector<double> currentFftInverse_;
    std::mutex fftDataMutex_;

    LowPassFilter lpFilter;

    PitchShifter pitchShifter;

    uint64_t outFrame_;
    uint64_t inFrame_;
    uint64_t inBufPos = 0;
    
    float pitchShift_ = 1.0;

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
    , inputBuffer_(std::make_shared<ringbuffer>())
    , processedBuffer_(std::make_shared<ringbuffer>())
    , fft_(1024)
    , pitchShifter(parent.audioFormat_.sample_rate, 1024, 8.0)
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
            pitchShifter.grainSize_,
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
            pitchShifter.grainSize_,
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

    //DBGOUT("Start FFT Handler");
    //fftThread_ = std::thread(&PortAudioLayerImpl::processFFT, this);
}



void
PortAudioLayer::PortAudioLayerImpl::processFFT()
{
    while (running_) {
        std::unique_lock<std::mutex> lk(fftMutex_);
        fftCv_.wait(lk, [this] {
            return (inBufPos >= fft_.getWindowSize()) || !running_;
        });
        if (!running_) {
            continue;
        }

        auto grainSize = fft_.getWindowSize();

        static float magnitude;
        static float phase;
        static float tmp;
        static float window;
        static float real = 0;
        static float imag = 0;
        static float freqPerBin;
        static float expectedFreq;
        static long qpd;
        static long index;
        static long latency;
        static long stepSize;
        static long fftFrameSize2;
        static long overlap;
        static int i, k;
        if (inBufPos == 0) {
            fftFrameSize2 = grainSize / 2;
            stepSize = grainSize / overlap;
            freqPerBin = SAMPLERATE / (double)grainSize;
            expectedFreq = M_2_PI * (double)stepSize / (double)grainSize;
            latency = grainSize - stepSize;
            inBufPos = latency;
        }

        AudioSample sample;
        while (i <= grainSize) {
            if (!inputBuffer_.get()->try_pop(sample))
                continue;
            pitchShifter.srcBuffer_[inBufPos - latency] = sample * INT16TOFLOAT;
            processedBuffer_.get()->try_push(
                pitchShifter.dstBuffer_[inBufPos - latency]
            );
            inBufPos++;
        }
    }
}

int
PortAudioLayer::PortAudioLayerImpl::paOutputCallback(PortAudioLayer& parentLayer,
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

    auto pitchShiftRatio = pitchShift_;

    auto elapsedTimePerFrame = 1.0 / static_cast<double>(parentLayer.audioFormat_.sample_rate);
    for (; outFrame_ < endOutFrame; outFrame_++) {
        auto t = (outFrame_)* elapsedTimePerFrame;

        auto grainSize = pitchShifter.grainSize_;

        static float magnitude;
        static float phase;
        static float tmp;
        static float window;
        static float real = 0;
        static float imag = 0;
        static float freqPerBin;
        static float expectedFreq;
        static long qpd;
        static long index;
        static long latency;
        static long stepSize;
        static long fftFrameSize2;
        static long overlap = 8;
        static int i, k;
        if (inBufPos == 0) {
            fftFrameSize2 = grainSize / 2;
            stepSize = grainSize / overlap;
            freqPerBin = SAMPLERATE / (double)grainSize;
            expectedFreq = M_2_PI * (double)stepSize / (double)grainSize;
            latency = grainSize - stepSize;
            inBufPos = latency;
        }

        AudioSample sample;
        inputBuffer_.get()->try_pop(sample);
        pitchShifter.srcBuffer_[inBufPos] = sample * INT16TOFLOAT;
        auto processed = pitchShifter.dstBuffer_[inBufPos - latency];
        auto int16Data = static_cast<AudioSample>(processed * FLOATTOINT16);
        *out++ = int16Data;
        *out++ = int16Data;
        inBufPos++;

        if (inBufPos >= grainSize) {
            inBufPos = latency;

            for (k = 0; k < grainSize; k++) {
                window = -.5*cos(M_2_PI * (double)k / (double)grainSize) + .5;
                FFT_R(pitchShifter.fftData_, k) = pitchShifter.srcBuffer_[k] * window;
                FFT_I(pitchShifter.fftData_, k) = 0.;
            }

            fft_.computeA(&pitchShifter.fftData_[0]);
            for (k = 0; k <= fftFrameSize2; k++) {
                real = FFT_R(pitchShifter.fftData_, k);
                imag = FFT_I(pitchShifter.fftData_, k);

                magnitude = 2. * sqrt(real*real + imag*imag);
                phase = atan2(imag, real);

                /* compute phase difference */
                tmp = phase - pitchShifter.lastPhase_[k];
                pitchShifter.lastPhase_[k] = phase;

                /* subtract expected phase difference */
                tmp -= (double)k * expectedFreq;

                /* map delta phase into +/- Pi interval */
                qpd = tmp / M_PI;
                if (qpd >= 0) qpd += qpd & 1;
                else qpd -= qpd & 1;
                tmp -= M_PI * (double)qpd;

                /* get deviation from bin frequency from the +/- Pi interval */
                tmp = overlap * tmp / (M_2_PI);

                /* compute the k-th partials' true frequency */
                tmp = (double)k*freqPerBin + tmp*freqPerBin;

                /* store magnitude and true frequency in analysis arrays */
                pitchShifter.analysisMag_[k] = magnitude;
                pitchShifter.analysisFrq_[k] = tmp;
            }

            memset(&pitchShifter.synthesisMag_, 0, grainSize * sizeof(float));
            memset(&pitchShifter.synthesisFrq_, 0, grainSize * sizeof(float));
            for (k = 0; k <= fftFrameSize2; k++) {
                index = k*pitchShiftRatio;
                if (index <= fftFrameSize2) {
                    pitchShifter.synthesisMag_[index] = pitchShifter.analysisMag_[k];
                    pitchShifter.synthesisFrq_[index] = pitchShifter.analysisFrq_[k] * pitchShiftRatio;
                }
            }

            for (k = 0; k <= fftFrameSize2; k++) {
                /* get magnitude and true frequency from synthesis arrays */
                magnitude = pitchShifter.synthesisMag_[k];
                tmp = pitchShifter.synthesisFrq_[k];

                /* subtract bin mid frequency */
                tmp -= (double)k*freqPerBin;

                /* get bin deviation from freq deviation */
                tmp /= freqPerBin;

                /* take overlappingSamples into account */
                tmp = M_2_PI * tmp / overlap;

                /* add the overlap phase advance back in */
                tmp += (double)k*expectedFreq;

                /* accumulate delta phase to get bin phase */
                pitchShifter.sumPhase_[k] += tmp;
                phase = pitchShifter.sumPhase_[k];

                /* get real and imag part and re-interleave */
                FFT_R(pitchShifter.fftData_, k) = magnitude*cos(phase);
                FFT_I(pitchShifter.fftData_, k) = magnitude*sin(phase);
            }

            /* zero negative frequencies */
            for (k = grainSize / 2; k < grainSize; k++) {
                FFT_R(pitchShifter.fftData_, k) = 0;
                FFT_I(pitchShifter.fftData_, k) = 0;
            }

            /* do inverse transform */
            fft_.computeInverseA(&pitchShifter.fftData_[0]);

            /* do windowing and add to output accumulator */
            double crossfadeSize = grainSize - stepSize;
            for (k = 0; k < grainSize; k++) {
                window = -.5 * cos(M_2_PI * (double)k / (double)grainSize) + .5;
                pitchShifter.outputAccum_[k] += 2. * window * FFT_R(pitchShifter.fftData_, k) / (fftFrameSize2 * overlap);
            }
            for (k = 0; k < stepSize; k++) {
                pitchShifter.dstBuffer_[k] = pitchShifter.outputAccum_[k];
            }

            /* shift accumulator */
            std::memmove(&pitchShifter.outputAccum_[0], &pitchShifter.outputAccum_[0] + stepSize, grainSize * sizeof(double));

            /* move input FIFO */
            for (k = 0; k < latency; k++) {
                pitchShifter.srcBuffer_[k] = pitchShifter.srcBuffer_[k + stepSize];
            }
        }
    }

    return paContinue;
}

int
PortAudioLayer::PortAudioLayerImpl::paInputCallback(PortAudioLayer& parentLayer,
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
    double freq = 440.0f;
    auto elapsedTimePerFrame = 1.0 / static_cast<double>(parentLayer.audioFormat_.sample_rate);
    auto periodFrames = static_cast<int>((1.0 / freq) / elapsedTimePerFrame);

    for (; inFrame_ < endInFrame; inFrame_++) {
        // generate sine
        /*auto t = (inFrame_)* elapsedTimePerFrame;
        auto value = sin(M_2_PI * t * freq);
        auto int16Data = static_cast<AudioSample>(value * gain * FLOATTOINT16);
        inputBuffer_.get()->try_push(int16Data);*/

        // mic input
        inputBuffer_.get()->try_push(*in++);
    }

    return paContinue;
}