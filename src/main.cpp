#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <memory>

#include "app.h"
#include "pa_layer.h"
#include "soundfile.h"
#include "filter.h"
#include "Log.hpp"
#include "common.h"

#include "fft.h"
#include "pitchshift.h"

void pitchShift(    float pitchShiftRatio,
                    uint32_t grainSize_,
                    fft_wrapper& fft_,
                    long overlappingSamples,
                    long nSamples,
                    float sampleRate,
                    float* source,
                    float* out)
{
    constexpr const uint32_t maxFrameLength = 8192;
    std::array<float, maxFrameLength> srcRingBuffer_{};
    std::array<float, maxFrameLength> dstRingBuffer_{};
    std::array<kiss_fft_cpx, maxFrameLength> fftData_{};
    std::array<float, maxFrameLength / 2 + 1> lastPhase_{};
    std::array<float, maxFrameLength / 2 + 1> sumPhase_{};
    std::array<float, 2 * maxFrameLength> outputAccum_{};
    std::array<float, maxFrameLength> analysisFrq_{};
    std::array<float, maxFrameLength> analysisMag_{};
    std::array<float, maxFrameLength> synthesisFrq_{};
    std::array<float, maxFrameLength> synthesisMag_{};
    long inBufPos = 0;
    double magnitude, phase, tmp, window, real=0, imag=0;
    double freqPerBin, expectedFreq;
    long qpd, index, latency, stepSize, fftFrameSize2;
    int i, k;

    /* set up some handy variables */
    fftFrameSize2 = grainSize_ / 2;
    stepSize = grainSize_ / overlappingSamples;
    freqPerBin = sampleRate / (double)grainSize_;
    expectedFreq = 2.*M_PI*(double)stepSize / (double)grainSize_;
    latency = grainSize_ - stepSize;
    if (inBufPos == 0) {
        inBufPos = latency;
    }

    /* main processing loop */
    for (i = 0; i < nSamples; i++) {
        /* As long as we have not yet collected enough data just read in */
        srcRingBuffer_[inBufPos] = source[i];
        out[i] = dstRingBuffer_[inBufPos - latency];
        inBufPos++;

        /* now we have enough data for processing */
        if (inBufPos >= grainSize_) {
            inBufPos = latency;

            for (k = 0; k < grainSize_; k++) {
                window = -.5*cos(M_2_PI * (double)k / (double)grainSize_) + .5;
                FFT_R(fftData_, k) = srcRingBuffer_[k] * window;
                FFT_I(fftData_, k) = 0.;
            }

            /* ***************** ANALYSIS ******************* */
            /* do transform */
            fft_.computeA(&fftData_[0]);
            /* this is the analysis step */
            for (k = 0; k <= fftFrameSize2; k++) {

                /* de-interlace FFT buffer */
                real = FFT_R(fftData_, k);
                imag = FFT_I(fftData_, k);

                /* compute magnitude and phase */
                magnitude = 2. * sqrt(real*real + imag*imag);
                phase = atan2(imag, real);

                /* compute phase difference */
                tmp = phase - lastPhase_[k];
                lastPhase_[k] = phase;

                /* subtract expected phase difference */
                tmp -= (double)k*expectedFreq;

                /* map delta phase into +/- Pi interval */
                qpd = tmp / M_PI;
                if (qpd >= 0) qpd += qpd & 1;
                else qpd -= qpd & 1;
                tmp -= M_PI*(double)qpd;

                /* get deviation from bin frequency from the +/- Pi interval */
                tmp = overlappingSamples*tmp / (M_2_PI);

                /* compute the k-th partials' true frequency */
                tmp = (double)k*freqPerBin + tmp*freqPerBin;

                /* store magnitude and true frequency in analysis arrays */
                analysisMag_[k] = magnitude;
                analysisFrq_[k] = tmp;

            }

            /* ***************** PROCESSING ******************* */
            /* this does the actual pitch shifting */
            memset(&synthesisMag_, 0, grainSize_ * sizeof(float));
            memset(&synthesisFrq_, 0, grainSize_ * sizeof(float));
            for (k = 0; k <= fftFrameSize2; k++) {
                index = k*pitchShiftRatio;
                if (index <= fftFrameSize2) {
                    synthesisMag_[index] += analysisMag_[k];
                    synthesisFrq_[index] = analysisFrq_[k] * pitchShiftRatio;
                }
            }

            /* ***************** SYNTHESIS ******************* */
            /* this is the synthesis step */
            for (k = 0; k <= fftFrameSize2; k++) {

                /* get magnitude and true frequency from synthesis arrays */
                magnitude = synthesisMag_[k];
                tmp = synthesisFrq_[k];

                /* subtract bin mid frequency */
                tmp -= (double)k*freqPerBin;

                /* get bin deviation from freq deviation */
                tmp /= freqPerBin;

                /* take overlappingSamples into account */
                tmp = M_2_PI * tmp / overlappingSamples;

                /* add the overlap phase advance back in */
                tmp += (double)k*expectedFreq;

                /* accumulate delta phase to get bin phase */
                sumPhase_[k] += tmp;
                phase = sumPhase_[k];

                /* get real and imag part and re-interleave */
                FFT_R(fftData_, k) = magnitude*cos(phase);
                FFT_I(fftData_, k) = magnitude*sin(phase);
            }

            /* zero negative frequencies */
            for (k = grainSize_ / 2; k < grainSize_; k++) {
                FFT_R(fftData_, k) = 0;
                FFT_I(fftData_, k) = 0;
            }

            /* do inverse transform */
            fft_.computeInverseA(&fftData_[0]);

            /* do windowing and add to output accumulator */
            double crossfadeSize = grainSize_ - stepSize;
            for (k = 0; k < grainSize_; k++) {
                window = -.5 * cos(M_2_PI * (double)k / (double)grainSize_) + .5;
                outputAccum_[k] += 2. * window * FFT_R(fftData_, k) / (fftFrameSize2 * overlappingSamples);
            }
            for (k = 0; k < stepSize; k++) {
                dstRingBuffer_[k] = outputAccum_[k];
            }

            /* shift accumulator */
            std::memmove(&outputAccum_[0], &outputAccum_[0] + stepSize, grainSize_ * sizeof(double));

            /* move input FIFO */
            for (k = 0; k < latency; k++) {
                srcRingBuffer_[k] = srcRingBuffer_[k + stepSize];
            }
        }
    }
}

int main(int argc, char* argv[])
{
    uint16_t numChannels;
    uint32_t sampleRate;
    uint16_t bytesPerSample;
    std::vector<float> source, out;

    DBGOUT("loading file...");
    auto nSamples = ReadWaveFile("in.wav", source, numChannels, sampleRate, bytesPerSample);
    if (!nSamples) {
        return 0;
    }

    uint32_t frameSize = 1024;
    out.resize(nSamples);
    fft_wrapper fft_(frameSize);

    PitchShifter pitchShifter(sampleRate, frameSize, 4);

    std::chrono::time_point<std::chrono::steady_clock> start;
    auto pitchShiftRatio = 1.0;

    start = std::chrono::high_resolution_clock::now();
    pitchShifter.pitchShift(pitchShiftRatio, nSamples, &source[0], &out[0]);
    //pitchShift(pitchShiftRatio, 4, fft_, frameSize, nSamples, sampleRate, &source[0], &out[0]);
    DBGOUT("saving file...");
    WriteWaveFile("out2.wav", out, numChannels, sampleRate, bytesPerSample);
    std::cout << "Elapsed time: " << std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count() << " s\n";

    ("done.");

    //App app(512, 512);
    //app.execute(argc, argv);

#ifdef _WIN32
    system("pause");
#endif
    return 0;
}
