#include "pitchshift.h"
#include "common.h"
#include "Log.hpp"

PitchShifter::PitchShifter()
{
}

PitchShifter::PitchShifter(const uint32_t& sampleRate, const uint32_t& grainSize, const uint32_t& overlap)
    : grainSize_(grainSize)
    , sampleRate_(sampleRate)
    , overlap_(overlap)
    , fft_(grainSize)
{
}

PitchShifter::~PitchShifter()
{
}

void
PitchShifter::pitchShift(
    float pitchShiftRatio,
    long nSamples,
    float* source,
    float* out)
{
    this->pitchShift(pitchShiftRatio, grainSize_ , overlap_, nSamples, sampleRate_, source, out);
}

void
PitchShifter::pitchShift(
    float pitchShiftRatio,
    uint32_t grainSize,
    long overlap,
    long nSamples,
    float sampleRate,
    float* source,
    float* out)
{
    static long inBufPos = 0;
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
    static int i, k;

    grainSize_ = grainSize;
    fftFrameSize2 = grainSize / 2;
    stepSize = grainSize_ / overlap;
    freqPerBin = sampleRate / (double)grainSize;
    expectedFreq = M_2_PI * (double)stepSize / (double)grainSize;
    latency = grainSize - stepSize;
    if (inBufPos == 0) {
        inBufPos = latency;
    }

    /* main processing loop */
    for (i = 0; i < nSamples; i++) {
        /* As long as we have not yet collected enough data just read in */
        srcBuffer_[inBufPos] = source[i];
        out[i] = dstBuffer_[inBufPos - latency];
        inBufPos++;

        /* now we have enough data for processing */
        if (inBufPos >= grainSize) {
            inBufPos = latency;

            for (k = 0; k < grainSize; k++) {
                window = -.5*cos(M_2_PI * (double)k / (double)grainSize) + .5;
                FFT_R(fftData_, k) = srcBuffer_[k] * window;
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
                analysisMag_[k] = magnitude;
                analysisFrq_[k] = tmp;

            }

            /* ***************** PROCESSING ******************* */
            /* this does the actual pitch shifting */
            memset(&synthesisMag_, 0, grainSize * sizeof(float));
            memset(&synthesisFrq_, 0, grainSize * sizeof(float));
            for (k = 0; k <= fftFrameSize2; k++) {
                index = k*pitchShiftRatio;
                if (index <= fftFrameSize2) {
                    synthesisMag_[index] = analysisMag_[k];
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
                tmp = M_2_PI * tmp / overlap;

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
            for (k = grainSize / 2; k < grainSize; k++) {
                FFT_R(fftData_, k) = 0;
                FFT_I(fftData_, k) = 0;
            }

            /* do inverse transform */
            fft_.computeInverseA(&fftData_[0]);

            /* do windowing and add to output accumulator */
            double crossfadeSize = grainSize - stepSize;
            for (k = 0; k < grainSize; k++) {
                window = -.5 * cos(M_2_PI * (double)k / (double)grainSize) + .5;
                outputAccum_[k] += 2. * window * FFT_R(fftData_, k) / (fftFrameSize2 * overlap);
            }
            for (k = 0; k < stepSize; k++) {
                dstBuffer_[k] = outputAccum_[k];
            }

            /* shift accumulator */
            std::memmove(&outputAccum_[0], &outputAccum_[0] + stepSize, grainSize * sizeof(double));

            /* move input FIFO */
            for (k = 0; k < latency; k++) {
                srcBuffer_[k] = srcBuffer_[k + stepSize];
            }
        }
    }
}