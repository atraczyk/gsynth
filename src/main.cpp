#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>

#include "app.h"
#include "pa_layer.h"
#include "soundfile.h"
#include "filter.h"
#include "Log.hpp"
#include "common.h"

#include "fft.h"

/****************************************************************************
*
* NAME: smbPitchShift.cpp
* VERSION: 1.2
* HOME URL: http://blogs.zynaptiq.com/bernsee
* KNOWN BUGS: none
*
* SYNOPSIS: Routine for doing pitch shifting while maintaining
* duration using the Short Time Fourier Transform.
*
* DESCRIPTION: The routine takes a pitchShift factor value which is between 0.5
* (one octave down) and 2. (one octave up). A value of exactly 1 does not change
* the pitch. numSampsToProcess tells the routine how many samples in indata[0...
* numSampsToProcess-1] should be pitch shifted and moved to outdata[0 ...
* numSampsToProcess-1]. The two buffers can be identical (ie. it can process the
* data in-place). fftFrameSize defines the FFT frame size used for the
* processing. Typical values are 1024, 2048 and 4096. It may be any value <=
* MAX_FRAME_LENGTH but it MUST be a power of 2. osamp is the STFT
* oversampling factor which also determines the overlap between adjacent STFT
* frames. It should at least be 4 for moderate scaling ratios. A value of 32 is
* recommended for best quality. sampleRate takes the sample rate for the signal
* in unit Hz, ie. 44100 for 44.1 kHz audio. The data passed to the routine in
* indata[] should be in the range [-1.0, 1.0), which is also the output range
* for the data, make sure you scale the data accordingly (for 16bit signed integers
* you would have to divide (and multiply) by 32768).
*
* COPYRIGHT 1999-2015 Stephan M. Bernsee <s.bernsee [AT] zynaptiq [DOT] com>
*
* 						The Wide Open License (WOL)
*
* Permission to use, copy, modify, distribute and sell this software and its
* documentation for any purpose is hereby granted without fee, provided that
* the above copyright notice and this license appear in all source copies.
* THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT EXPRESS OR IMPLIED WARRANTY OF
* ANY KIND. See http://www.dspguru.com/wol.htm for more information.
*
*****************************************************************************/

#include <string.h>
#include <math.h>
#include <stdio.h>

#define M_PI 3.14159265358979323846
#define MAX_FRAME_LENGTH 8192

void smbFft(float *fftBuffer, long fftFrameSize, long sign);
double smbAtan2(double x, double y);


// -----------------------------------------------------------------------------------------------------------------


void smbPitchShift(float pitchShift, long numSampsToProcess, long fftFrameSize, long osamp, float sampleRate, float *indata, float *outdata)
/*
Routine smbPitchShift(). See top of file for explanation
Purpose: doing pitch shifting while maintaining duration using the Short
Time Fourier Transform.
Author: (c)1999-2015 Stephan M. Bernsee <s.bernsee [AT] zynaptiq [DOT] com>
*/
{

    static float gInFIFO[MAX_FRAME_LENGTH];
    static float gOutFIFO[MAX_FRAME_LENGTH];
    static float gFFTworksp[2 * MAX_FRAME_LENGTH];
    static float gLastPhase[MAX_FRAME_LENGTH / 2 + 1];
    static float gSumPhase[MAX_FRAME_LENGTH / 2 + 1];
    static float gOutputAccum[2 * MAX_FRAME_LENGTH];
    static float gAnaFreq[MAX_FRAME_LENGTH];
    static float gAnaMagn[MAX_FRAME_LENGTH];
    static float gSynFreq[MAX_FRAME_LENGTH];
    static float gSynMagn[MAX_FRAME_LENGTH];
    static long gRover = false, gInit = false;
    double magn, phase, tmp, window, real, imag;
    double freqPerBin, expct;
    long i, k, qpd, index, inFifoLatency, stepSize, fftFrameSize2;

    /* set up some handy variables */
    fftFrameSize2 = fftFrameSize / 2;
    stepSize = fftFrameSize / osamp;
    freqPerBin = sampleRate / (double)fftFrameSize;
    expct = 2.*M_PI*(double)stepSize / (double)fftFrameSize;
    inFifoLatency = fftFrameSize - stepSize;
    if (gRover == false) gRover = inFifoLatency;

    /* initialize our static arrays */
    if (gInit == false) {
        memset(gInFIFO, 0, MAX_FRAME_LENGTH * sizeof(float));
        memset(gOutFIFO, 0, MAX_FRAME_LENGTH * sizeof(float));
        memset(gFFTworksp, 0, 2 * MAX_FRAME_LENGTH * sizeof(float));
        memset(gLastPhase, 0, (MAX_FRAME_LENGTH / 2 + 1) * sizeof(float));
        memset(gSumPhase, 0, (MAX_FRAME_LENGTH / 2 + 1) * sizeof(float));
        memset(gOutputAccum, 0, 2 * MAX_FRAME_LENGTH * sizeof(float));
        memset(gAnaFreq, 0, MAX_FRAME_LENGTH * sizeof(float));
        memset(gAnaMagn, 0, MAX_FRAME_LENGTH * sizeof(float));
        gInit = true;
    }

    /* main processing loop */
    for (i = 0; i < numSampsToProcess; i++) {

        /* As long as we have not yet collected enough data just read in */
        gInFIFO[gRover] = indata[i];
        outdata[i] = gOutFIFO[gRover - inFifoLatency];
        gRover++;

        /* now we have enough data for processing */
        if (gRover >= fftFrameSize) {
            DBGOUT("processing... i: %d", i);
            gRover = inFifoLatency;

            /* do windowing and re,im interleave */
            for (k = 0; k < fftFrameSize; k++) {
                window = -.5*cos(2.*M_PI*(double)k / (double)fftFrameSize) + .5;
                gFFTworksp[2 * k] = gInFIFO[k] * window;
                gFFTworksp[2 * k + 1] = 0.;
            }


            /* ***************** ANALYSIS ******************* */
            /* do transform */
            smbFft(gFFTworksp, fftFrameSize, -1);

            /* this is the analysis step */
            for (k = 0; k <= fftFrameSize2; k++) {

                /* de-interlace FFT buffer */
                real = gFFTworksp[2 * k];
                imag = gFFTworksp[2 * k + 1];

                /* compute magnitude and phase */
                magn = 2.*sqrt(real*real + imag*imag);
                phase = atan2(imag, real);

                /* compute phase difference */
                tmp = phase - gLastPhase[k];
                gLastPhase[k] = phase;

                /* subtract expected phase difference */
                tmp -= (double)k*expct;

                /* map delta phase into +/- Pi interval */
                qpd = tmp / M_PI;
                if (qpd >= 0) qpd += qpd & 1;
                else qpd -= qpd & 1;
                tmp -= M_PI*(double)qpd;

                /* get deviation from bin frequency from the +/- Pi interval */
                tmp = osamp*tmp / (2.*M_PI);

                /* compute the k-th partials' true frequency */
                tmp = (double)k*freqPerBin + tmp*freqPerBin;

                /* store magnitude and true frequency in analysis arrays */
                gAnaMagn[k] = magn;
                gAnaFreq[k] = tmp;

            }

            /* ***************** PROCESSING ******************* */
            /* this does the actual pitch shifting */
            memset(gSynMagn, 0, fftFrameSize * sizeof(float));
            memset(gSynFreq, 0, fftFrameSize * sizeof(float));
            for (k = 0; k <= fftFrameSize2; k++) {
                index = k*pitchShift;
                if (index <= fftFrameSize2) {
                    gSynMagn[index] += gAnaMagn[k];
                    gSynFreq[index] = gAnaFreq[k] * pitchShift;
                }
            }

            /* ***************** SYNTHESIS ******************* */
            /* this is the synthesis step */
            for (k = 0; k <= fftFrameSize2; k++) {

                /* get magnitude and true frequency from synthesis arrays */
                magn = gSynMagn[k];
                tmp = gSynFreq[k];

                /* subtract bin mid frequency */
                tmp -= (double)k*freqPerBin;

                /* get bin deviation from freq deviation */
                tmp /= freqPerBin;

                /* take osamp into account */
                tmp = 2.*M_PI*tmp / osamp;

                /* add the overlap phase advance back in */
                tmp += (double)k*expct;

                /* accumulate delta phase to get bin phase */
                gSumPhase[k] += tmp;
                phase = gSumPhase[k];

                /* get real and imag part and re-interleave */
                gFFTworksp[2 * k] = magn*cos(phase);
                gFFTworksp[2 * k + 1] = magn*sin(phase);
            }

            /* zero negative frequencies */
            for (k = fftFrameSize + 2; k < 2 * fftFrameSize; k++) gFFTworksp[k] = 0.;

            /* do inverse transform */
            smbFft(gFFTworksp, fftFrameSize, 1);

            /* do windowing and add to output accumulator */
            for (k = 0; k < fftFrameSize; k++) {
                window = -.5*cos(2.*M_PI*(double)k / (double)fftFrameSize) + .5;
                gOutputAccum[k] += 2.*window*gFFTworksp[2 * k] / (fftFrameSize2*osamp);
            }
            for (k = 0; k < stepSize; k++) gOutFIFO[k] = gOutputAccum[k];

            /* shift accumulator */
            memmove(gOutputAccum, gOutputAccum + stepSize, fftFrameSize * sizeof(float));

            /* move input FIFO */
            for (k = 0; k < inFifoLatency; k++) gInFIFO[k] = gInFIFO[k + stepSize];
        }
    }
}

// -----------------------------------------------------------------------------------------------------------------


void smbFft(float *fftBuffer, long fftFrameSize, long sign)
/*
FFT routine, (C)1996 S.M.Bernsee. Sign = -1 is FFT, 1 is iFFT (inverse)
Fills fftBuffer[0...2*fftFrameSize-1] with the Fourier transform of the
time domain data in fftBuffer[0...2*fftFrameSize-1]. The FFT array takes
and returns the cosine and sine parts in an interleaved manner, ie.
fftBuffer[0] = cosPart[0], fftBuffer[1] = sinPart[0], asf. fftFrameSize
must be a power of 2. It expects a complex input signal (see footnote 2),
ie. when working with 'common' audio signals our input signal has to be
passed as {in[0],0.,in[1],0.,in[2],0.,...} asf. In that case, the transform
of the frequencies of interest is in fftBuffer[0...fftFrameSize].
*/
{
    float wr, wi, arg, *p1, *p2, temp;
    float tr, ti, ur, ui, *p1r, *p1i, *p2r, *p2i;
    long i, bitm, j, le, le2, k;

    for (i = 2; i < 2 * fftFrameSize - 2; i += 2) {
        for (bitm = 2, j = 0; bitm < 2 * fftFrameSize; bitm <<= 1) {
            if (i & bitm) j++;
            j <<= 1;
        }
        if (i < j) {
            p1 = fftBuffer + i; p2 = fftBuffer + j;
            temp = *p1; *(p1++) = *p2;
            *(p2++) = temp; temp = *p1;
            *p1 = *p2; *p2 = temp;
        }
    }
    for (k = 0, le = 2; k < (long)(log(fftFrameSize) / log(2.) + .5); k++) {
        le <<= 1;
        le2 = le >> 1;
        ur = 1.0;
        ui = 0.0;
        arg = M_PI / (le2 >> 1);
        wr = cos(arg);
        wi = sign*sin(arg);
        for (j = 0; j < le2; j += 2) {
            p1r = fftBuffer + j; p1i = p1r + 1;
            p2r = p1r + le2; p2i = p2r + 1;
            for (i = j; i < 2 * fftFrameSize; i += le) {
                tr = *p2r * ur - *p2i * ui;
                ti = *p2r * ui + *p2i * ur;
                *p2r = *p1r - tr; *p2i = *p1i - ti;
                *p1r += tr; *p1i += ti;
                p1r += le; p1i += le;
                p2r += le; p2i += le;
            }
            tr = ur*wr - ui*wi;
            ui = ur*wi + ui*wr;
            ur = tr;
        }
    }
}


// -----------------------------------------------------------------------------------------------------------------

/*

12/12/02, smb

PLEASE NOTE:

There have been some reports on domain errors when the atan2() function was used
as in the above code. Usually, a domain error should not interrupt the program flow
(maybe except in Debug mode) but rather be handled "silently" and a global variable
should be set according to this error. However, on some occasions people ran into
this kind of scenario, so a replacement atan2() function is provided here.

If you are experiencing domain errors and your program stops, simply replace all
instances of atan2() with calls to the smbAtan2() function below.

*/


double smbAtan2(double x, double y)
{
    double signx;
    if (x > 0.) signx = 1.;
    else signx = -1.;

    if (x == 0.) return 0.;
    if (y == 0.) return signx * M_PI / 2.;

    return atan2(x, y);
}


// -----------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------

double
princArg(double phaseIn)
{
    double a = phaseIn / M_2_PI;
    return phaseIn - round(a) * M_2_PI;
}

void
interpolate(uint32_t bufSize,
            std::vector<double>& buffer,
            double pitchShiftRatio)
{
    std::vector<double> bufferCpy { buffer };
    auto factor = 1.0 / pitchShiftRatio;
    auto x1 = 0;
    
    for (int i = 0; i < bufSize - 1; ++i) {
        auto y1 = bufferCpy.at(i);
        auto x2 = x1 + factor;
        auto y2 = bufferCpy.at(i + 1);

        for (int j = 0; j < (floor(factor) + 1); ++j) {
            auto xt = std::min((int)bufSize - 1, x1 + j);
            auto yt = (y2 - y1) / (x2 - x1) * (xt - x1) + y1;
            buffer.at(xt) = yt;
        }
        x1 = x2;
    }
}

template<typename T>
void
fftShift(T buffer)
{
    std::rotate(
        buffer.begin(),
        buffer.begin() + buffer.size() / 2,
        buffer.end()
    );
}

void pitchShift1(float pitchShiftRatio, long nSamples, float sampleRate, const std::vector<float>& source, std::vector<float>& out)
{
    LowPassFilter lpFilter(sampleRate, 4000.0, 1.0, 8.0);

    double bufSizeMultiplier = 2.0;
    uint32_t bufSize = 1536 * bufSizeMultiplier;
    uint32_t hopSize = 256 * bufSizeMultiplier;
    uint32_t grainSize = 512 * bufSizeMultiplier;
    uint32_t halfGrainSize = grainSize / 2;
    uint16_t grainsPerBuf = static_cast<double>(bufSize) / hopSize;
    auto nBuffers = static_cast<int>(ceil(static_cast<double>(nSamples) / bufSize));
    auto nGrains = grainsPerBuf * nBuffers;
    DBGOUT("samples: %d, buffers: %d, grains: %d", nSamples, nBuffers, nGrains);

    DBGOUT("preparing buffers...");
    std::vector<std::vector<double>> buffers(nBuffers, std::vector<double>());
    std::vector<std::vector<double>> grains(nGrains);
    std::vector<std::vector<double>> phases(nGrains, std::vector<double>(grainSize));
    std::vector<std::vector<double>> magnitudes(nGrains, std::vector<double>(grainSize));
    std::vector<std::vector<kiss_fft_cpx>> rawffts(nGrains);

    ///////////////////////////////////////////////////
    std::vector<std::vector<double>> gLastPhase(
        nGrains,
        std::vector<double>(grainSize, 0)
    );
    std::vector<std::vector<double>> gSumPhase(
        nGrains,
        std::vector<double>(grainSize, 0)
    );
    std::vector<std::vector<double>> gAnaFreq(
        nGrains,
        std::vector<double>(grainSize, 0)
    );
    std::vector<std::vector<double>> gAnaMagn(
        nGrains,
        std::vector<double>(grainSize, 0)
    );
    std::vector<std::vector<double>> gSynFreq(
        nGrains,
        std::vector<double>(grainSize, 0)
    );
    std::vector<std::vector<double>> gSynMagn(
        nGrains,
        std::vector<double>(grainSize, 0)
    );
    
    double magn, phase, tmp, window, real, imag;
    double freqPerBin, expct;
    long i, k, qpd, index, inFifoLatency, stepSize, fftFrameSize2;
    long osamp = 8;

    /* set up some handy variables */
    fftFrameSize2 = grainSize / 2;
    stepSize = grainSize / osamp;
    freqPerBin = sampleRate / (double)grainSize;
    expct = M_2_PI * (double)stepSize / (double)grainSize;

    DBGOUT("serializing to buffers...");
    for (int i = 0; i < nBuffers; ++i) {
        size_t fromPos = i * bufSize;
        size_t toPos = (i < nBuffers - 1) ? i * bufSize + bufSize : i * bufSize + (nSamples - (i * bufSize));
        auto from = source.begin() + fromPos;
        auto to = source.begin() + toPos;
        std::copy_n(
            source.begin() + fromPos,
            toPos - fromPos,
            std::back_inserter(buffers.at(i))
        );
    }

    DBGOUT("splitting into grains...");
    for (int i = 0; i < nBuffers; ++i) {
        size_t bufPos = 0;
        for (int j = 0; j < grainsPerBuf; ++j) {
            auto grainIndex = i * grainsPerBuf + j;
            if (j < grainsPerBuf - 1 && i < nBuffers - 1) {
                std::copy_n(
                    buffers.at(i).begin() + bufPos,
                    grainSize,
                    std::back_inserter(grains.at(grainIndex))
                );
                bufPos += hopSize;
                continue;
            }
            uint32_t truncatedBufSize = (buffers.at(i).size() - bufPos);
            std::copy_n(
                buffers.at(i).begin() + bufPos,
                std::min(truncatedBufSize, grainSize),
                std::back_inserter(grains.at(grainIndex))
            );
            if (i < nBuffers - 1 && truncatedBufSize < grainSize) {
                std::copy_n(
                    buffers.at(i + 1).begin(),
                    std::min(
                        (size_t)(grainSize - truncatedBufSize),
                        buffers.at(i + 1).size()
                    ),
                    std::back_inserter(grains.at(grainIndex))
                );
            }
        }
    }

    DBGOUT("analyzing...");
    fft_wrapper fft(sampleRate, grainSize);
    for (int i = 0; i < nGrains; ++i) {
        // windowing
        auto dgrainSize = static_cast<double>(grainSize);
        double window;
        int windowType = 0;
        for (double j = 0; j < grains.at(i).size(); ++j) {
            if (windowType == 0) { // Hann
                window = 0.5 * (1.0 - cos(2.0 * M_PI * j / (dgrainSize - 1)));
            }
            else if(windowType == 1) { // Hamming
                window = 0.54 - 0.46 * cos(2.0 * M_PI * j / (dgrainSize - 1));
            }
            else if (windowType == 2) { // Welch
                auto f = (((j - (dgrainSize - 1) / 2)) / ((dgrainSize - 1) / 2));
                window = 1.0 - f * f;
            }
            else if (windowType == 3) { // Blackman-Harris
                auto f = M_PI * j / (dgrainSize - 1);
                window =    0.35875 -
                            0.48829 * cos(2 * f) +
                            0.14128 * cos(4 * f) -
                            0.01168 * cos(6 * f);
            }
            grains.at(i).at(j) *= window;
        }

        // fft shift
        fftShift(grains.at(i));

        fft.setData(FFTDirection::In, grains.at(i));
        
        // compute fft
        auto data = fft.computeStft();
        rawffts.at(i) = fft.getRawOutput();

        for (int j = 0; j < grainSize; ++j) {
            phases.at(i).at(j) = data.at(j).phase; 
            magnitudes.at(i).at(j) = data.at(j).amplitude;
        }

        for (int k = 0; k <= grains.at(i).size() / 2; k++) {
            /* de-interlace FFT buffer */
            real = rawffts[i][k].r;
            imag = rawffts[i][k].i;

            /* compute magnitude and phase */
            magn = 2. * sqrt(real * real + imag * imag);
            phase = atan2(imag, real);

            /* compute phase difference */
            tmp = phase - gLastPhase[i][k];
            gLastPhase[i][k] = phase;

            /* subtract expected phase difference */
            tmp -= (double)k * expct;

            /* map delta phase into +/- Pi interval */
            qpd = tmp / M_PI;
            if (qpd >= 0) {
                qpd += qpd & 1;
            }
            else {
                qpd -= qpd & 1;
            }
            tmp -= M_PI * (double)qpd;

            /* get deviation from bin frequency from the +/- Pi interval */
            tmp = osamp * tmp / M_2_PI;

            /* compute the k-th partials' true frequency */
            tmp = (double)k * freqPerBin + tmp * freqPerBin;

            /* store magnitude and true frequency in analysis arrays */
            gAnaMagn[i][k] = magn;
            gAnaFreq[i][k] = tmp;

        }
    }
    
    // frequency domain processing
    DBGOUT("processing...");
    //if (0) {
    //    //std::vector<double> real(grainSize, 0);
    //    //std::vector<double> imag(grainSize, 0);
    //    std::vector<double> omega(grainSize, 0);
    //    std::vector<std::vector<double>> newPhases(
    //        nGrains,
    //        std::vector<double>(grainSize, 0)
    //    );
    //    std::vector<std::vector<double>> deltaPhi(
    //        nGrains - 1,
    //        std::vector<double>(grainSize, 0)
    //    );

    //    // unwrap phases
    //    for (int j = 0; j < grainSize; ++j) {
    //        omega.at(j) = (M_2_PI * hopSize * j) / grainSize;
    //    }
    //    for (int i = 0; i < (nGrains - 1); ++i) {
    //        for (int j = 0; j < grainSize; ++j) {
    //            if (i == 0) {
    //                deltaPhi.at(i).at(j) = omega[j] + princArg(phases[i][j] - omega[j]);
    //            }
    //            else {
    //                deltaPhi.at(i).at(j) = omega[j] + princArg(phases[i][j] - phases[i - 1][j] - omega[j]);
    //            }
    //        }
    //    }

    //    // compute new phases
    //    for (int j = 0; j < grainSize; ++j) {
    //        newPhases[0][j] = phases[0][j];
    //    }
    //    for (int i = 1; i < nGrains; ++i) {
    //        for (int j = 0; j < grainSize; ++j) {
    //            newPhases[i][j] = princArg(newPhases[i - 1][j] + deltaPhi[i - 1][j] * pitchShiftRatio);
    //        }
    //    }
    //}

    /* ***************** PROCESSING ******************* */
    /* this does the actual pitch shifting */
    for (int i = 0; i < nGrains; ++i) {
        std::fill_n(gSynMagn[i].begin(), gSynMagn[i].size(), 0);
        std::fill_n(gSynFreq[i].begin(), gSynFreq[i].size(), 0);
        for (k = 0; k <= fftFrameSize2; k++) {
            index = k * pitchShiftRatio;
            if (index <= fftFrameSize2) {
                gSynMagn[i][index] += gAnaMagn[i][k];
                gSynFreq[i][index] = gAnaFreq[i][k] * pitchShiftRatio;
            }
        }
    }

    DBGOUT("re-synthesizing...");
    std::vector<double> reals(grainSize, 0);
    std::vector<double> imags(grainSize, 0);

    for (int i = 0; i < nGrains; ++i) {
        for (k = 0; k < grainSize; k++) {

            /* get magnitude and true frequency from synthesis arrays */
            magn = gSynMagn[i][k];
            tmp = gSynFreq[i][k];

            /* subtract bin mid frequency */
            tmp -= (double)k*freqPerBin;

            /* get bin deviation from freq deviation */
            tmp /= freqPerBin;

            /* take osamp into account */
            tmp = 2.*M_PI*tmp / osamp;

            /* add the overlap phase advance back in */
            tmp += (double)k*expct;

            /* accumulate delta phase to get bin phase */
            gSumPhase[i][k] += tmp;
            phase = gSumPhase[i][k];

            /* get real and imag part and re-interleave */
            reals[k] =  magn*cos(phase);
            imags[k] = magn*sin(phase);
        }

        fft.setData(FFTDirection::Out, reals, imags);

        // generate real/imag
        for(int j = 0; j < grainSize; ++j){
            //real[j] = magnitudes[i][j] * cos(newPhases[i][j]);
            //imag[j] = magnitudes[i][j] * sin(newPhases[i][j]);
        }

        //fft.setData(FFTDirection::Out , real, imag);

        // compute ifft
        grains.at(i) = fft.computeInverseStft();
        
        // fft shift (reverse)
        fftShift(grains.at(i));

        //interpolate(grainSize, grains.at(i), pitchShiftRatio);
    }

    // apply crossfade window to grains
    DBGOUT("applying crossfade window to grains...");
    for (int i = 0; i < nGrains; ++i) {
        auto thisGrainSize = grains.at(i).size();
        double crossfadeSize = thisGrainSize - hopSize;
        for (int j = 0; j < thisGrainSize; ++j) {
            double mult = 1.0;
            if (j < crossfadeSize) {
                mult = j / (crossfadeSize);
            }
            else if (crossfadeSize < j && j < thisGrainSize - crossfadeSize) {
                mult = 1.0;
            }
            else if (thisGrainSize - crossfadeSize < j) {
                mult = (-j / crossfadeSize) + (thisGrainSize / crossfadeSize);
            }
            grains.at(i).at(j) *= mult;
        }
    }

    // gain adjust
    double gain = 1.5;
    for (int i = 0; i < nGrains; ++i) {
        for (int j = 0; j < grainSize; ++j) {
            grains.at(i).at(j) *= gain;
        }
    }

    // overlap-add grains to output buffer
    DBGOUT("overlap-adding grains to output buffer...");
    size_t bufPos = 0;
    for (int i = 0; i < nGrains; ++i) {
        if (i < nGrains - 1 && bufPos + grains.at(i).size() < out.size()) {
            for (int j = 0; j < grains.at(i).size() - 1; ++j) {
                out.at(bufPos + j) += grains.at(i).at(j);
            }
            bufPos += hopSize;
            continue;
        }
        for (int k = 0; k < out.size() - bufPos - 1; ++k) {
            out.at(bufPos + k) += grains.at(i).at(k);
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
    auto nSamples = ReadWaveFile("in2.wav", source, numChannels, sampleRate, bytesPerSample);
    if (!nSamples) {
        return 0;
    }

    out.resize(nSamples);
    
    int method = 1;
    if (method == 0) {
        smbPitchShift(0.8, nSamples, 1024, 4, sampleRate, &source[0], &out[0]);
        DBGOUT("saving file...");
        WriteWaveFile("out0.wav", out, numChannels, sampleRate, bytesPerSample);
    }
    else if (method == 1) {
        pitchShift1(1.8, nSamples, sampleRate, source, out);
        DBGOUT("saving file...");
        WriteWaveFile("out1.wav", out, numChannels, sampleRate, bytesPerSample);
    }

    DBGOUT("done.");

    //App app(512, 512);
    //app.execute(argc, argv);

#ifdef _WIN32
    //system("pause");
#endif
    return 0;
}
