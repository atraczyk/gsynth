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

int main(int argc, char* argv[])
{
    uint16_t numChannels;
    uint32_t sampleRate;
    uint16_t bytesPerSample;
    std::vector<float> source;

    DBGOUT("loading file...");
    auto nSamples = ReadWaveFile("in2.wav", source, numChannels, sampleRate, bytesPerSample);
    if (!nSamples) {
        return 0;
    }

    LowPassFilter lpFilter(sampleRate, 4000.0, 1.0, 8.0);

    double bufSizeMultiplier = 2.0;
    uint32_t bufSize = 1536 * bufSizeMultiplier;
    uint32_t hopSize = 384 * bufSizeMultiplier;
    uint32_t grainSize = 512 * bufSizeMultiplier;
    uint32_t halfGrainSize = grainSize / 2;
    uint16_t grainsPerBuf = static_cast<double>(bufSize) / hopSize;
    auto nBuffers = static_cast<int>(ceil(static_cast<double>(nSamples) / bufSize));
    auto nGrains = grainsPerBuf * nBuffers;
    DBGOUT("samples: %d, buffers: %d, grains: %d", nSamples, nBuffers, nGrains);

    DBGOUT("preparing buffers...");
    std::vector<float> out(nSamples, 0);
    std::vector<std::vector<double>> buffers(nBuffers, std::vector<double>());
    std::vector<std::vector<double>> grains(nGrains);
    std::vector<std::vector<double>> phases(nGrains, std::vector<double>(grainSize));
    std::vector<std::vector<double>> magnitudes(nGrains, std::vector<double>(grainSize));
    std::vector<std::vector<kiss_fft_cpx>> rawffts(nGrains);

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
                    grainSize - truncatedBufSize,
                    std::back_inserter(grains.at(grainIndex))
                );
            }
        }
    }

    DBGOUT("analyzing...");
    fft_wrapper fft(sampleRate, grainSize);
    for (int i = 0; i < nGrains; ++i) {
        // lp filter
        for (int j = 0; j < grains.at(i).size(); ++j) {
            grains.at(i).at(j) = lpFilter.processSample(grains.at(i).at(j));
        }

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

        fft.setInput(&grains.at(i)[0]);
        
        // compute fft
        auto data = fft.computeStft();
        rawffts.at(i) = fft.getRawOutput();

        for (int j = 0; j < grainSize; ++j) {
            phases.at(i).at(j) = data.at(j).phase; 
            magnitudes.at(i).at(j) = data.at(j).amplitude;
        }
    }
    
    // frequency domain processing
    DBGOUT("processing...");
    double pitchShiftRatio = 1.0 / 2.0;
    std::vector<double> real(grainSize, 0);
    std::vector<double> imag(grainSize, 0);
    std::vector<double> omega(grainSize, 0);
    std::vector<std::vector<double>> newPhases(
        nGrains,
        std::vector<double>(grainSize, 0)
    );
    std::vector<std::vector<double>> deltaPhi(
        nGrains - 1,
        std::vector<double>(grainSize, 0)
    );

    // unwrap phases
    for (int j = 0; j < grainSize; ++j) {
        omega.at(j) = (M_2_PI * hopSize * j) / grainSize;
    }
    for(int i = 0; i < (nGrains - 1); ++i){
        for (int j = 0; j < grainSize; ++j) {
            if (i == 0){
                deltaPhi.at(i).at(j) = omega[j] + princArg(phases[i][j] - omega[j]);
            } else {
                deltaPhi.at(i).at(j) = omega[j] + princArg(phases[i][j] - phases[i-1][j] - omega[j]);
            }
        }
    }

    // compute new phases
    for(int j = 0; j < grainSize; ++j){
        newPhases[0][j] = phases[0][j];
    }
    for(int i = 1; i < nGrains; ++i){
        for(int j = 0; j < grainSize; ++j){
            newPhases[i][j] = princArg(newPhases[i-1][j] + deltaPhi[i-1][j] * pitchShiftRatio);
        }
    }

    DBGOUT("re-synthesizing...");
    for (int i = 0; i < nGrains; ++i) {
        // generate real/imag
        for(int j = 0; j < grainSize; ++j){
            real[j] = magnitudes[i][j] * cos(newPhases[i][j]);
            imag[j] = magnitudes[i][j] * sin(newPhases[i][j]);
        }

        fft.setOutput(&real[0], &imag[0]);

        // compute ifft
        grains.at(i) = fft.computeInverseStft();
        
        // fft shift (reverse)
        fftShift(grains.at(i));

        interpolate(grainSize, grains.at(i), pitchShiftRatio);
    }

    // apply crossfade window to grains
    DBGOUT("applying crossfade window to grains...");
    for (int i = 0; i < nGrains; ++i) {
        auto thisGrainSize = grains.at(i).size();
        for (int j = 0; j < thisGrainSize; ++j) {
            double mult = 1.0;
            auto crossfadeSize = thisGrainSize * 0.25;
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
    double gain = 2.0;
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

    // save file
    DBGOUT("saving file...");
    WriteWaveFile("out.wav", out, numChannels, sampleRate, bytesPerSample);

    DBGOUT("done.");

    //App app(512, 512);
    //app.execute(argc, argv);

#ifdef _WIN32
    system("pause");
#endif
    return 0;
}
