#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>

#include "app.h"
#include "pa_layer.h"
#include "soundfile.h"
#include "Log.hpp"

#include "fft.h"

double
princArg(double phaseIn)
{
    double a = phaseIn / M_2_PI;
    return phaseIn - round(a) * M_2_PI;
}

std::vector<double>
interpolate(uint32_t bufSize,
            std::vector<double> buffer,
            double pitchShiftRatio)
{
    std::vector<double> outBuffer(bufSize, 0);
    auto factor = 1.0 / pitchShiftRatio;
    auto x1 = 0;
    
    for (int i = 0; i < bufSize; ++i) {
        auto y1 = buffer.at(i);
        auto x2 = x1 + factor;
        auto y2 = buffer.at(i + 1);

        for (int j = 0; j < (floor(factor) + 1); ++i) {
            auto xt = x1 + j;
            auto yt = (y2 - y1) / (x2 - x1) * (xt - x1) + y1;
            outBuffer.at(xt) = yt;
        }
        x1 = x2;
    }
}

int main(int argc, char* argv[])
{
    uint16_t numChannels;
    uint32_t sampleRate;
    uint16_t bytesPerSample;
    std::vector<float> source;

    DBGOUT("loading file...");
    auto nSamples = ReadWaveFile("in.wav", source, numChannels, sampleRate, bytesPerSample);
    if (!nSamples) {
        return 0;
    }

    uint32_t bufSize = 1536;
    uint32_t hopSize = 384;
    uint32_t grainSize = 512;
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
    //out.resize(nSamples);
    //~ buffers.resize(nBuffers);
    //~ grains.resize(nGrains);
    //~ phases.resize(nGrains);
    //~ magnitudes.resize(nGrains);
    //~ for (auto& phase: phases) {
        //~ phase.resize(grainSize);
    //~ }
    //~ for (auto& magnitude: magnitudes) {
        //~ magnitude.resize(grainSize);
    //~ }
    //~ std::fill_n(out.begin(), out.size(), 0);
    //~ for (auto& buf : buffers) {
        //~ std::fill_n(buf.begin(), buf.size(), 0);
    //~ }

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
        // fft shift
        DBGOUT("fft shift...");
        std::rotate(grains.at(i).begin(),
                    grains.at(i).begin() + halfGrainSize,
                    grains.at(i).end());
        fft.setInput(&grains.at(i)[0]);
        
        // compute fft
        DBGOUT("computing fft...");
        auto data = fft.computeStft();
        //fft.shift(); // ?
        for (int j = 0; j < halfGrainSize; ++j) {
            phases.at(i).at(j) = data.at(j).phase; 
            magnitudes.at(i).at(j) = data.at(j).amplitude;
            phases.at(i).at(grainSize - j - 1) = data.at(j).phase;
            magnitudes.at(i).at(grainSize - j - 1) = data.at(j).amplitude;
        }
    }
    
    // frequency domain processing
    DBGOUT("processing...");
    double pitchShiftRatio = 1.0;
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
        for (int j = 0; j < (nGrains); ++j) {
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

        // compute ifft
        DBGOUT("compute ifft...");
        //fft.shift(); // ?
        grains.at(i) = fft.computeInverseStft();
        
        // fft shift (reverse)
        DBGOUT("fft shift (reverse)...");
        std::rotate(grains.at(i).begin(),
                    grains.at(i).begin() + halfGrainSize,
                    grains.at(i).end());
    }

    //~ // apply crossfade window to grains
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
