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
    std::vector<float> source, out;
    std::vector<std::vector<double>> buffers;
    std::vector<std::vector<double>> grains;
    std::vector<double> bufferOverlap;
    uint32_t bufSize = 1536;
    uint32_t hopSize = 384;
    uint32_t grainSize = 512;
    uint16_t grainsPerBuf = static_cast<double>(bufSize) / hopSize;
    std::vector<std::vector<double>> phases;
    std::vector<std::vector<double>> magnitudes;

    // load file
    DBGOUT("loading file...");
    auto nSamples = ReadWaveFile("in.wav", source, numChannels, sampleRate, bytesPerSample);
    if (!nSamples) {
        return 0;
    }
    auto nBuffers = static_cast<int>(ceil(static_cast<double>(nSamples) / bufSize));
    auto nGrains = grainsPerBuf * nBuffers;
    DBGOUT("samples: %d, buffers: %d, grains: %d", nSamples, nBuffers, nGrains);

    // prepare buffers
    DBGOUT("preparing buffers...");
    out.resize(nSamples);
    buffers.resize(nBuffers);
    grains.resize(nGrains);
    phases.resize(nGrains);
    magnitudes.resize(nGrains);
    for (auto& phase: phases) {
        phase.resize(grainSize);
    }
    for (auto& magnitude: magnitudes) {
        magnitude.resize(grainSize);
    }
    std::fill_n(out.begin(), out.size(), 0);
    for (auto& buf : buffers) {
        std::fill_n(buf.begin(), buf.size(), 0);
    }

    // serialize to buffers
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

    // split into grains
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

    // process
    DBGOUT("processing...");
    fft_wrapper fft(sampleRate, grainSize);
    std::vector<double> real(grainSize, 0);
    std::vector<double> imag(grainSize, 0);
    std::vector<double> newPhase(grainSize, 0);
    std::vector<double> omega(grainSize, 0);
    std::vector<std::vector<double>> deltaPhi(
        nGrains - 1,
        std::vector<double>(grainSize, 0)
    );
    for (int j = 0; j < grainSize; ++j) {
        omega.at(j) = (M_2_PI * hopSize * j) / grainSize;
    }
    
    for (int i = 0; i < nGrains; ++i) {
        auto halfGrainSize = grainSize / 2;
        
        // fft shift
        DBGOUT("fft shift...");
        //~ std::rotate(grains.at(i).begin(),
                    //~ grains.at(i).begin() + halfGrainSize,
                    //~ grains.at(i).end());
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

        // frequency domain processing
        // 1 unwrap phases
        

        for (int j = 0; j < (nGrains); ++j) {
            
        }
        // 2 compute new phases
        // 3 generate real/imag

        // compute ifft
        DBGOUT("compute ifft...");
        //fft.shift(); // ?
        grains.at(i) = fft.computeInverseStft();
        
        // fft shift (reverse)
        DBGOUT("fft shift (reverse)...");
        //~ std::rotate(grains.at(i).begin(),
                    //~ grains.at(i).begin() + halfGrainSize,
                    //~ grains.at(i).end());
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
