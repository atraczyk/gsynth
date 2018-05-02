#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>

#include "app.h"
#include "pa_layer.h"
#include "soundfile.h"
#include "Log.hpp"

#include "fft.h"

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
