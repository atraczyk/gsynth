#include <stdio.h>
#include <memory.h>
#include <inttypes.h>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <math.h>

#ifdef __unix
#define fopen_s(pFile,filename,mode) ((*(pFile))=fopen((filename),  (mode)))==NULL
#endif

#include "soundfile.h"
#include "Log.hpp"

inline void FloatToPCM(unsigned char *PCM, const float& in, size_t numBytes)
{
    // 8 bit is unsigned
    if (numBytes == 1) {
        PCM[0] = static_cast<unsigned char>((in * 0.5f + 0.5f) * 255.0f);
        return;
    }

    // casting to double because floats can't exactly store 0x7fffffff, but doubles can.
    // Details of that: https://blog.demofox.org/2017/11/21/floating-point-precision/
    uint32_t data;
    if (in < 0.0f)
        data = uint32_t(double(in) * double(0x80000000));
    else
        data = uint32_t(double(in) * double(0x7fffffff));

    switch (numBytes) {
    case 4: PCM[3] = ((data >> 24) & 0xFF); PCM[2] = ((data >> 16) & 0xFF); PCM[1] = ((data >> 8) & 0xFF); PCM[0] = (data & 0xFF); break;
    case 3: PCM[2] = ((data >> 24) & 0xFF); PCM[1] = ((data >> 16) & 0xFF); PCM[0] = ((data >> 8) & 0xFF); break;
    case 2: PCM[1] = ((data >> 24) & 0xFF); PCM[0] = ((data >> 16) & 0xFF); break;
    }
}

inline void PCMToFloat(float& out, const unsigned char *PCM, size_t numBytes)
{
    // 8 bit is unsigned
    if (numBytes == 1) {
        out = (float(PCM[0]) / float(255.0f)) * 2.0f - 1.0f;
        return;
    }

    uint32_t data = 0;
    switch (numBytes) {
    case 4: data = (uint32_t(PCM[3]) << 24) | (uint32_t(PCM[2]) << 16) | (uint32_t(PCM[1]) << 8) | uint32_t(PCM[0]); break;
    case 3: data = (uint32_t(PCM[2]) << 24) | (uint32_t(PCM[1]) << 16) | (uint32_t(PCM[0]) << 8); break;
    case 2: data = (uint32_t(PCM[1]) << 24) | (uint32_t(PCM[0]) << 16); break;
    }

    // casting to double because floats can't exactly store 0x7fffffff, but doubles can.
    // Details of that: https://blog.demofox.org/2017/11/21/floating-point-precision/
    if (data & 0x80000000)
        out = float(double(int32_t(data)) / double(0x80000000));
    else
        out = float(double(data) / double(0x7fffffff));
}

// numBytes can be 1, 2, 3, or 4.
// Coresponding to 8 bit, 16 bit, 24 bit, and 32 bit audio.
bool WriteWaveFile(const char *fileName, std::vector<float>& dataFloat, uint16_t numChannels, uint32_t sampleRate, uint16_t numBytes)
{
    if (dataFloat.size() == 0) {
        return false;
    }
    std::vector<unsigned char> data;
    data.resize(dataFloat.size() * numBytes);
    for (size_t i = 0; i < dataFloat.size(); ++i)
        FloatToPCM((unsigned char*)&data[i*numBytes], dataFloat[i], numBytes);

    uint32_t dataSize = (uint32_t)data.size();
    uint16_t bitsPerSample = numBytes * 8;

    //open the file if we can
    FILE *File = nullptr;
    fopen_s(&File, fileName, "w+b");
    if (!File) {
        DBGOUT("[-----ERROR-----] Could not open %s for writing.", fileName);
        return false;
    }

    SMinimalWaveFileHeader waveHeader;

    //fill out the main chunk
    memcpy(waveHeader.m_chunkID, "RIFF", 4);
    waveHeader.m_chunkSize = dataSize + 36;
    memcpy(waveHeader.m_format, "WAVE", 4);

    //fill out sub chunk 1 "fmt "
    memcpy(waveHeader.m_subChunk1ID, "fmt ", 4);
    waveHeader.m_subChunk1Size = 16;
    waveHeader.m_audioFormat = 1;
    waveHeader.m_numChannels = numChannels;
    waveHeader.m_sampleRate = sampleRate;
    waveHeader.m_byteRate = sampleRate * numChannels * bitsPerSample / 8;
    waveHeader.m_blockAlign = numChannels * bitsPerSample / 8;
    waveHeader.m_bitsPerSample = bitsPerSample;

    //fill out sub chunk 2 "data"
    memcpy(waveHeader.m_subChunk2ID, "data", 4);
    waveHeader.m_subChunk2Size = dataSize;

    //write the header
    fwrite(&waveHeader, sizeof(SMinimalWaveFileHeader), 1, File);

    //write the wave data itself
    fwrite(&data[0], dataSize, 1, File);

    //close the file and return success
    fclose(File);
    DBGOUT("%s saved.", fileName);
    return true;
}

bool ReadFileIntoMemory(const char *fileName, std::vector<unsigned char>& data)
{
    //open the file if we can
    FILE *file = nullptr;
    fopen_s(&file, fileName, "rb");
    if (!file) {
        DBGOUT("[-----ERROR-----]Could not open %s for reading.", fileName);
        return false;
    }

    // get the file size and resize the vector to hold the data
    fseek(file, 0, SEEK_END);
    data.resize(ftell(file));

    // read the file into the vector
    fseek(file, 0, SEEK_SET);
    fread(&data[0], 1, data.size(), file);

    // return success
    fclose(file);
    return true;
}

uint64_t ReadWaveFile(const char *fileName, std::vector<float>& data, uint16_t& numChannels, uint32_t& sampleRate, uint16_t& numBytes)
{
    // read the whole file into memory if we can
    std::vector<unsigned char> fileData;
    if (!ReadFileIntoMemory(fileName, fileData))
        return 0;
    size_t fileIndex = 0;

    //make sure the main chunk ID is "RIFF"
    if ((fileData.size() < fileIndex + 4) || memcmp(&fileData[fileIndex], "RIFF", 4)) {
        DBGOUT("[-----ERROR-----]%s is an invalid input file. (1)", fileName);
        return 0;
    }
    fileIndex += 4;

    //get the main chunk size
    uint32_t chunkSize;
    if (fileData.size() < fileIndex + 4) {
        DBGOUT("[-----ERROR-----]%s is an invalid input file. (2)", fileName);
        return 0;
    }
    chunkSize = *(uint32_t*)&fileData[fileIndex];
    fileIndex += 4;

    //make sure the format is "WAVE"
    if ((fileData.size() < fileIndex + 4) || memcmp(&fileData[fileIndex], "WAVE", 4)) {
        DBGOUT("[-----ERROR-----]%s is an invalid input file. (3)", fileName);
        return 0;
    }
    fileIndex += 4;

    size_t chunkPosFmt = -1;
    size_t chunkPosData = -1;
    while (chunkPosFmt == -1 || chunkPosData == -1) {
        // get a chunk id and chunk size if we can
        if (fileData.size() < fileIndex + 8) {
            DBGOUT("[-----ERROR-----]%s is an invalid input file. (4)", fileName);
            return 0;
        }

        // get the chunk id if we can
        const unsigned char* chunkID = (unsigned char*)&fileData[fileIndex];
        fileIndex += 4;
        chunkSize = *(uint32_t*)&fileData[fileIndex];
        fileIndex += 4;

        //if we hit a fmt
        if (!memcmp(chunkID, "fmt ", 4)) {
            chunkPosFmt = (long)(fileIndex - 8);
        }
        //else if we hit a data
        else if (!memcmp(chunkID, "data", 4)) {
            chunkPosData = (long)(fileIndex - 8);
        }

        //skip to the next chunk
        fileIndex += chunkSize;
    }

    //we'll use this handy struct to load in 
    SMinimalWaveFileHeader waveData;

    //load the fmt part if we can
    fileIndex = chunkPosFmt;
    if (fileData.size() < fileIndex + 24) {
        DBGOUT("[-----ERROR-----]%s is an invalid input file. (5)", fileName);
        return 0;
    }
    memcpy(&waveData.m_subChunk1ID, &fileData[fileIndex], 24);
    fileIndex += 24;

    //load the data part if we can
    fileIndex = chunkPosData;
    if (fileData.size() < fileIndex + 8) {
        DBGOUT("[-----ERROR-----]%s is an invalid input file. (6)", fileName);
        return 0;
    }
    memcpy(&waveData.m_subChunk2ID, &fileData[fileIndex], 8);
    fileIndex += 8;

    //verify a couple things about the file data
    if (waveData.m_audioFormat != 1 ||
        waveData.m_numChannels < 1 ||
        waveData.m_numChannels > 2 ||
        waveData.m_bitsPerSample > 32 ||
        waveData.m_bitsPerSample % 8 != 0 ||
        waveData.m_blockAlign > 8) {
        DBGOUT("[-----ERROR-----]%s is an invalid input file. (7)", fileName);
        return 0;
    }

    //figure out how many samples and blocks there are total in the source data
    size_t bytesPerSample = waveData.m_blockAlign / waveData.m_numChannels;
    size_t numSourceSamples = waveData.m_subChunk2Size / bytesPerSample;

    //allocate space for the source samples
    data.resize(numSourceSamples);

    //read in the source samples at whatever sample rate / number of channels it might be in
    if (fileData.size() < fileIndex + numSourceSamples * bytesPerSample) {
        DBGOUT("[-----ERROR-----]%s is an invalid input file. (8)", fileName);
        return 0;
    }

    for (size_t nIndex = 0; nIndex < numSourceSamples; ++nIndex) {
        PCMToFloat(data[nIndex], &fileData[fileIndex], bytesPerSample);
        fileIndex += bytesPerSample;
    }

    //return our data
    numChannels = waveData.m_numChannels;
    sampleRate = waveData.m_sampleRate;
    numBytes = waveData.m_bitsPerSample / 8;

    DBGOUT("%s loaded - channels: %d, sampleRate: %d, bytesPerSample: %d, numBytes: %d",
        fileName, numChannels, sampleRate, bytesPerSample, numBytes);
    return numSourceSamples;
}

static float CubicHermite(float A, float B, float C, float D, float t)
{
    float a = -A / 2.0f + (3.0f*B) / 2.0f - (3.0f*C) / 2.0f + D / 2.0f;
    float b = A - (5.0f*B) / 2.0f + 2.0f*C - D / 2.0f;
    float c = -A / 2.0f + C / 2.0f;
    float d = B;

    return a*t*t*t + b*t*t + c*t + d;
}

inline float SampleChannelFractional(const std::vector<float>& input, float sampleFloat, uint16_t channel, uint16_t numChannels)
{
    size_t sample = size_t(sampleFloat);
    float sampleFraction = sampleFloat - std::floor(sampleFloat);

    size_t sampleIndexNeg1 = (sample > 0) ? sample - 1 : sample;
    size_t sampleIndex0 = sample;
    size_t sampleIndex1 = sample + 1;
    size_t sampleIndex2 = sample + 2;

    sampleIndexNeg1 = sampleIndexNeg1 * numChannels + channel;
    sampleIndex0 = sampleIndex0 * numChannels + channel;
    sampleIndex1 = sampleIndex1 * numChannels + channel;
    sampleIndex2 = sampleIndex2 * numChannels + channel;

    sampleIndexNeg1 = std::min(sampleIndexNeg1, input.size() - 1);
    sampleIndex0 = std::min(sampleIndex0, input.size() - 1);
    sampleIndex1 = std::min(sampleIndex1, input.size() - 1);
    sampleIndex2 = std::min(sampleIndex2, input.size() - 1);

    return CubicHermite(
        input[sampleIndexNeg1],
        input[sampleIndex0],
        input[sampleIndex1],
        input[sampleIndex2],
        sampleFraction);
}
