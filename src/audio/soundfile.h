#pragma once

#include <vector>

//this struct is the minimal required header data for a wav file
struct SMinimalWaveFileHeader
{
    //the main chunk
    unsigned char   m_chunkID[4];
    uint32_t        m_chunkSize;
    unsigned char   m_format[4];

    //sub chunk 1 "fmt "
    unsigned char   m_subChunk1ID[4];
    uint32_t        m_subChunk1Size;
    uint16_t        m_audioFormat;
    uint16_t        m_numChannels;
    uint32_t        m_sampleRate;
    uint32_t        m_byteRate;
    uint16_t        m_blockAlign;
    uint16_t        m_bitsPerSample;

    //sub chunk 2 "data"
    unsigned char   m_subChunk2ID[4];
    uint32_t        m_subChunk2Size;

    //then comes the data!
};

enum class ECrossFade
{
    None,
    In,
    Out,
};

inline void FloatToPCM(unsigned char *PCM, const float& in, size_t numBytes);
inline void PCMToFloat(float& out, const unsigned char *PCM, size_t numBytes);
bool WriteWaveFile(const char *fileName, std::vector<float>& dataFloat, uint16_t numChannels, uint32_t sampleRate, uint16_t numBytes);
bool ReadFileIntoMemory(const char *fileName, std::vector<unsigned char>& data);
uint64_t ReadWaveFile(const char *fileName, std::vector<float>& data, uint16_t& numChannels, uint32_t& sampleRate, uint16_t& numBytes);
static float CubicHermite(float A, float B, float C, float D, float t);
inline float SampleChannelFractional(const std::vector<float>& input, float sampleFloat, uint16_t channel, uint16_t numChannels);