#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "app.h"
#include "pa_layer.h"
#include "soundfile.h"
#include "Log.hpp"

#ifdef _WIN32
#define NOMINMAX
#include <Objbase.h> // CoInitialize/CoUninitialize
#endif // _WIN32

int main(int argc, char* argv[])
{
    // COM init not done by Portaudio?
#ifdef _WIN32
    CoInitialize(0);
#endif // _WIN32

    App app(512, 512);
    app.execute(argc, argv);

    uint16_t numChannels;
    uint32_t sampleRate;
    uint16_t numBytes;
    std::vector<float> source, out;
    ReadWaveFile("legend1.wav", source, numChannels, sampleRate, numBytes);

    //system("pause");

#ifdef _WIN32
    CoUninitialize();
#endif // _WIN32
    return 0;
}
