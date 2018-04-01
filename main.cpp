#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "pa_layer.h"
#include "Log.hpp"

#include "kiss_fft.h"

#ifndef M_PI
#define M_PI 3.14159265358979324
#endif

#define SR 22050
#define N (int)(SR*0.025)
#define FR1 880
#define FR2 666

void TestFft(const char* title, const kiss_fft_cpx in[N], kiss_fft_cpx out[N])
{
    kiss_fft_cfg cfg;

    printf("%s\n", title);

    if ((cfg = kiss_fft_alloc(N, 0/*is_inverse_fft*/, NULL, NULL)) != NULL)
    {
        size_t i;

        kiss_fft(cfg, in, out);
        free(cfg);

        /*for (i = 0; i < N; i++)
            printf(" in[%02i]=%+f, %+f  out[%02i]=%+f, %+f M[%02i]=%+f\n",
            i, in[i].r, in[i].i,
            i, out[i].r, out[i].i,
            i, sqrt((out[i].r * out[i].r) + (out[i].i * out[i].i)));*/
    }
    else
    {
        printf("not enough memory?\n");
        exit(-1);
    }
}

int main(void)
{
    /*PortAudioLayer pa;

    pa.startStream();
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    pa.stopStream();*/

    kiss_fft_cpx in[N], out[N];
    size_t i;

    for (i = 0; i < N; i++)
        in[i].r = in[i].i = 0;
    TestFft("Zeroes (complex)", in, out);

    for (i = 0; i < N; i++)
        in[i].r = 1, in[i].i = 0;
    TestFft("Ones (complex)", in, out);

    auto elapsedTimePerFrame = 1.0 / static_cast<double>(SR);
    for (i = 0; i < N; i++) {
        auto t = i * elapsedTimePerFrame;
        in[i].r = sin(2 * M_PI * t * FR1);
        //in[i].r = 0.75 * sin(2 * M_PI * t * FR1) + 0.25 * sin(2 * M_PI * t * FR2);
        in[i].i = 0;
    }
    TestFft("SineWave (complex)", in, out);
    double mag[N];
    double maxMag = 0;
    int maxIndex = 0;
    for (i = 0; i < N / 2; i++) {
        mag[i] = sqrt((out[i].r * out[i].r) + (out[i].i * out[i].i));
        if (mag[i] >= maxMag) {
            maxMag = mag[i];
            maxIndex = i;
        }
        //DBGOUT("i: %d, m: %0.4f", i, mag[i]);
    }
    DBGOUT("maxi: %d, maxm: %0.4f", maxIndex, maxMag);
    // weight the index
    auto width = 3;
    auto fundValue = static_cast<double>(maxIndex);
    for (i = 1; i <= width; i++) {
        if (maxIndex - i >= 0)
            fundValue = fundValue - (mag[maxIndex - i] / mag[maxIndex]);
        if (maxIndex + i < N)
            fundValue = fundValue + (mag[maxIndex + i] / mag[maxIndex]);
    }
    auto fundFreq = (fundValue / N) * SR;
    DBGOUT("fundamental freq: %0.4f", fundFreq);
    //DBGOUT("fundamental freq: %0.4f, diff: %0.4f", fundFreq, fundFreq - FR1);

    system("pause");
    return 0;
}
