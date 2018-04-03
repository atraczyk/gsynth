#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "pa_layer.h"
#include "Log.hpp"

#include "fft.h"

#define SR 44100

void test()
{
    auto out_format = AudioFormat(SR, 2);
    auto in_format = AudioFormat(SR, 1);
    PortAudioLayer pa(out_format, in_format);

    pa.startStream();
    std::this_thread::sleep_for(std::chrono::milliseconds(20000));
    pa.stopStream();

}

int main(void)
{
    test();

    system("pause");
    return 0;
}
