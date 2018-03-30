#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "pa_layer.h"
#include "Log.hpp"

int main(void)
{
    PortAudioLayer pa;

    pa.startStream();
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    pa.stopStream();

    system("pause");
    return 0;
}
