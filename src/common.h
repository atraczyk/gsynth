#pragma once

#include <cmath>

#ifndef M_PI
#define M_PI   3.14159265358979323846264338327950288
#endif
#ifndef M_PI_2
#define M_PI_2  1.57079632679489661923132169163975144
#endif
#ifndef M_2PI
#define M_2PI   6.28318530717958647692528676655900576
#endif

enum SynthType {
    resynth = 0,
    ifft = 1,
    sample = 2
};

#define SYNTHTYPE SynthType::resynth
