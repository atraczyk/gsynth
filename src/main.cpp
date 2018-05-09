#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <memory>

#include "app.h"
#include "pa_layer.h"
#include "soundfile.h"
#include "filter.h"
#include "Log.hpp"
#include "common.h"

#include "fft.h"
#include "pitchshift.h"

int main(int argc, char* argv[])
{
    App app(512, 512);
    app.execute(argc, argv);

#ifdef _WIN32
    system("pause");
#endif
    return 0;
}
