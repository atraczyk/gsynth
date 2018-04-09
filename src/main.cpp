#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "app.h"
#include "pa_layer.h"
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

    //system("pause");

#ifdef _WIN32
    CoUninitialize();
#endif // _WIN32
    return 0;
}
