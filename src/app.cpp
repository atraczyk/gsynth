#include "app.h"
#include "Log.hpp"

#ifdef _WIN32
#include <conio.h>
#else
#include <unistd.h>
#include <termios.h>
char getch()
{
    char buf=0;
    struct termios old={0};
    fflush(stdout);
    if(tcgetattr(0, &old)<0)
        perror("tcsetattr()");
    old.c_lflag&=~ICANON;
    old.c_lflag&=~ECHO;
    old.c_cc[VMIN]=1;
    old.c_cc[VTIME]=0;
    if(tcsetattr(0, TCSANOW, &old)<0)
        perror("tcsetattr ICANON");
    if(read(0,&buf,1)<0)
        perror("read()");
    old.c_lflag|=ICANON;
    old.c_lflag|=ECHO;
    if(tcsetattr(0, TCSADRAIN, &old)<0)
        perror ("tcsetattr ~ICANON");
    printf("%c\n",buf);
    return buf;
 }
#endif

App::App()
{
}

App::App(int windowWidth, int windowHeight)
    : windowWidth_(windowWidth)
    , windowHeight_(windowHeight)
    , paLayer(AudioFormat(SAMPLERATE, 2), AudioFormat(SAMPLERATE, 1))
{
}

bool App::initVideo()
{
    if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_EVENTS) != 0) {
        SDL_Log("Failed to initialize SDL: %s", SDL_GetError());
        return false;
    }

    if ((window_ = SDL_CreateWindow(
        "gsynth - spectrogram",
        SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
        windowWidth_, windowHeight_, SDL_WINDOW_OPENGL)
        ) == nullptr) {
        SDL_Log("Unable to create SDL Window: %s", SDL_GetError());
        return false;
    }

    glContext_ = SDL_GL_CreateContext(window_);

    if (SDL_GL_SetSwapInterval(1) < 0) {
        SDL_Log("Warning: Unable to set VSync! SDL Error: %s", SDL_GetError());
    }

    return true;
}

void App::render()
{
    auto frequencyData = paLayer.getFrequencyData();
    //DBGOUT("FREQ.SIZE: %d", frequencyData.size());

    SDL_GL_MakeCurrent(window_, glContext_);

    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    glLineWidth(2.5);
    glColor4f(0.1f, 0.2f, 0.7f, 1.0f);

    glBegin(GL_LINES);

    bool linearPlot = false;
    double h;
    auto barSize = 2.0f / frequencyData.size();
    for (int i = 0; i < frequencyData.size(); i++) {
        if (linearPlot) {
            h = frequencyData.at(i).second * 2.0f;
        } else {
            auto dB = 20 * log10(frequencyData.at(i).second);
            h = (1.0f - std::min(1.0, dB / -60.0f)) * 2.0f;
        }
        auto x = -1.0f + barSize + barSize * i;
        glVertex3f(x, -1.0f, 0.0f);
        glVertex3f(x, -1.0f + h, 0.0f);
    }

    glEnd();

    SDL_GL_SwapWindow(window_);
}

void App::cleanup()
{
    SDL_GL_DeleteContext(glContext_);

    if (window_) {
        SDL_DestroyWindow(window_);
        window_ = nullptr;
    }

    SDL_Quit();

    paLayer.stopStream();
}

void App::onEvent(SDL_Event* event)
{
    if (event->type == SDL_APP_WILLENTERBACKGROUND ||
        event->type == SDL_KEYUP && event->key.keysym.sym == SDLK_ESCAPE) {
        running_ = false;
    }
}

void App::execute(int argc, char* argv[])
{
    paLayer.startStream();

    auto videoInitialized = initVideo();

    if (!videoInitialized) {
        getch();
        cleanup();
        return;
    }

    SDL_Event event;

    while (running_) {
        while (SDL_PollEvent(&event) != 0) {
            onEvent(&event);

            if (event.type == SDL_QUIT) {
                running_ = false;
            }
        }

        render();

        SDL_Delay(1);
    }

    cleanup();
    return;
}
