#pragma once

#ifdef _WIN32
#include <SDL.h>
#include <SDL_opengl.h>
#else
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#endif // _WIN32

#undef main

#include "pa_layer.h"

#define SAMPLERATE 48000

class App {
private:
    App();

    void onEvent(SDL_Event* Event);
    bool initVideo();
    void render();
    void cleanup();

    bool running_ = true;

    SDL_Window* window_ = NULL;
    SDL_GLContext glContext_;

    int windowWidth_;
    int windowHeight_;

    PortAudioLayer paLayer;

public:
    App(int windowWidth = 512, int windowHeight = 512);

    void execute(int argc, char* argv[]);
};
