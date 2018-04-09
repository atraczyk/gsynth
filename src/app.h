#pragma once

#include <SDL2/SDL.h>
#undef main
#include <SDL2/SDL_opengl.h>

#include "pa_layer.h"

#define SAMPLERATE 44100

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
