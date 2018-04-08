#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "pa_layer.h"
#include "Log.hpp"

#include <SDL.h>
#undef main
#include <SDL_opengl.h>

#define SAMPLERATE 44100

static bool quitting = false;
static float r = 0.0f;
static SDL_Window *window = NULL;
static SDL_GLContext gl_context;

void render()
{

    SDL_GL_MakeCurrent(window, gl_context);

    r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

    glClearColor(r, 0.4f, 0.1f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    SDL_GL_SwapWindow(window);

}

int SDLCALL watch(void *userdata, SDL_Event* event)
{

    if (event->type == SDL_APP_WILLENTERBACKGROUND) {
        quitting = true;
    }

    return 1;
}

int main(void)
{
    auto out_format = AudioFormat(SAMPLERATE, 2);
    auto in_format = AudioFormat(SAMPLERATE, 1);
    PortAudioLayer pa(out_format, in_format);

    pa.startStream();

    //std::this_thread::sleep_for(std::chrono::milliseconds(20000));

    if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_EVENTS) != 0) {
        SDL_Log("Failed to initialize SDL: %s", SDL_GetError());
        return 1;
    }

    window = SDL_CreateWindow("gsynth - spectrogram", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, 512, 512, SDL_WINDOW_OPENGL);

    gl_context = SDL_GL_CreateContext(window);

    SDL_AddEventWatch(watch, NULL);


    while (!quitting) {

        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                quitting = true;
            }
        }

        render();
        SDL_Delay(2);

    }

    SDL_DelEventWatch(watch, NULL);
    SDL_GL_DeleteContext(gl_context);
    SDL_DestroyWindow(window);
    SDL_Quit();

    pa.stopStream();

    system("pause");
    return 0;
}
