#include "app.h"
#include "Log.hpp"

#ifdef _WIN32
#include <conio.h>
#else
#include <ncurses.h>
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
    auto _fftDataBlob = paLayer.getFftData();

    SDL_GL_MakeCurrent(window_, glContext_);

    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    glLineWidth(2.5);
    glColor4f(0.1f, 0.2f, 0.7f, 1.0f);

    glBegin(GL_LINES);

    bool linearPlot = false;
    float h;
    auto barSize = 2.0f / _fftDataBlob.size();
    for (int i = 0; i < _fftDataBlob.size(); i++) {
        if (linearPlot) {
            h = _fftDataBlob.at(i).amplitude * 2.0f;
        } else {
            auto dB = 20 * log10(_fftDataBlob.at(i).amplitude);
            h = (1.0f - std::min(1.0f, dB / -60.0f)) * 2.0f;
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
    
    float pitchShift = 0.0;

    while (running_) {
        while (SDL_PollEvent(&event) != 0) {
            onEvent(&event);

            switch( event.type ) {
                case SDL_QUIT:
                    running_ = false;
                    break;
                case SDL_KEYDOWN:
                case SDL_KEYUP:
                    switch( event.key.keysym.sym ){
                        case SDLK_ESCAPE:
                            running_ = false;
                            break;
                        case SDLK_UP:
                            pitchShift += 1;
                            std::cout << "pitchShift " << pitchShift << "\n";
                            paLayer.setPitchShift(pitchShift);
                            break;
                        case SDLK_DOWN:
                            pitchShift -= 1;
                            std::cout << "pitchShift " << pitchShift << "\n";
                            paLayer.setPitchShift(pitchShift);
                            break;
                        default:
                            break;
                    }
                default:
                    break;
            }
        }

        render();

        SDL_Delay(1);
    }

    cleanup();
    return;
}
