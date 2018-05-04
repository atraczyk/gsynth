cmake_minimum_required(VERSION 3.10)
project(gsynth)

set(CMAKE_CXX_STANDARD 11)

include_directories(src)
include_directories(src/audio)
include_directories(src/audio/portaudio)
include_directories(src/fft)
include_directories(src/fft/kiss_fft)

add_executable(gsynth
        src/audio/portaudio/pa_layer.cpp
        src/audio/portaudio/pa_layer.h
        src/audio/audiolayer.cpp
        src/audio/audiolayer.h
        src/audio/ringbuffer.h
        src/audio/soundfile.cpp
        src/audio/soundfile.h
        src/fft/kiss_fft/kiss_fft.c
        src/fft/kiss_fft/kiss_fft.h
        src/fft/kiss_fft/kiss_fft_guts.h
        src/fft/fft.cpp
        src/fft/fft.h
        src/app.cpp
        src/app.h
        src/common.h
        src/filter.cpp
        src/filter.h
        src/Log.hpp
        src/main.cpp)
