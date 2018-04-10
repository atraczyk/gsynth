# gsynth (WIP)

A real-time granular voice synthesizer designed for pitch shifting with visual spectrogram display, supporting three modes:

* sinusoidal re-synthesis
* inverse fft
* sample based 

### Installing

## Dependencies

# Linux

Ubuntu 16.04

```
sudo apt install portaudio19-dev libsdl2-dev
git clone https://atraczyk@bitbucket.org/atraczyk/gsynth.git && cd gsynth
make
```

# Windows
 
* download Portaudio source: http://www.portaudio.com/download.html
* download SDL2 binaries: https://www.libsdl.org/download-2.0.php
* extract them to gsynth\
* open gsynth.sln and build Debug/x64 configuration (should build portaudio.lib as a dependency)

## Built With

# Linux

Ubuntu 16.04 - gcc 7.2.0

# Windows

Visual Studio 2017

## TODO

* use grains shorter than the frame buffer
* overlap and crossfade delayed grains
* use low pass filter before computing fft
* try to low pass filter the fft by ignoring high-frequency bins
* only use hanning window and band limit when re-synthesizing
* sample based reconstruction (linear/cubic interpolation)
