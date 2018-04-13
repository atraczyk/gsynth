# GSYNTH (WIP)

A real-time granular synthesizer / phase vocoder designed for pitch shifting.
Will support three re-synthesis modes:

* sinusoidal additive (computationally expensive / simple processing)
* inverse fft (fast / difficult processing)
* sample based reconstruction

Also outputs frequency spectrum to SDL2 surface.

## Building

```
git clone https://bitbucket.org/atraczyk/gsynth.git
```

### Linux

Debian/Ubuntu deps
```
sudo apt install portaudio19-dev libsdl2-dev libncurses-dev
```

Fedora deps
```
sudo dnf install portaudio-devel SDL2-devel ncurses-devel
```

then
```
make
```

### Windows
 
* download Portaudio source: http://www.portaudio.com/download.html
* download SDL2 binaries: https://www.libsdl.org/download-2.0.php
* extract them to gsynth\
* open gsynth.sln and build Debug/x64 configuration (should build portaudio.lib as a dependency)

## Tested builds With

### Linux

Ubuntu 16.04 	- gcc 7.2.0
Fedora 27 		- gcc 7.3.1

### Windows

Windows 10 1709 (16299.371) Visual Studio 2017

## TODO

* export wav files and use gtkglext ui for debug
* implement phase incoherence adjustment
* use grains shorter than the frame buffer
* overlap and add (crossfade) delayed grains
* sample based reconstruction (linear/cubic interpolation)
