IDIRS=-I../src -Isrc -Isrc/audio -Isrc/audio/portaudio -Isrc/fft -Isrc/fft/kiss_fft
CC=g++
CFLAGS=-I$(IDIRS) -std=c++11 -D_DEBUGPRINTF

LIBS=-lSDL2 -lm -lpthread -lportaudio -lGL -lncurses

gsynth:
	$(CC) -o gsynth \
	src/main.cpp \
	src/app.cpp \
	src/filter.cpp \
	src/audio/audiolayer.cpp \
	src/audio/soundfile.cpp \
	src/audio/portaudio/pa_layer.cpp \
	src/fft/fft.cpp \
	src/fft/kiss_fft/kiss_fft.c \
	$(CFLAGS) $(LIBS)

clean:
	-@rm -rf *.o && rm gsynth 2>/dev/null || true
