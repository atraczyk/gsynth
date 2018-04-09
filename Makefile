IDIRS=-I../src -Isrc -Isrc/audio -Isrc/audio/portaudio -Isrc/fft -Isrc/fft/kiss_fft
CC=g++
CFLAGS=-I$(IDIRS)

LIBS=-lSDL2main -lSDL2 -lm -lpthread -lportaudio -lGL

gsynth:
	$(CC) -o gsynth \
	src/main.cpp \
	src/app.cpp \
	src/audio/ringbuffer.cpp \
	src/audio/audiolayer.cpp \
	src/audio/portaudio/pa_layer.cpp \
	src/fft/fft.cpp \
	src/fft/kiss_fft/kiss_fft.c \
	$(CFLAGS) $(LIBS)

clean:
	-@rm -rf *.o && rm gsynth 2>/dev/null || true
