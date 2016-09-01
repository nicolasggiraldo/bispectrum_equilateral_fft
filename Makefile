CC=gcc
CFLAGS=-c -O3 -Wall -I/home/$(USER)/local/include/ -Iinclude
LFLAGS=-lm -L/home/$(USER)/local/lib -lfftw3 -lgsl -lgslcblas
SOURCES=main.c read_parameters.c readBinaryFile.c window_functions.c
#SOURCES=main.c read_parameters.c window_functions.c
OBJECTS=$(SOURCES:%.c=lib/%.o)
PROGRAM=bispectrum_equilateral_fft

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@.x $(LFLAGS)

lib/%.o:src/%.c
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm ./*/*~ *~ ./*/*# *# lib/*.o *.x
