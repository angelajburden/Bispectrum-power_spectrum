# Makefile for bispectrum
CC = g++
CFLAGS = -Wall
PKG_CONFIG_PATH=/gpfs/apps/hpc/Libs/FFTW/2.1.5
export PKG_CONFIG_PATH

INCLUDE = -I$(PKG_CONFIG_PATH)/include
DEPS = header.h
OBJ = calc_bispectrum.o bispectrum_functions.o util.o integration.o io_gadget_slim.o plotting_functions.o
LIBFFTW = -L$(PKG_CONFIG_PATH)/lib
LIBS = $(LIBFFTW) $(INCLUDE) -lrfftw -lfftw -lm

%.o: %.c $(DEPS)
	$(CC)  -c -o $@ $< $(CFLAGS) $(LIBS) 

bi_sp3: $(OBJ)
	$(CC) -g -o $@ $^ $(CFLAGS) $(LIBS)
clean:
	\rm calc_bispectrum.o bispectrum_functions.o util.o integration.o io_gadget_slim.o plotting_functions.o bi_sp3