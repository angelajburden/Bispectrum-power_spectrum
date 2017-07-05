# Makefile for pwer spectrum
CC = g++
CFLAGS = -Wall
PKG_CONFIG_PATH=/gpfs/apps/hpc/Libs/FFTW/2.1.5
export PKG_CONFIG_PATH

INCLUDE = -I$(PKG_CONFIG_PATH)/include
DEPS = header.h
OBJ = calc_power.o power_functions.o util.o integration.o io_gadget_slim.o plotting_functions.o
LIBFFTW = -L$(PKG_CONFIG_PATH)/lib
LIBS = $(LIBFFTW) $(INCLUDE) -lrfftw -lfftw -lm

%.o: %.c $(DEPS)
	$(CC)  -c -o $@ $< $(CFLAGS) $(LIBS) 

power_sp: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)
clean:
	\rm calc_power.o power_functions.o util.o integration.o io_gadget_slim.o plotting_functions.o power_sp