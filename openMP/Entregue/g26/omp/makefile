CC = gcc
CFLAGS = -lm -fopenmp -O2

all: parsim 

parsim: parsim-omp.c
	$(CC) -o parsim-omp parsim-omp.c $(CFLAGS)

