CC = gcc
CFLAGS = -lm -fopenmp -O2 
all: parsim_mpi 

parsim: parsim_mpi.c
	$(CC) -o parsim_mpi parsim_mpi.c $(CFLAGS)

