CC = mpicc
CFLAGS = -lm -fopenmp -O2 
all: parsim_mpi 

parsim: parsim_mpi.c
	$(CC) -o parsim-mpi parsim-mpi.c $(CFLAGS)

