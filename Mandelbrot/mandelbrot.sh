#!/bin/bash
#$ -q class8i
#$ -pe mpi 16
#$ -N mandelbrot_ms
#$ -R y

date
hostname
echo -e "\n\n"

# Module load boost
module load boost/1.57.0

# Module load OpenMPI
module load openmpi-1.8.3/gcc-4.9.2

# Run
mpirun -np 8 ./mandelbrot_ms 1000 1000
