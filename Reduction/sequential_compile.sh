 #Module load Cuda Compilers and GCC
module load  cuda/5.0
module load  gcc/4.4.3

#compile
nvcc sequential.cu timer.c -o sequential
