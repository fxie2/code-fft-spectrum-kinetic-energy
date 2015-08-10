#!/bin/bash
if [ "$1" = "s"  ] 
then
   h5fc -O3 fft.f90  main_ke_mem.f90 -o fftke
elif [ "$1" = "omp" ] 
then
   h5fc -O3 -fopenmp fft.f90 main_ke_mem.f90 -o fftke
fi
   
