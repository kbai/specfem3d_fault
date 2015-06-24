#!/bin/bash
#./configure FC=ifort MPIFC=mpiifort  --with-cuda
./configure FC=ifort MPIFC=mpiifort --with-cuda 
make -j 9  
