#!/bin/bash

module load openmpi/3.1.2-pgi-20.7
module load pgi/20.7
module load cmake

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/khoi/Modules/libraries/libcint/gcc-8.2.0/lib64
export LIBCINT_PATH=/home/khoi/Modules/libraries/libcint/gcc-8.2.0

export OMP_NUM_THREADS=2

