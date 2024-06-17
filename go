#!/bin/bash

module load openmpi/3.1.2-pgi-20.7
module load pgi/20.7
export PGI_LIBS=/export/apps/CentOS7/pgi/2020/Linux_x86_64/20.7/
module load cmake

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/export/zimmerman/ndmeier/SlaterGPU/src/libio
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/export/zimmerman/ndmeier/SlaterGPU/src/libcintw

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/export/zimmerman/khoidang/eris/libcint/build_zimintel_pgi/libcint/lib64/
export LIBCINT_PATH=~khoidang/eris/libcint/build_pgi/libcint/

export OMP_NUM_THREADS=2

alias r="~/integrate/gpu/zest.exe"
alias r2="~/integrate/gpu/zest.exe2"
alias mpirun="mpirun --bind-to none --mca btl ^openib"
#alias m="mpirun -n 1 --bind-to socket --mca btl ^openib /home/paulzim/integrate/test.exe "
