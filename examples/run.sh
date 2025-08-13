#!/bin/bash
#SBATCH -N 1 --ntasks-per-node=8 --gres=gpu:1 --mem=60gb --time=2-00:00:00
#SBATCH -e log.err -o log.out
#SBATCH -p zimA10 --job-name=sgpu-test

##* Move this file to input directory and submit to SLURM

# source ../../env.set.local0
eval "$(pixi shell-hook)"
export OMP_NUM_THREADS=2

time sgpu.exe > hf.out
