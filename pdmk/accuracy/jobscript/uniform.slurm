#!/bin/bash

# set name of job
#SBATCH --job-name=dmk_accuracy_uniform
#SBATCH --output=job.%J.out
#SBATCH --error=job.%J.err

#SBATCH -p ccm
#SBATCH -C genoa
#SBATCH -N1 -n 1

module load gcc fftw openmpi

module list

export OMP_NUM_THREADS=64

cd /mnt/home/xgao1/work/PDMK4MC_Benchmark/pdmk/accuracy
julia --project=/mnt/home/xgao1/work/PDMK4MC_Benchmark/pdmk uniform_sys.jl