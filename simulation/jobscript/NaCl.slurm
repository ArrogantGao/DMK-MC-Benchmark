#!/bin/bash

# set name of job
#SBATCH --job-name=NaCl_mc
#SBATCH --output=job.%J.out
#SBATCH --error=job.%J.err

#SBATCH -p ccm
#SBATCH -N1

module load gcc fftw openmpi

module list

export OMP_NUM_THREADS=64

cd /mnt/home/xgao1/work/PDMK4MC_Benchmark/simulation/1_1_salt
julia -t auto --project=/mnt/home/xgao1/work/PDMK4MC_Benchmark/simulation NaCl_mc.jl > NaCl_mc.log 2>&1