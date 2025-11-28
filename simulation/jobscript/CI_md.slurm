#!/bin/bash

# set name of job
#SBATCH --job-name=charge_inv_md
#SBATCH --output=job.%J.out
#SBATCH --error=job.%J.err

#SBATCH -p ccm
#SBATCH -C genoa
#SBATCH -N1

module load gcc fftw openmpi

module list

export OMP_NUM_THREADS=64

cd /mnt/home/xgao1/work/PDMK4MC_Benchmark/simulation/colloid
julia -t auto --project=/mnt/home/xgao1/work/PDMK4MC_Benchmark/simulation md.jl