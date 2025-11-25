#!/bin/bash

# set name of job
#SBATCH --job-name=fmmmmc_benchmark
#SBATCH --output=job.%J.out
#SBATCH --error=job.%J.err

#SBATCH -N1 -n 1
#SBATCH -C genoa
#SBATCH -p ccm
#SBATCH --time=48:00:00

module load python gcc fftw openmpi

source /mnt/home/xgao1/venv/fmmmc/bin/activate

module list

export OMP_NUM_THREADS=1
export MPLCONFIGDIR=$(pwd)/.cache

cd /mnt/home/xgao1/work/PDMK4MC_Benchmark/fmmmc/runtime
mpirun -n 1 python /mnt/home/xgao1/work/PDMK4MC_Benchmark/fmmmc/scripts/vary_number_nonuniform.py fmmmc_pbc_local pbc --local