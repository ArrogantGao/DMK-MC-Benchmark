# PDMK4MC_Benchmark

This repository is used to benchmark the performance of PDMK4MC's C++ implementation at [PDMK4MC](https://github.com/ArrogantGao/PDMK4MC), which is wrapped in Julia.

The results are organized as follows:
- `fmmmc/` is for benchmarking the FMM-based MC algorithm in https://doi.org/10.1016/j.jcp.2020.110099, with their own code. I slightly modified the number of particles and the distribution.
- `pdmk/` is for benchmarking our method, `pdmk/accuracy/` includes the accuracy test and the corresponding runtime results, and `runtime/` is purely for runtime comparison with that of `fmmmc/`.
- `simulation/` contains the simulation code for the colloid system and the NaCl system, with the `Molly.jl` as backend of both MC and MD simulations, `Molly.jl` is slightly modified with the scripts in `MollyExt/`.
- `plots/` contains the scripts to plot the results.

## fmmmc

One have to create a python environment with the requirements.txt in `fmmmc/requirements.txt`. The core deps are [`ppmd`](https://github.com/ppmd/ppmd), [`coulomb-mc`](https://github.com/ppmd/coulomb_mc), and [`coulomb-kmc`](https://github.com/ppmd/coulomb_kmc).

I only used the scripts `fmmmc/scripts/vary_number_uniform.py` and `fmmmc/scripts/vary_number_nonuniform.py` to generate the runtime results.

## pdmk

One have to clone the PDMK4MC repo and compile the C++ code locally, then in julia REPL, run `] dev /path/to/PDMK4MC/julia/` to add the package to the Julia environment. Also, notice to change the `MPIPreferences` to `OpenMPI` to be compatible with the C++ code.

Then install `ParticleMeshEwald.jl` by running `] add https://github.com/flatironinstitute/ParticleMeshEwald.jl` in the julia REPL, it is a lightweight implementation of the Particle Mesh Ewald method by me.

## simulation

Similarly, one have to manually install `PDMK4MC` and `ParticleMeshEwald.jl` in the simulation code.