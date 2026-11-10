#!/bin/bash
#SBATCH -J gts_rand
#SBATCH --nodes=1
#SBATCH --ntasks=32                 # must equal mzetamax (8) x npartdom (4) in INPUT.d
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=03:00:00
##SBATCH -A <account>

# Example Slurm script.
# Load the modules that were used to build GTS, for example:
#   module load <compiler> <mpi> <hdf5> <netcdf> <fftw>
# and copy (or link) the GTS executable into this directory as ./GTS.

export OMP_NUM_THREADS=1
mkdir -p snapshots_hdf5/AllPlanes

srun ./GTS >& job.out
