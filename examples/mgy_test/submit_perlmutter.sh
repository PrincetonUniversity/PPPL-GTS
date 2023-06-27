#!/bin/bash
#SBATCH -A mp19
#SBATCH -C gpu
#SBATCH -q debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=2
#SBATCH --gpus-per-node=4
#SBATCH --time=00:30:00
#SBATCH --gpu-bind=none


module list
export OMP_NUM_THREADS=1
#export PGI_ACC_TIME=1
#export PGI_ACC_DEBUG=1
#export PGI_ACC_NOTIFY=16
#export PGI_ACC_SYNCHRONOUS=1
mkdir -p snapshots_hdf5/AllPlanes
#pgprof -o output.%h.%p.%q{SLURM_PROCID} ./GTS >& job.out
export SLURM_CPU_BIND="cores"
srun ./GTS >& job.out
#pgprof ./GTS >& job.out

