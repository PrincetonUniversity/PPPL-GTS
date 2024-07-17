#!/bin/bash
#SBATCH -A pppl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --ntasks-per-socket=16
#SBATCH --cpus-per-task=4
#SBATCH --gpus-per-node=4
#SBATCH --time=00:30:00
#SBATCH --gpu-bind=none
#SBATCH --reservation=test

module load cudatoolkit/10.2 pgi/19.9 openmpi/pgi-20.4/4.0.4/64 hdf5/pgi-20.4/openmpi-4.0.4/1.10.6 fftw/gcc/openmpi-4.0.4/3.3.8
module list
export OMP_NUM_THREADS=1
# export PGI_ACC_TIME=1
#export PGI_ACC_DEBUG=1
#export PGI_ACC_NOTIFY=16
#export PGI_ACC_SYNCHRONOUS=1
mkdir -p snapshots_hdf5/AllPlanes
#pgprof -o output.%h.%p.%q{SLURM_PROCID} ./GTS >& job.out
export SLURM_CPU_BIND="cores"
srun ./MPS_wrapper.sh ./select_gpu_device.sh ./GTS >& job.out
# srun ./GTS >& job.out

