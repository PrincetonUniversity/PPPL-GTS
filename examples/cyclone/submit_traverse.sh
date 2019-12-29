#!/bin/bash
#SBATCH -A pppl
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=32
#SBATCH --ntasks-per-socket=16
#SBATCH --cpus-per-task=4
#SBATCH --gpus-per-node=4
#SBATCH --time=24:00:00

export OMP_NUM_THREADS=1
export PGI_ACC_POOL_THRESHOLD=100
#export PGI_ACC_TIME=1
#export PGI_ACC_DEBUG=1
#export PGI_ACC_MEM_MANAGE=0
#export PGI_ACC_POOL_ALLOC=0
#export PGI_ACC_NOTIFY=2

srun ./GTS >& job.out
