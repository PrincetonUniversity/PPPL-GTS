#!/bin/bash -l

#SBATCH -p debug  

#SBATCH -N 24

#SBATCH -t 00:30:00 

#SBATCH -J test  

#SBATCH --constraint=haswell
 
module load craype-hugepages2M

export OMP_NUM_THREADS=1

srun -n 768 ./GTS >& job.out
