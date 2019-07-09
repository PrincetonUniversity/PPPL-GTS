#!/bin/bash -l

#SBATCH -p debug  
#SBATCH -N 6
#SBATCH -t 00:30:00 
#SBATCH -J test  
#SBATCH --constraint=haswell
 
module load craype-hugepages2M

export OMP_NUM_THREADS=1

if [[ -z "$iloop" ]]
then
   iloop=0
fi

echo "iloop = $iloop"
if (("$iloop" <= 0))
then
   echo "Stop submit job!"
else
   echo "Submit another job!"
   let b=$iloop-1
   sbatch -d afterany:$SLURM_JOB_ID --export=iloop=$b submit_cori_loop.sh
fi

srun -n 192 ./GTS >& job.out

#running example
#sbatch --export=iloop=20 submit_cori_loop.sh

