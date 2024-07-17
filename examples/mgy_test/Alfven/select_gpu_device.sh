#!/bin/bash
# select_gpu_device wrapper script
num_GPUs_per_Node=4
echo "num_GPUs_per_Nodes is set to be $num_GPUs_per_Node" 
export CUDA_VISIBLE_DEVICES=$(( SLURM_LOCALID % num_GPUs_per_Node ))
exec $*