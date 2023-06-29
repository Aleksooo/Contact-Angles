#!/bin/bash
#SBATCH --job-name=l5g19
#SBATCH --partition=gpu_debug
#SBATCH --nodes=1
#SBATCH --gpus=1
#SBATCH --cpus-per-task=8
module load apps/gromacs-2019.3
export OMP_NUM_THREADS=8

srun -n 1 gmx_mpi grompp -p gromos -f nvt_steep -c w_4162 -o w_4162
srun -n 1 gmx_mpi mdrun -s -o -x -c -e -g -v -deffnm w_4162
srun -n 1 gmx_mpi grompp -p gromos -f nvt_short -c w_4162 -o w_4162
srun -n 1 gmx_mpi mdrun -s -o -x -c -e -g -v -deffnm w_4162
