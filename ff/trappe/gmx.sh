#!/bin/bash
#SBATCH --job-name=q6
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00-06

module load apps/gromacs-2019.3
export OMP_NUM_THREADS=8

srun -n 1 gmx_mpi grompp -p trappe -f nvt_qua_steep -c qua_ow_400_12000_h60 -o qua_ow_400_12000_h60
srun -n 1 gmx_mpi mdrun -s -o -x -c -e -g -v -deffnm qua_ow_400_12000_h60
rm ./*pdb
srun -n 1 gmx_mpi grompp -p trappe -f nvt_qua_short -c qua_ow_400_12000_h60 -o qua_ow_400_12000_h60
srun -n 1 gmx_mpi mdrun -s -o -x -c -e -g -v -deffnm qua_ow_400_12000_h60
srun -n 1 gmx_mpi grompp -p trappe -f nvt_qua -c qua_ow_400_12000_h60 -o qua_ow_400_12000_h60
srun -n 1 gmx_mpi mdrun -s -o -x -c -e -g -v -deffnm qua_ow_400_12000_h60
rm ./#*#
