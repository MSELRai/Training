#!/bin/bash
#SBATCH --job-name=MYCP2K
#SBATCH --output=MYCP2K.o%j
#SBATCH --error=MYCP2K.e%j
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --partition=wholenode
#SBATCH --time=00:30:00
#SBATCH --account=##########
#SBATCH --exclusive
ulimit -s unlimited
export OMP_NUM_THREADS=1
echo "SLURM_NTASKS: " $SLURM_NTASKS

exe=/anvil/projects/x-che140141/Programs/cp2k/2022.1/exe/Linux-gnu-x86_64/cp2k.psmp

mpirun -np $SLURM_NTASKS $exe -i cp2k.inp -o cp2k.out
