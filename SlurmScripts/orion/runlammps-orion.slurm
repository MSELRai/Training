#!/bin/bash
#SBATCH --job-name=Ar
#SBATCH --output=Ar.o%j
#SBATCH --error=Ar.e%j
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --partition=400p48h
#SBATCH --qos=funded
#SBATCH --time=01:00:00
#SBATCH --mail-user=USER@msstate.edu
#SBATCH --mail-type=ALL
#SBATCH --account=######-######
#SBATCH --exclusive
ulimit -s unlimited
export OMP_NUM_THREADS=1
echo "SLURM_NTASKS: " $SLURM_NTASKS

module purge
module load intel/2021.2
module load impi/2021.2
module load contrib
module load lammps-msel/29Sep2021

date 
echo "My LAMMPS Simulation"
srun lmp -in in.lammps
date

