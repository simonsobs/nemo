#!/bin/sh
#SBATCH --nodes=18
#SBATCH --ntasks-per-node=16
#SBATCH --mem=64000
#SBATCH --time=03:59:00

source ~/.bashrc
time mpiexec nemo S18d_202003.yml -M -n
