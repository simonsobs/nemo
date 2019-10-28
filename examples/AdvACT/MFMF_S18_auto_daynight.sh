#!/bin/sh
#SBATCH --nodes=18
#SBATCH --ntasks-per-node=16
#SBATCH --mem=64000
#SBATCH --time=23:59:00

source ~/.bashrc
time mpiexec nemo MFMF_S18_auto_daynight.yml -M -n
