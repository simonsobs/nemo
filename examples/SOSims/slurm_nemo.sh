#!/bin/sh
#SBATCH --nodes=27
#SBATCH --ntasks-per-node=20
#SBATCH --mem=64000
#SBATCH --time=03:30:00

source ~/.bashrc
time mpiexec nemo MFMF_SOSim_3freq_tiles.yml -M
