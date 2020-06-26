#!/bin/sh
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=20
#SBATCH --mem=64000
#SBATCH --time=03:30:00

#source ~/.bashrc
source /home/mjh/SETUP_CONDA.sh
time mpiexec nemoMass MFMF_SOSim_3freq_tiles.yml -M
