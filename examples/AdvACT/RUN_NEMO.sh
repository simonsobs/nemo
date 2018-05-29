#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=64000
#SBATCH --time=12:00:00

source ~/.bashrc
time mpirun nemo MF_AdvACT_multiScale_tileDeck_hybrid.par
