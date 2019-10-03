#!/bin/sh
#SBATCH --nodes=7
#SBATCH --ntasks-per-node=20
#SBATCH --mem=64000
#SBATCH --time=23:59:00

source ~/.bashrc
time mpiexec nemo PS_S16_f150_auto.yml -M
