#!/bin/sh
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=20
#SBATCH --mem=64000
#SBATCH --time=23:59:00

source ~/.bashrc
time mpiexec nemoCosmo mocks_small/mockCatalog_1.fits MFMF_SOSim_3freq_small/selFn

