#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=64000
#SBATCH --time=23:59:00

source ~/.bashrc
srun nemoCosmo MFMF_SOSim_3freq_tiles/selFn/config.yml mocks_tiles/mockCatalog_3.fits MFMF_SOSim_3freq_tiles/selFn 1000 20 -o cosmoResults_tiles/mockCatalog3

