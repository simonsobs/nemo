#!/bin/sh
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=20
#SBATCH --mem=64000
#SBATCH --time=23:59:00

source /home/mjh/SETUP_CONDA.sh
time mpiexec nemoCosmo mocks_small/mockCatalog_1.fits MFMF_SOSim_3freq_small/selFn -c cobayaConf_SOSims_fixedScalingRelation.yml -S 5 -o small_cosmoResults1

