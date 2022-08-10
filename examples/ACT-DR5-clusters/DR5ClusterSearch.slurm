#!/bin/sh
#SBATCH --nodes=15
#SBATCH --ntasks-per-node=20
#SBATCH --mem=64000
#SBATCH --time=01:59:00
#SBATCH --output=DR5ClusterSearch.log
#SBATCH --error=DR5ClusterSearch.err

source $HOME/SETUP_CONDA.sh

time mpiexec nemo DR5ClusterSearch.yml -M
