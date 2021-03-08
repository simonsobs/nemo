#!/bin/sh
#SBATCH --nodes=18
#SBATCH --ntasks-per-node=16
#SBATCH --mem=64000
#SBATCH --time=03:59:00
#SBATCH --output=DR5ClusterSearch.log
#SBATCH --error=DR5ClusterSearch.err

time mpiexec nemo DR5ClusterSearch.yml -M
