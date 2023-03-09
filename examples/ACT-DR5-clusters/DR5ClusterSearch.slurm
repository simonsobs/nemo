#!/bin/sh
#SBATCH --nodes=15
#SBATCH --ntasks-per-node=20
#SBATCH --mem=64000
#SBATCH --time=04:59:00
#SBATCH --output=DR5ClusterSearch.log
#SBATCH --error=DR5ClusterSearch.err

time mpiexec nemo DR5ClusterSearch.yml -M
