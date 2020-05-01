#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32000
#SBATCH --time=02:59:00

source ~/.bashrc
export OMP_NUM_THREADS=20
export NTHREADS=20
time python createSurveyMasks_v7_post10mJy.py
