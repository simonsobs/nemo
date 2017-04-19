#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=64000
#SBATCH --time=10:00:00

nemo advact_s16_f150_pa2_v2.par | tee advact_s16_f150_pa2_v2.out
