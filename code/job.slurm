#!/usr/bin/env bash
#SBATCH --time=80:00:00
#SBATCH --mem=80gb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --job-name=coastline
#SBATCH --mail-user=chunlidai1@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -o out1   
cd $SLURM_SUBMIT_DIR
#module load matlab/r2020b
module load pkgsrc  #gdalinfo
module load matlab
exec matlab < CoastTileMonoMain.m
