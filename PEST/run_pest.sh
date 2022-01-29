#!/bin/bash -e
#SBATCH --job-name=pest
#SBATCH --account=uoa00014
#SBATCH --time=24:0:0
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100M
#SBATCH --mail-user=g.bogle@auckland.ac.nz
#SBATCH --mail-type=END
module load intel/2017a
# the next line will run the seven simulations and return once they are all finished
# each simulation will write its own standard output file named like batch-1.out
# this is the line that will go in the PEST batch file eventually
pest drm_monolayer.pst
