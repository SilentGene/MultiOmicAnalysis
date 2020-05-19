#!/bin/bash

#SBATCH --job-name=iqtree
#SBATCH -p cloud,physical,bigmem

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6

#SBATCH -t 03:00:00
#SBATCH --mail-user=heyu.lin@student.unimelb.edu.au
#SBATCH --mail-type=ALL

fas_in=$1
cpus=6

module load iqtree/1.6.7-intel-2016.u3
iqtree-mpi -s $fas_in".phy" -bb 1000 -alrt 1000 -nt $cputs
