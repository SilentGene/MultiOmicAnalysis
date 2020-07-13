#!/bin/bash

#SBATCH --job-name=iqtree
#SBATCH -p cloud,physical,bigmem

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6

#SBATCH -t 03:00:00
#SBATCH --mail-user=your@email
#SBATCH --mail-type=ALL

phy_in=$1
cpus=6

module load iqtree/1.6.7-intel-2016.u3
iqtree-mpi -s $phy_in -bb 1000 -alrt 1000 -nt $cpus
