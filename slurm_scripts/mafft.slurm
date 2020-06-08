#!/bin/bash

#SBATCH --job-name=mafft
#SBATCH -p cloud

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6

#SBATCH -t 00:10:00
#SBATCH --mail-user=your@email
#SBATCH --mail-type=ALL

fas_in=$1
fas_out=$2
cpus=6

module load MAFFT/7.273-intel-2016.u3-with-extensions
mafft --localpair --quiet --maxiterate 1000 --thread $cpus $fas_in > $fas_out