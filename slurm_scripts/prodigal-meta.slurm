#!/bin/bash

#SBATCH --job-name=prodigal
#SBATCH -p cloud,physical

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

#SBATCH -t 6:00:00
#SBATCH --mail-user=your@email
#SBATCH --mail-type=FAIL

module load Prodigal/2.60-intel-2016.u3

input=$1
output=$2

prodigal -q -p meta \
	-i $input -f gff -o $output/${input%.*}.gff \
	-a $output/${input%.*}.faa  \
	-d $output/${input%.*}.ffn

#-q: Run quietly

