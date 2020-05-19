#!/bin/bash

#SBATCH --job-name="binning"
#SBATCH -p cloud,physical

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=50G

#SBATCH -t 2-00:00:00
#SBATCH --mail-user=heyu.lin@student.unimelb.edu.au
#SBATCH --mail-type=ALL

# cpus are important for mapping. Consider to use the same number of CPUs with the number of samples.

#gunzip clean_data/*.gz
#rename .fq .fastq clean_data/*.fq

contigs=$1

# the input reads name must end with .fastq (not .fq)
reads=$2/*fastq

# the output dir will be created automatically
out_dir=$3
cpu=6
memory=50

. /home/heyu/miniconda2/etc/profile.d/conda.sh
conda activate metawrap-env
metawrap binning -a $contigs -o $out_dir -t $cpu -m 50 --metabat1 --metabat2 --maxbin2 $reads
conda deactivate