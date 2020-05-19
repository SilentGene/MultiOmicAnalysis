#!/bin/bash

#SBATCH --job-name="metaWRAP-reassembleBins"
#SBATCH -p bigmem,physical

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=450G

#SBATCH -t 3-00:00:00
#SBATCH --mail-user=heyu.lin@student.unimelb.edu.au
#SBATCH --mail-type=ALL

# SPAdes uses 512 Mb per thread for buffers

clean_data=$1
output=$2

cat $clean_data/*1.fq.gz > 'clean_1.fq.gz'
cat $clean_data/*2.fq.gz > 'clean_2.fq.gz'
gunzip *.gz

cpus=8
mem=450

. /home/heyu/miniconda2/etc/profile.d/conda.sh
conda activate metawrap-env
metawrap reassemble_bins -o $output -t $cpus -m $mem \
    -1 'clean_1.fq' -2 'clean_2.fq' \
    -c 90 -x 5 -l 2000 -b $clean_data
conda deactivate
