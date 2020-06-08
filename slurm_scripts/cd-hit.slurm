#!/bin/bash

#SBATCH --job-name="cd-hit"
#SBATCH -p cloud,physical

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G

#SBATCH -t 01:00:00
#SBATCH --mail-user=your@email
#SBATCH --mail-type=ALL

#cd-hit is for Clustering peptide sequences

in_file=$1
out_file=$2
cpus=2
memory=10000

module load CD-HIT/4.6.8-spartan_intel-2017.u2


cd-hit -i $in_file -o $out_file -T $cpus -M $memory -c 0.99 -n 5 -d 0 -aS 0.5 -aL 0.5

#-c 0.99, means 99% identity, is the clustering threshold
#-n 5 is the word size
#-d 0 use sequence name in fasta header till the first white space
#-M 16000, to use 16GB RAM
#-T 8, to use 8 threads
# -S 100 length difference cutoff in amino acid, default 999999,if set to 60, the length difference between the shorter sequences and the representative of the cluster can not be bigger than 60
# -aS 0.9 alignment coverage for the shorter sequence, default 0.0, if set to 0.9, the alignment must covers 90% of the sequence
# -aL 0.9 alignment coverage for the longer sequence, default 0.0, if set to 0.9, the alignment must covers 90% of the sequence