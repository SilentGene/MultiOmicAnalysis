#!/bin/bash

#SBATCH --job-name="binrefine"
#SBATCH -p physical,bigmem

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=500G

#SBATCH -t 2-10:00:00
#SBATCH --mail-user=heyu.lin@student.unimelb.edu.au
#SBATCH --mail-type=ALL

# note: each pplacer thread uses >40GB
# 4 combination ways for 3 binners, so at least 4 cpus are recommanded

binnerA=$1/metabat1_bins
binnerB=$1/metabat2_bins
binnerC=$1/maxbin2_bins
out_dir=$2

cpu=10
memory=500

. /home/heyu/miniconda2/etc/profile.d/conda.sh
conda activate metawrap-env
metawrap bin_refinement -o $out_dir -t $cpu -m $memory -A $binnerA -B $binnerB -C $binnerC -c 70 -x 10
conda deactivate