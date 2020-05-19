#!/bin/bash

#SBATCH --job-name=megahit
#SBATCH -p physical,bigmem

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G

#SBATCH -t 12:00:00
#SBATCH --mail-user=heyu.lin@student.unimelb.edu.au
#SBATCH --mail-type=ALL

input_dir=$1
out_dir=$2
cpu_num=8
memory=10e10

cat $input_dir/*clean_1.fq.gz > all_1.fq.gz
cat $input_dir/*clean_2.fq.gz > all_2.fq.gz

forward_reads=all_1.fq.gz
reverse_reads=all_2.fq.gz

# megahit will create the output folder automatically

module load MEGAHIT/1.1.1-intel-2016.u3
megahit -1 $forward_reads -2 $reverse_reads -t $cpu_num -o $out_dir -m $memory
