#!/bin/bash

#SBATCH --job-name="checkM"
#SBATCH -p bigmem

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=300G

#SBATCH -t 10:00:00
#SBATCH --mail-user=your@email
#SBATCH --mail-type=ALL

module load Python/2.7.13-intel-2017.u2
module load HMMER/3.1b2-intel-2017.u2
module load Prodigal/2.60-intel-2016.u3
module load pplacer/1.1

# the input folder should contain metagenomic bins which end in '.fa'
in_dir=$1
out_dir=$2
cpu=6

# create the output folder if it doesn't exist.
if [ ! -d $out_dir  ];then
  mkdir $out_dir
else
  echo $out_dir exist
fi

# note: pplacer required more than 40GB memory for 1 cpu
checkm lineage_wf -t $cpu -x fa $in_dir $out_dir
checkm qa -o 2 -t $cpu --tab_table -f $out_dir/checkm.out $out_dir/lineage.ms $out_dir
