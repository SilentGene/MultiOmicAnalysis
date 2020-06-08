#!/bin/bash

#SBATCH --job-name="GTDBTk"
#SBATCH -p physical,bigmem

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=500G

#SBATCH -t 1-00:00:00
#SBATCH --mail-user=your@email
#SBATCH --mail-type=ALL

in_dir=$1
out_dir=$2
cpu=3

if [ ! -d $out_dir  ];then
  mkdir $out_dir
else
  echo $out_dir exist
fi

. /home/heyu/miniconda2/etc/profile.d/conda.sh
conda activate gtdbtk-env
#please just use 3 cpus when you only have 500G RAM
gtdbtk classify_wf --genome_dir $in_dir --out_dir $out_dir --cpus $cpu -x fa
conda deactivate