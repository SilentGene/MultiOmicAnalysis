#!/bin/bash

#SBATCH --job-name="htseq"
#SBATCH -p cloud,physical

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G

#SBATCH -t 1-10:00:00
#SBATCH --mail-user=your@email
#SBATCH --mail-type=ALL

bam_dir=$1
gff_file=$2
out_dir=$3

# create output dir
if [ ! -d $out_dir  ]; then
   mkdir $out_dir
fi

. /home/heyu/miniconda2/etc/profile.d/conda.sh

conda activate htseq-env

for f in $bam_dir/*markdup.bam; do
    bn=`basename $f`
	htseq-count --stranded=no -r pos -t CDS -f bam $f $gff_file > $out_dir/${bn%%.*}.count
done

conda deactivate

# For stranded=no, a read is considered overlapping with a feature regardless of 
# whether it is mapped to the same or the opposite strand as the feature. 
# For stranded=yes and single-end reads, the read has to be mapped to the same 
# strand as the feature. For paired-end reads, the first read has to be on the 
# same strand and the second read on the opposite strand. For stranded=reverse, 
# these rules are reversed