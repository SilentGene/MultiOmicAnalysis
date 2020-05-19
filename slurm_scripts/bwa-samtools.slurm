#!/bin/bash

#SBATCH --job-name="bwa-samtools"
#SBATCH -p cloud,physical

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

#SBATCH -t 4:00:00
#SBATCH --mail-user=heyu.lin@student.unimelb.edu.au
#SBATCH --mail-type=FAIL

####### Desc
# Mapping reads to genes/contigs and calculate the number of mapped reads

####### Usage
# $ sbatch bwa-samtools.slurm gene.fna clean_data

########### Configuration Start ###########
cpus=4

# input
input_fa=$1
fq_dir=$2

# ouptput
out_stats_dir=$3
########### Configuration End ###########

module load BWA/0.7.17-intel-2016.u3
module load SAMtools/1.9-intel-2018.u4

# create output dir
if [ ! -d $out_stats_dir  ]; then
   mkdir $out_stats_dir
fi

# index
input_name=${input_fa##*/}
bwa index -p $out_stats_dir/$input_name.idx $input_fa

# mapping
find $fq_dir -name "*_1.fq.gz" |
while read reads1; do
    reads2=${reads1:0:-7}'2.fq.gz'
    fn=${reads1##*/}
    base=${fn%%_*}
    echo "Starting processing $base"

    # get sam & -> bam
    bwa mem -t $cpus -v 2 $out_stats_dir/$input_name.idx $reads1 $reads2 | \
        samtools view -@ $cpus -bS - > $out_stats_dir/$base.bam
    # sort bam
    samtools sort -@ $cpus -o $out_stats_dir/$base.sort.bam $out_stats_dir/$base.bam

    # sort bam
    samtools index $out_stats_dir/$base.sort.bam

    # calculate mapped reads
    samtools idxstats $out_stats_dir/$base.sort.bam > $out_stats_dir/$base.stats

    # cleanning
    rm $out_stats_dir/$base.bam
    rm $out_stats_dir/$base.sort.bam
    rm $out_stats_dir/$base.sort.bam.bai
done
