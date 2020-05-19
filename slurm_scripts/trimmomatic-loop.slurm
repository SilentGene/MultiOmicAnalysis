#!/bin/bash

#SBATCH --job-name="Trimmomatic"
#SBATCH -p cloud,physical

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=50G

#SBATCH -t 2-00:00:00
#SBATCH --mail-user=heyu.lin@student.unimelb.edu.au
#SBATCH --mail-type=ALL

cpus=6
input=$1
output=$2
# the reads files in the input should look like: aaa.1.fastq.gz, aaa.2.fastq.gz

module load Trimmomatic/0.36-Java-1.8.0_71

find $input -name *1.fastq.gz |
while read reads1
	do
		reads2=${reads1:0:-10}'2.fastq.gz'
		filename=${reads1##*/}
		samplename=${filename%%_*}
		java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -phred33 -threads $cpus \
		$reads1 $reads2 \
		$output/$samplename'_clean_1.fastq.gz' $output/$samplename'_unpaired_1.fastq.gz' \
		$output/$samplename'_clean_2.fastq.gz' $output/$samplename'_unpaired_2.fastq.gz' \
		ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
	done

# $reads1 and $reads2 will be the paths of a pair of reads. eg. project/raw_data/SRR3724315_1.fastq.gz
# $filename will be the full name of the file. eg. SRR3724315_1.fastq.gz
# $samplename will be a part of the full file name. eg. SRR3724315