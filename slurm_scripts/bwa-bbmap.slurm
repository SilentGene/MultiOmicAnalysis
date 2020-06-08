#!/bin/bash

#SBATCH --job-name="bwa-bbmap"
#SBATCH -p cloud,physical

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=60G

#SBATCH -t 8:00:00
#==SBATCH --mail-user=your@email
#==SBATCH --mail-type=ALL


input_fa=$1
input_reads_dir=$2
out_dir=$3
cpus=6

index_name=$input_fa.idx

module load BWA/0.7.17-intel-2016.u3
module load BBMap/36.62-intel-2016.u3-Java-1.8.0_71

if [ ! -d $OUT_DIR  ]; then
	mkdir $OUT_DIR
fi

bwa index -p $index_name $input_fa

find $input_reads_dir -name "*_1.f*q.gz" |
while read reads1
	do
		suffix=${reads1##*_}
		suffix=${suffix:1}
		reads2=${reads1%_*}'_2'$suffix
		fn=${reads1##*/}
		base=${fn%_*}
		echo "STARTING "$base

		if [ ! -f $out_dir'/'$base'_bbmap.out' ]; then
		# if the result doesn't exist, then run the following steps
			##############
			# bwa
			##############
			bwa mem -t $cpus -v 2 $index_name $reads1 $reads2 > $out_dir'/'$base'_mapping.sam'
			##############
			# bbmap
			##############
			pileup.sh usejni=t in=$out_dir'/'$base'_mapping.sam' out=$out_dir'/'$base'_bbmap.out'
			if [ -f $out_dir'/'$base'_bbmap.out' ];then
			# if the final result has been generated, the sam file can be deleted
				rm $out_dir'/'$base'_mapping.sam'
			else
			# if the result still doesn't exist after the pileip.sh step, there should be something wrong that occured.
				echo '!!! Warning: something went wrong with '$base' mapping'
			fi
		else
			echo $out_dir'/'$base'_bbmap.out was found already existed, skip this step'
		fi
	done
