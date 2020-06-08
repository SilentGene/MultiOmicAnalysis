#!/bin/bash

#SBATCH --job-name=picard
#SBATCH -p physical

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G

#SBATCH -t 1-00:00:00
#SBATCH --mail-user=your@email
#SBATCH --mail-type=ALL

module load picard/2.18.27
bam_dir=$1

for f in $bam_dir/*.sort.bam; do
	java -Xms20g -Xmx90g -jar /usr/local/easybuild/software/picard/2.18.27/picard.jar MarkDuplicates \
	INPUT=$f \
	OUTPUT=${f%.*}.markdup.bam \
	METRICS_FILE=${f%.*}.markdup.metrics \
	AS=TRUE VALIDATION_STRINGENCY=LENIENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 REMOVE_DUPLICATES=TRUE
	rm $f
done

# picard mannual: http://broadinstitute.github.io/picard/command-line-overview.html