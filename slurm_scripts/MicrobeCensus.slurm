#!/bin/bash

#SBATCH --job-name="bwa-bbmap-mc"
#SBATCH -p cloud,physical,bigmem

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=60G

#SBATCH -t 10:00:00
#SBATCH --mail-user=your@email
#SBATCH --mail-type=ALL

input_reads_dir=$1
out_dir=$2
cpus=6

module load Biopython/1.67-intel-2016.u3-Python-3.5.2
module load numpy/1.12.1-intel-2017.u2-Python-3.5.2

mkdir $out_dir

find $input_reads_dir -name "*_1.f*q.gz" |
while read reads1
    do
        suffix=${reads1##*_}
        suffix=${suffix:1}
        reads2=${reads1%_*}'_2'$suffix
        fn=${reads1##*/}
        base=${fn%_*}

        echo "STARTING "$base
        if [ ! -f $out_dir'/'$base'_AGS' ]; then
            python ~/software/MicrobeCensus/scripts/run_microbe_census.py \
                -t $cpus $reads1,$reads2 $out_dir'/'$base'_AGS'
        else
            $out_dir'/'$base'_AGS was found already existed, skip this step'
        fi
    done