#!/bin/bash

#SBATCH --job-name=prokka
#SBATCH -p cloud,physical,bigmem

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3

#SBATCH -t 3:00:00
#SBATCH --mail-user=heyu.lin@student.unimelb.edu.au
#SBATCH --mail-type=FAIL


# module load Prokka/1.13-intel-2016.u3-BioPerl-1.6.924
# module load Java/1.8.0_152

cpus=3
input_fa=$1
outdir=$2

bn=`basename $input_fa`
mkdir $outdir

module load Prokka/1.13-intel-2016.u3-BioPerl-1.6.924

prokka --outdir $outdir/$bn.out --prefix ${bn%%.*} --rnammer --cpus $cpus --locustag ${bn%%.*} $input
