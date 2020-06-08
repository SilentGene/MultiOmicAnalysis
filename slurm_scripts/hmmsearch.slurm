#!/bin/bash

#SBATCH --job-name=hmmsearch
#SBATCH -p cloud,physical

#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1

#SBATCH -t 1:00:00
#=SBATCH --mail-user=your@email
#=SBATCH --mail-type=ALL

module load HMMER/3.2.1-spartan_intel-2017.u2

query=$1
bn=`basename $query`
output_dir=$2
cpus=2
db=/home/heyu/db/HgcA.hmm

hmmsearch --tblout $output_dir/$bn.hit.tab -E 1e-5 --cpu $cpus --notextw $db $query > /dev/null

esl-sfetch --index $query
grep -v "^#" $output_dir/$bn.hit.tab | awk '{print $1}' | esl-sfetch -f $query - > $query.hit.faa
