[![Generic badge](https://img.shields.io/badge/DOI-10.1038%2Fs41396--020--00889--4-blue)](https://doi.org/10.1038/s41396-020-00889-4)

# Genome-centric Multi-omic Analysis Workflow

This repository contains a collection of code and scripts used in the paper **Mercury methylation by metabolically versatile and cosmopolitan marine bacteria** (DOI: [10.1038/s41396-020-00889-4](https://doi.org/10.1038/s41396-020-00889-4)) by Lin et al.. 

The links below the sub-headings lead to the scripts needed for the corresponding steps. Most of the scripts were developed for running on the [SLURM](https://slurm.schedmd.com/) workload manager. All code and scripts were created for and tested on [Spartan](https://dashboard.hpc.unimelb.edu.au/) HPC at The University of Melbourne. You may download and adapt the scripts to suit your own requirements.

## 1. Software used in this workflow

### Software that has been integrated into Spartan system

- [Perl](https://www.perl.org/)
- [Python3](https://www.python.org/)

- [Trimmomatic](https://github.com/timflutre/trimmomatic)

- [MEGAHIT](https://github.com/voutcn/megahit)

- [Prokka](https://github.com/tseemann/prokka)

- [Prodigal](https://github.com/hyattpd/Prodigal)

- [HMMER](http://hmmer.org/)

- [CD-HIT](https://github.com/weizhongli/cdhit)

- [MAFFT](https://mafft.cbrc.jp/alignment/software/)

- [IQ-TREE](http://www.iqtree.org/)

- [bwa](http://bio-bwa.sourceforge.net/)
- [BBMap](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/)

- [Pandas](https://pandas.pydata.org/)

- [Numpy](https://numpy.org/)

### Software that needs to be installed manually 

- [metaWRAP](https://github.com/bxlab/metaWRAP)
- [CheckM](https://ecogenomics.github.io/CheckM/)

- [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk)

- [MicrobeCensus](https://github.com/snayfach/MicrobeCensus)

- [Picard](https://broadinstitute.github.io/picard/)



> Take the five samples from "SI047 S3" station as an example.

## 2. Raw Data Download

 ```bash
ids="SRR3724469;SRR3724456;SRR3724482;SRR3724508;SRR3724533"  # SRA accessions for the 5 metagenomic samples
urls="SI047_ftp_urls.txt"  # output
out_dir="SI047_raw_data"  # output

# Using ENA API to retreive data urls
for str in ${ids//;/ } ; do echo 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/'${str:0:6}'/00'${str:0-1}'/'$str'/'$str'_1.fastq.gz' >> $urls; echo 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/'${str:0:6}'/00'${str:0-1}'/'$str'/'$str'_2.fastq.gz' >> $urls; done

# Download raw fastq files to $out_dir
wget -b -q -i $urls -P $out_dir
 ```

## 3. Read trimming

 [trimmomatic-loop.slurm](slurm_scripts/trimmomatic-loop.slurm) 

```bash
input="SI047_raw_data"
output="SI047_clean_data"
sbatch trimmomatic-loop.slurm $input $output
```

## 4. Assembly

### MEGAHIT

 [megahit-coassembly.slurm](slurm_scripts/megahit-coassembly.slurm) 

```bash
input="SI047_clean_data"
output="SI047_megahit"
sbatch megahit-coassembly.slurm $input $output
```

### Remove small contigs

 [remove_small_seqs.pl](other_scripts/remove_small_seqs.pl) 

```bash
input="SI047_megahit/final.contigs.fa"
output="SI047_megahit/SI047.2k.fa"
perl remove_small_seqs.pl 2000 $input > $output
```

## 5. Binning

* `metaWRAP` conda env is required

### Initial Binning

 [metaWRAP-binning.slurm](slurm_scripts/metaWRAP-binning.slurm) 

```bash
input_fa="SI047_megahit/SI047.2k.fa"
input_fq="SI047_clean_data"
output="SI047_binning_out"
sbatch metaWRAP-binning $input_fa $input_fq $output  # Using --metabat1 --metabat2 --maxbin2
```

### Bins Refinement

 [metaWRAP-BinRefinement.slurm](slurm_scripts/metaWRAP-BinRefinement.slurm) 

```bash
input="SI047_binning_out"
output="SI047_bins_refined"
sbatch metaWRAP-BinRefinement.slurm $input $output  # Integrating the output from 3 binners
```

### Bins Reassembly

 [metaWRAP-ReassembleBins.slurm](slurm_scripts/metaWRAP-ReassembleBins.slurm) 

```bash
# Pick up some bins of interest to process reassembly
mkdir SI047_bins_interest
cp SI047_bins_refined/metawrap_70_10_bins/bin5.fa SI047_bins_interest
cp SI047_bins_refined/metawrap_70_10_bins/bin12.fa SI047_bins_interest

# Running metaWRAP
input_reads_dir="SI047_bins_refined/metawrap_70_10_bins"
input_bins_dir="SI047_bins_interest"
output="SI047_bins_reassembled"
sbatch metaWRAP-ReassembleBins.slurm $input_reads_dir $input_bins_dir $output
```

### Final Bins

```bash
mkdir SI047_final_bins
cp SI047_bins_refined/metawrap_70_10_bins/*fa SI047_final_bins
cp -rf SI047_bins_reassembled/reassembled_bins/*fa SI047_final_bins
```

## 6. Bins Annotation

 [prokka.slurm](slurm_scripts/prokka.slurm) 

```bash
input="SI047_final_bins"
output="SI047_final_bins_prokka"
for fa in $input/*.fa;
	sbatch prokka.slurm $fa $output
done
```

## 7. Bins Quality Assessment & Classification

### CheckM

 [checkm.slurm](slurm_scripts/checkm.slurm) 

```bash
input="SI047_final_bins"
output="SI047_final_bins_checkm"
sbatch checkm.slurm $input $output
```

### GTDB-Tk

 [gtdbtk.slurm](slurm_scripts/gtdbtk.slurm) 

```bash
input="SI047_final_bins"
output="SI047_final_bins_gtdbtk"
sbatch gtdbtk.slurm $input $output
```

## 8. Searching for target genes

### HMM Search

> The hmm database for HgcA proteins is needed for this step. `HgcA.hmm` is available upon request.

 [hmmsearch.slurm](slurm_scripts/hmmsearch.slurm) 

```bash
input="SI047_final_bins_prokka"
output="SI047_HgcA_search"
find $input -name *.faa |
while read faa
	do
		sbatch hmmsearch.slurm $faa $output
	done
```

### Reduce Redundancy

 [cd-hit.slurm](slurm_scripts/cd-hit.slurm) 

```bash
cat SI047_HgcA_search/*faa > SI047_hgcA_all.faa  # need some outgroups
input="SI047_hgcA_all.faa"
output="SI047_hgcA.faa"
sbatch cd-hit.slurm $input $output
```

## 9. Phylogenetic Tree

### Alignment

 [mafft.slurm](slurm_scripts/mafft.slurm) 

```bash
input="SI047_hgcA.faa"
output="SI047_hgcA.aln.faa"
sbatch mafft.slurm $input $output
```
### Maximum Likelihood (ML) Tree

 [fasta2relaxedphylip.py](other_scripts/fasta2relaxedphylip.py) 

 [iqtree.slurm](slurm_scripts/iqtree.slurm) 

```bash
# Change Fasta to Phylip format
python fasta2relaxedphylip.py -i SI047_hgcA.aln.faa -o SI047_hgcA.aln.phy

# Make Tree
input="SI047_hgcA.aln.phy"
sbatch mafft-iqtree.slurm $input
```

## 10. Gene Abundance in metagenomic datasets

 [bwa-bbmap.slurm](slurm_scripts/bwa-bbmap.slurm) 

 [MicrobeCensus.slurm](slurm_scripts/MicrobeCensus.slurm) 

```bash
input_fa="15hgcA.fna"
input_metaG_fq="SI047_clean_data"
out_dir_metaG="hgcA_metaG_abundance_SI047"
out_dir_MC="SI047_MC"
sbatch bwa-bbmap.slurm $input_fa $input_metaG_fq $out_dir
sbatch MicrobeCensus.slurm $input_metaG_fq $out_dir_MC
```
## 11. PRKM calculation for metatranscriptomic datasets

### Preparation

* `htseq` conda env is required

[prodigal.slurm](slurm_scripts/prodigal.slurm)

[bwa-samtools.slurm](slurm_scripts/bwa-samtools.slurm)

[picard.slurm](slurm_scripts/picard.slurm)

[htseq.slurm](slurm_scripts/htseq.slurm)

```bash
# ORF prediction
input_fa="SI047_megahit/SI047.2k.fa"
out_dir_prodigal="SI047_prodigal"
sbatch prodigal.slurm $input_fa $out_dir_prodigal

# mapping reads with bwa
input_fa="SI047_megahit/SI047.2k.fa"
input_metaT_fq="SI047_metaT_clean_data" # derived from a similar procedure to metagenomic clean data
out_dir="SI047_metaT_mapping"
sbatch bwa-samtools.slurm $input_fa $input_metaT_fq $out_dir  # mapping and sort
sbatch picard.slurm $out_dir  # removing duplicates

# counting mapped reads per gene
input_bams_dir="SI047_metaT_mapping"
input_gff="SI047_prodigal/SI047.gff"
output_count="SI047_reads_count"
sbatch htseq.slurm $input_bams_dir $input_gff $output_count

# calculating gene lengths
input_gff="SI047_prodigal/SI047.gff"
cut -f4,5,9 $input_gff | sed 's/gene_id //g' | gawk '{print $3,$2-$1+1}' | tr ' ' '\t' > ${input_gff%.*}.gl.txt
```

### RPKM table

RPKM: Reads per kilo base per million mapped reads

> **Formula**
>
> RPKM =   C / ( (L/1000) * (N/1,000,000) ) = (10^9 * C)/(N * L)
>
> - C - number of reads mapped to a gene sequence
> - L - gene length in base-pairs for a gene
> - N - total number of mapped reads of a sample

[RPKM_cal.py](other_scripts/RPKM_cal.py)

[Pandas](https://pandas.pydata.org/)

[Numpy](https://numpy.org/)

```bash
module load pandas/0.23.4-intel-2016.u3-Python-3.5.2  # loading Python3 and Pandas module
module load numpy/1.12.1-intel-2017.u2-Python-3.5.2 # loading Numpy module

count_dir="SI047_reads_count"
gene_len="SI047_prodigal/SI047.gl.txt"
RPKM_out="SI047.rpkm.tsv"
python RPKM_cal.py -c $count_dir -l $gene_len -o $RPKM_out
```



## Copyright

Heyu Lin [heyu.lin@student.unimelb.edu.au](mailto:heyu.lin@student.unimelb.edu.au)

School of Earth Sciences, The University of Melbourne

Please cite the article if the scripts are helpful in your research.

Lin, H., Ascher, D.B., Myung, Y. *et al.* Mercury methylation by metabolically versatile and cosmopolitan marine bacteria. *ISME J* (2021). [https://doi.org/10.1038/s41396-020-00889-4](https://doi.org/10.1038/s41396-020-00889-4)

