# Genome-centric Multi-omic Analysis Workflow

This repository contains a collection of code and scripts used in the paper "A novel and environmentally widespread *Marinimicrobia* clade methylates mercury in low-oxygen marine environments". 

The links below the sub-headings lead to the scripts needed for the corresponding steps. Most of the scripts were developed for running on the [SLURM](https://slurm.schedmd.com/) workload manager. All code and scripts were created for and tested on [Spartan](https://dashboard.hpc.unimelb.edu.au/) HPC at The University of Melbourne.

## 1. Software used in this workflow

[Perl](https://www.perl.org/)

[Python3](https://www.python.org/)

[Trimmomatic](https://github.com/timflutre/trimmomatic)

[MEGAHIT](https://github.com/voutcn/megahit)

[metaWRAP](https://github.com/bxlab/metaWRAP)

[SPAdes](https://github.com/ablab/spades)

[Prokka](https://github.com/tseemann/prokka)

[CheckM](https://ecogenomics.github.io/CheckM/)

[GTDB-Tk](https://github.com/Ecogenomics/GTDBTk)

[HMMER](http://hmmer.org/)

[CD-HIT](https://github.com/weizhongli/cdhit)

[MAFFT](https://mafft.cbrc.jp/alignment/software/)

[IQ-TREE](http://www.iqtree.org/)

[bwa](http://bio-bwa.sourceforge.net/)

[BBMap](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/)

[MicrobeCensus](https://github.com/snayfach/MicrobeCensus)



> Take the five samples from SI047 station as an example.

## 2. Raw Data Download

 ```bash
ids="SRR3724469;SRR3724456;SRR3724482;SRR3724508;SRR3724533"  # SRA accessions for the 5 metagenomic samples
urls="SI047_ftp_urls.txt"  # output

# Using ENA API to retreive data urls
for str in ${ids//;/ } ; do echo 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/'${str:0:6}'/00'${str:0-1}'/'$str'/'$str'_1.fastq.gz' >> $urls; echo 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/'${str:0:6}'/00'${str:0-1}'/'$str'/'$str'_2.fastq.gz' >> $urls; done

wget -b -q -i $urls  # download
 ```

## 3. Read trimming

 [trimmomatic-loop.slurm](slurm_scripts\trimmomatic-loop.slurm) 

```bash
input="SI047_raw_data"
output="SI047_clean_data"
sbatch trimmomatic-loop.slurm $input $output
```

## 4. Assembly

### MEGAHIT

 [megahit-coassembly.slurm](slurm_scripts\megahit-coassembly.slurm) 

```bash
input="SI047_clean_data"
output="SI047_megahit"
sbatch megahit-coassembly.slurm $input $output
```

### Remove small contigs

 [remove_small_seqs.pl](other_scripts\remove_small_seqs.pl) 

```bash
input="SI047_megahit/final.contigs.fa"
output="SI047_megahit/SI047.2k.fa"
perl remove_small_seqs.pl 2000 $input > $output
```

## 5. Binning

### Initial Binning

 [metaWRAP-binning.slurm](slurm_scripts\metaWRAP-binning.slurm) 

```bash
input_fa="SI047_megahit/SI047.2k.fa"
input_fq="SI047_clean_data"
output="SI047_binning_out"
sbatch metaWRAP-binning $input_fa $input_fq $output  # Using --metabat1 --metabat2 --maxbin2
```

### Bins Refinement

 [metaWRAP-BinRefinement.slurm](slurm_scripts\metaWRAP-BinRefinement.slurm) 

```bash
input="SI047_binning_out"
output="SI047_bins_refined"
sbatch metaWRAP-BinRefinement.slurm $input $output  # Integrating the output from 3 binners
```

### Bins Reassembly

 [metaWRAP-ReassembleBins.slurm](slurm_scripts\metaWRAP-ReassembleBins.slurm) 

```bash
input="SI047_bins_refined/metawrap_70_10_bins"
output="SI047_bins_reassembled"
sbatch metaWRAP-ReassembleBins.slurm $input $output
```

### Final Bins

```bash
mkdir SI047_final_bins
cp SI047_bins_refined/metawrap_70_10_bins/*fa SI047_final_bins
cp -rf SI047_bins_reassembled/reassembled_bins/*fa
```

## 6. Bins Annotation

 [prokka.slurm](slurm_scripts\prokka.slurm) 

```bash
input="SI047_final_bins"
output="SI047_final_bins_prokka"
for fa in $input/*.fa;
	sbatch prokka.slurm $fa $output
done
```

## 7. Bins Quality Assessment & Classification

### CheckM

 [checkm.slurm](slurm_scripts\checkm.slurm) 

```bash
input="SI047_final_bins"
output="SI047_final_bins_checkm"
sbatch checkm.slurm $input $output
```

### GTDB-Tk

 [gtdbtk.slurm](slurm_scripts\gtdbtk.slurm) 

```bash
input="SI047_final_bins"
output="SI047_final_bins_gtdbtk"
sbatch gtdbtk.slurm $input $output
```

## 8. Searching for target genes

### HMM Search

> The hmm database for HgcA proteins is needed for this step. `HgcA.hmm` is available upon request.

 [hmmsearch.slurm](slurm_scripts\hmmsearch.slurm) 

```bash
input="SI047_final_bins_prokka"
output="SI047_HgcA_search"
find $input -name *.faa |
while read faa
	do
		sbatch hmmsearch.slurm $faa $output
	done
```

###Reduce Redundancy

 [cd-hit.slurm](slurm_scripts\cd-hit.slurm) 

```bash
cat SI047_HgcA_search/*faa > SI047_hgcA_all.faa  # need some outgroups
input="SI047_hgcA_all.faa"
output="SI047_hgcA.faa"
sbatch cd-hit.slurm $input $output
```

## 9. Phylogenetic Tree

### Alignment

 [mafft.slurm](slurm_scripts\mafft.slurm) 

```bash
input="SI047_hgcA.faa"
output="SI047_hgcA.aln.faa"
sbatch mafft.slurm $input $output
```
### ML Tree

 [fasta2relaxedphylip.py](other_scripts\fasta2relaxedphylip.py) 

 [iqtree.slurm](slurm_scripts\iqtree.slurm) 

```bash
# Change FastA to Phylip format
python fasta2relaxedphylip.py -i SI047_hgcA.aln.faa -o SI047_hgcA.aln.phy

# Make Tree
input="SI047_hgcA.aln.phy"
sbatch mafft-iqtree.slurm $input
```

## 10. Gene Abundance in metagenomic and metatranscriptomic datasets

 [bwa-bbmap.slurm](slurm_scripts\bwa-bbmap.slurm) 

 [MicrobeCensus.slurm](slurm_scripts\MicrobeCensus.slurm) 

 [bwa-samtools.slurm](slurm_scripts\bwa-samtools.slurm) 

```bash
# For metagenomic dataset
input_fa="15hgcA.fna"
input_metaG_fq="SI047_clean_data"
out_dir_metaG="hgcA_metaG_abundance_SI047"
out_dir_MC="SI047_MC"
sbatch bwa-bbmap.slurm $input_fa $input_metaG_fq $out_dir
sbatch MicrobeCensus.slurm $input_metaG_fq $out_dir_MC

# For metatranscriptomic dataset
input_fa="15hgcA.fna"
input_metaT_fq="SI047_metaT_clean_data"
out_dir_metaT="hgcA_metaT_abundance_SI047"
sbatch bwa-samtools.slurm $input_fa $input_metaT_fq $out_dir
# Count reads for metatranscriptomic data
echo $(zcat SRR3719718.fastq.gz | wc -l)/4 | bc
```

