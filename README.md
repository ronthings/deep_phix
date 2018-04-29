# Deep sequencing analysis scripts for ΦX174 analysis

This repo contains scripts used to analyse deepseq results from recent experiments in our lab.

## Summary
Script Name | Purpose
------------|--------
1_filter_reads.sh | quality and adapter trim (cutadapt)
2_subtract_spike.sh | map against spike-in, filter BAM to keep unmapped reads, convert back to FASTQs
3_map_reference.sh | map FASTQ to bacterial genome + ΦX174, filter to MAPQ>=20, convert to SAM, concatenate and index
4_call_snps.sh | freebayes call SNPs (naïve mode), filter against biases (vcffilter), tabulate (vcf2tsv from vcflib)

## Notes on setup
This is how I set up a conda environment with the required software and stored its config:

First go get Miniconda if you don't have conda on the system - here is the example for a linux machine:
```
curl -O https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod 755 Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh
sudo ln -s /home/ubuntu/miniconda3/etc/profile.d/conda.sh /etc/profile.d/conda.sh
conda update -n base conda
```

```
conda create -n phix cutadapt bwa samtools bedtools freebayes vcflib
conda env export -n phix > environment.yml
```
You can copy this exact environment using the YAML file in this repo as follows:
```
conda env create -f environment.yml
```

## Notes on Script 1
Programs to consider:
* [Sickle](https://github.com/najoshi/sickle)
* [Cutadapt](https://github.com/marcelm/cutadapt)
* [TrimGalore](https://github.com/FelixKrueger/TrimGalore)
* [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic)

I previously used sickle, but the repo doesn't appear to have recent updates. TrimGalore is a Java wrapper script (for FastQC and cutadapt) and is in version 0.n. As discussed below some of the settings are a little stringent and it seems better to work directly with cutadapt to keep a simple and transparent workflow. I have not tried trimmomatic in detail yet.

References:
* [Cutadapt manual](http://manpages.ubuntu.com/manpages/xenial/man1/cutadapt.1.html)
* [Cutadapt user guide](https://cutadapt.readthedocs.io/en/stable/guide.html)
* [TrimGalore user guide](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md)
* [Trimmomatic manual](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf)
* [FASTQ Q score encoding](https://en.wikipedia.org/wiki/FASTQ_format#Encoding)
* [FASTQ ASCII notes](http://drive5.com/usearch/manual/quality_score.html)

A note on adapters:
Illumina lists its adapter sequences [on this page](https://support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html). The link at the bottom of that page to the detailed adapter sequences document is broken but [here's a copy](http://www.ag.unr.edu/genomics/documents/Illumina_Adapter_Sequences.pdf). Illumina also offers this [visual explanation](https://support.illumina.com/bulletins/2016/04/adapter-trimming-why-are-adapter-sequences-trimmed-from-only-the--ends-of-reads.html), but this diagram appears to be inaccurate with respect to the orientation of the index 1 read. I have drawn my own to clarify:

![Image](../master/adapters.jpg?raw=true)

It looks common (e.g., in TrimGalore) to take the first 12 nucleotides from the CTGTCTCTTATA[CACATCT] sequence provided by Illumina for Nextera and Nextera XT. This corresponds to the (identical) reverse complement of the 3' end of the read 1 and read 2 primers. Therefore this sequence is assigned to a shell variable (```adapter```) in script 1.

Let's dissect the cutadapt command line:

### Quality Score Filter
```
cutadapt --quality-base=33 --quality-cutoff 30,30 \
```
This trims based on a Q score baseline across an interval from each end. We filter against Q30 (5' trim), and against Q30 (3' trim), with default ASCII encoding (declared explictly).


### Adapter Trimming
```
-a ${adapter} -A ${adapter} --error-rate=0.2 --overlap=3 \
```
We filter only against 3' adapter sequences (a/A). The forward strand is searched with ```-a``` and the reverse strand with ```-A```. The ```-A``` option may seem redundant but is required to avoid the legacy mode ([as explained here](https://cutadapt.readthedocs.io/en/stable/guide.html#filtering-paired-end-reads)).

Overlap refers to how many adapter bases need to be present in the read for them to be removed. TrimGalore chooses a very aggressive overlap of 1 base with the adapter before 3' trimming - that would almost certainly lead to false positive removals and result in trimming of 1/4 reads on average (those ending C = first base of adapter). I am going with cutadapt's default of 3 nucleotides because we can use a read placement bias filter to deal with errors caused at the ends of reads later. I have _increased_ the stringency by allow 0.2 (20%) mismatch. According to [this part of the manual](https://cutadapt.readthedocs.io/en/stable/guide.html#error-tolerance) the number of bases is rounded down to the nearest integer, so:

Multiplication | Equals | Rounded down
---------------|--------|-------------
0.2*12 | 2.4 | 2
0.2*10 | 2.0 | 2
0.2*9 | 1.8 | 1
0.2*5 | 1.0 | 1
0.2*4 | 0.8 | 0

This allows a maximum of 2 mismatches (for >=10 adapter nts). This goes down to 1 mismatch (5-9 nts) and 0 mismatches (<=4 nts). If a read ends with the bases CTG exactly, they will be removed.

### Discarding Ns and rejecting reads1
```
--trim-n --pair-filter=any --minimum-length=20 --cores=$NUMCPUS \
```
We remove Ns from reads and discard both reads in pair (compulsory) if at least one read (--pair-filter=any) is <20 nts. Also note that we use multiple cores on the machine for read processing.

### Outputs and inputs
```
-o TRIM/${id}_trimmed_R1.fastq.gz -p TRIM/${id}_trimmed_R2.fastq.gz \
FASTQ/${fwdrds} FASTQ/${rvsrds}
```
The paired end outputs are specified with ```-o``` and ```-p``` and written to the TRIM directory while forward and reverse input files are separate (not interleaved) and obtained from the FASTQ directory.

## Notes on 2_subtract_spike.sh
The second script maps against the spike-in in five steps:
 a. map trimmed reads against spike-in and store mapped as BAM
 b. also selects unmapped reads and converts back to FASTQ (discarding singletons)
 c. map trimmed reads against "resected" spike-in and store mapped reads
 d. map unmapped reads from b against "resected" spike-in, select unmapped reads again and convert -> unmapped FASTQ

## Notes on 3_map_reference.sh
The third script maps the unmapped FASTQ against reference which includes bacterium and phix
 a. it starts by identifying whether the sample was sequenced in _Salmonella_ or _Escherichia_
 b. maps against reference and stores mapped reads as indexed BAM (also subsets for phix)
 c. subsets output of b for phiX and uses FreeBayes in naive mode to create -> unfiltered VCF
 d. maps also against reference with resected phix
 e. stores only the phix subset from this and uses FreeBayes as above to create -> unfiltered VCF (resected)

## Notes on 4_call_snps.sh
This script is incomplete but contains command line for VCFfilter from vcflib, for creating filtered VCFs.

## Next steps (for Oye's pipeline)
1. re-index the resected VCFs
2. splice from them to repair the ends of the main VCFs
3. apply vcffilter to create -> filtered VCFs
4. select UNION of all positive sites in filtered VCFs
5. examine UNION list in original, unfiltered VCFs (to avoid false negatives)

## For Alex's pipeline
1. I think #1 and 2 in the existing pipeline can be discarded because they are handled by PEAR (read merge)
2. adjust step 3 of the existing pipeline to map (PEAR output) in singleton mode but keep resecting division
3. apply the same workflow exactly as Oye's pipeline (#1-5), viz. re-index and splice VCFs, apply filters, and step back to the unfiltered VCFs.

## Setup so far
I am separate VMs (B'ham for Oye, W'wick for Alex) for these pipelines which are pulling from my git repo to fetch scripts and references. The environment (meaning all software used) is defined by a conda YAML file. Data is already copied to drives. I've tested so far using Oye's VM running for two samples only.

## URL Dump
* [CLIMB](https://bryn.climb.ac.uk)
* [VCFfilter question](https://www.biostars.org/p/110670/)
* [Read Groups Readme (at GATK)](https://software.broadinstitute.org/gatk/documentation/article.php?id=6472)
* [Command line git manual](https://schacon.github.io/git/user-manual.html#sharing-development)

## To check later:
* [Video python data science stack](https://youtu.be/EBgUiuFXE3E)
* [Slides for above](http://clstaudt.me/wp-content/uploads/2016/07/PythonDataScienceEcosystem-Slides-slides.pdf)
