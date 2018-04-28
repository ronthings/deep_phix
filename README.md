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

## Notes on Script 2
The second script maps the trimmed reads against the spike-in reference genome.
It does this twice - the first time against the linear genome, the second time, against a resected version of the same.

## Future Plans
* Build yaml dependency file for conda - cutadapt is on conda.
* Check later - for VCFfilter:
https://www.biostars.org/p/110670/

## To check later:
* [Video python data science stack](https://youtu.be/EBgUiuFXE3E)
* [Slides for above](http://clstaudt.me/wp-content/uploads/2016/07/PythonDataScienceEcosystem-Slides-slides.pdf)
