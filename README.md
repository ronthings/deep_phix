# ΦX174 Deep Sequencing Analysis Pipeline
This repository contains scripts written by Ben Dickins and used to analyse deep sequencing results from the following work conducted in the Dickins laboratory:

Abbreviated Title | Working Title of Manuscript | First Author | Key Issues for Analysis
------------------|-----------------------------|--------------|------------------------
Elevated Mutation Rate Study | The effects of elevated mutation rate on evolutionary dynamics in a ssDNA phage | Alex Wilcox | overlapping read pairs, mapping over the origin
Host Switch Study | Signatures of adaptation to host switching in bacteriophage ΦX174 | Oyeronke Ayansola | subtraction of spike-in from alternate samples, mapping over the origin

## Summary
Here is a brief overview of the key scripts (in the bash_scripts/ folder) that define this pipeline:

Script Name | Purpose
------------|--------
0_config_and_run.sh | set global variables and run scripts 1-4
1_filter_reads.sh | quality and adapter trim (cutadapt)
2_subtract_spike.sh | map against spike-in (two ways) and convert unmapped reads to FASTQs
3_map_reference.sh | map FASTQ to bacterial genome + ΦX174, filter to MAPQ>=20, convert to indexed BAM
4_call_snps.sh | freebayes naïvely call SNPs, merge VCFs over origin (by coverage), anti-bias filter (vcffilter) and tabulate

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
conda create -n phix cutadapt bwa samtools bedtools freebayes vcflib pysam
conda env export -n phix > environment.yml
```
You can copy this exact environment using the YAML file in this repo as follows:
```
conda env create -f environment.yml
```

## Notes on Control Script (0_config_and_run.sh)
This script sets some global variables.

It calculates NUMCPUS = the number of CPUs on the system. This is currently set for Ubuntu, but you swap comment to activate for MacOS. Other variables include the locations of key FASTA and FASTQ files. Please check these carefully before running. The FASTA files are currently contained in this repository, but the FASTQ files will be elsewhere. The FASTA files are then indexed (only needs to happen once) and scripts 1 to 4 are executed (inheriting the variables set here).

## Notes on Trimming Script (1_filter_reads.sh)
This script searches for and removes adapters using cutadapt - however the settings are changed a little and will differ from other commonly used setting in being a little less stringent. For example a single C at the end of a read will not be removed because that could deplete coverage. A strand placement bias filter can be used later to handle this. Note that perfect matches for the first 3 (CTG) or 4 (CTGT) bases of the adapter occur in 132 and 25 places in the genome, respectively. For this analysis, please see the adapter_search/ folder of this repository.

Now follows a description of the rationale for using cutadapt as opposed to other programs. The following programs were considered:
* [Sickle](https://github.com/najoshi/sickle)
* [Cutadapt](https://github.com/marcelm/cutadapt)
* [TrimGalore](https://github.com/FelixKrueger/TrimGalore)
* [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic)

I previously used sickle, but at the time of writing the repo doesn't appear to have recent updates. TrimGalore is a Java wrapper script (for FastQC and cutadapt) and is in version 0.n. As discussed below some of the settings are a little stringent and it seems better to work directly with cutadapt to keep a simple and transparent workflow. I have not tried trimmomatic in detail yet.

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

This allows a maximum of 2 mismatches (for >=10 adapter nts). This goes down to 1 mismatch (5-9 nts) and 0 mismatches (3-4 nts). If a read ends with the bases CTG exactly, they will be removed.

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

## Notes on Spike-in Subtraction Script (2_subtract_spike.sh)
The second script maps against the spike-in in four steps:
1. map trimmed reads against spike-in and store mapped as BAM
2. also selects unmapped reads and converts back to FASTQ (discarding singletons)
3. map trimmed reads against "resected" spike-in and store mapped reads
4. map unmapped reads from step 2 against "resected" spike-in, select unmapped reads again and convert -> unmapped FASTQ

The effect of these steps is to generate coverage graphs for the spike-in (outstanding) and to subtract all reads that plausibly map to the spike-in (as well as their paired reads) from samples.

## Notes on Reference Mapping Script (3_map_reference.sh)
The third script maps the unmapped FASTQ against reference which includes bacterium and ΦX174
1. it starts by identifying whether the sample was sequenced in _Salmonella_ or _Escherichia_ (using ref_decoder.csv)
2. maps against reference, filters Q20 (should exclude multi-mappers?) and stores mapped reads as indexed BAM
3. subsets output of step 2 for the ΦX174 reference and indexes
4. maps also against reference with resected ΦX
5. subsets output of step 4 for the ΦX174 reference and indexes

BAM files generated here may be inspected. Coverage of bacterial chromosome, plasmid (if present) and phage can be examined using IGV.

## Notes on SNP Calling and Filtering Script (4_call_snps.sh)
The fourth script calls SNPs (currently simple SNPs only) and uses a series of python scripts to consolidate them:
1. uses FreeBayes in naïve mode to process ΦX174-subsetted BAM -> unfiltered VCF (e.g. A_phix_1.vcf)
2. FreeBayes is used again with the resected BAM -> unfiltered VCF (e.g. A_phix_2.vcf)
3. renumber the resected VCF: calls 2_recount_resected.py on (e.g.) A_phix_2.vcf to create (e.g.) A_phix_2_reindexed.vcf
4. merge two VCFs accounting for all variants and choosing greater coverage: calls 3_fuse_vcfs.py on (e.g.) A_phix_1.vcf and A_phix_2_reindexed.vcf to create (e.g.) A_phix_merged.vcf
5. apply vcf_parser.py script to (e.g.) A_phix_merged.vcf to tabulate by genome position with details of protein changes (to give, e.g., A_phix_unfiltered_table.tsv)
6. apply bias filters using vcffilter to (e.g.) A_phix_merged.vcf
7. apply vcf_parser.py script to output of step 6 also (to give, e.g., A_phix_filtered_table.tsv)

### To Do / Notes
* I carried out an idiot check on step 3 to see that positions in the 1 and 2 VCF are the same, not shifted (may be some uniques at edge of genome - handled in step 4). More checking may be advisable.
* I am a little uncertain of the validity of the filters applied in step 6. Are they too strict? See below.

### Dissection of freebayes command line
The following is used to refer to the reference genome (either regular or resected):
```
--fasta-reference
```

Next up, we look for all variants that pass threshold (unknown number of genome copies):
```
--pooled-continuous
```

Now, we allow even a single read different from reference to support a variant (filtering comes later), although we also specify that no fewer than 1% of the reads should support an alternative allele (this results in significant speedup - will use -@ later to handle false negatives):
```
--min-alternate-fraction 0.01 --min-alternate-count 1
```

Apply some filters on base and mapping qualities:
```
--min-mapping-quality 20 --min-base-quality 30
```

Exclude indels, multi-nucleotide events (must try without this) and complex events:
```
--no-indels --no-mnps --no-complex
```

### Dissection of vcffilter command line
Using [vcffilter](https://github.com/vcflib/vcflib#vcffilter):
```
vcffilter -f "SRP > 20" -f "SAP > 20" -f "EPP > 20" -f "QUAL > 30" -f "DP > 30"
```
SRP and SAP: strand biases
EPP: placement bias
QUAL: variant quality
DP: depth of coverage

Maybe we can stick to EPP, QUAL and DP as, according to [this comment](https://github.com/ekg/freebayes/issues/5#issuecomment-13016612), FreeBayes incorporates information about biases into its calling algorithm.

## Next steps (for Host Switching Study)
1. Consider -@ (--variant-input) and -l (--only-use-input-alleles) options in FreeBayes in order to restore alleles skipped in some samples.
2. Select UNION of all positive sites in filtered VCFs
3. Examine UNION list in original, unfiltered VCFs (to avoid false negatives) -- will need to re-run 3_map_reference.sh with -@ option

## Next steps (for Elevated Mutation Rate Study)
1. Steps 1 and 2 in the existing pipeline can be discarded because they are handled by PEAR (read merge)
2. Adjust step 3 of the existing pipeline to map (PEAR output) in singleton mode but keep resecting division

## Virtual Machine Setup
I am separate VMs (B'ham for Oye, W'wick for Alex) for these pipelines which are pulling from my git repo to fetch scripts and references. The environment (meaning all software used) is defined by a conda YAML file. Data is already copied to drives. I've tested so far using Oye's VM running for two samples only.

## URL Dump
* [CLIMB](https://bryn.climb.ac.uk)
* [VCFfilter question](https://www.biostars.org/p/110670/)
* [Read Groups Readme (at GATK)](https://software.broadinstitute.org/gatk/documentation/article.php?id=6472)
* [Command line git manual](https://schacon.github.io/git/user-manual.html#sharing-development)

## To check later:
* [Video python data science stack](https://youtu.be/EBgUiuFXE3E)
* [Slides for above](http://clstaudt.me/wp-content/uploads/2016/07/PythonDataScienceEcosystem-Slides-slides.pdf)
