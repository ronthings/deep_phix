# Deep sequencing analysis scripts for ΦX174 analysis

This repo contains scripts used to analyse deepseq results from recent experiments in our lab.

## Summary

Script Name | Purpose
------------|--------
1_filter_reads.sh | trim adapters (cutadapt), adaptive Q30 trim (sickle)
2_subtract_spike.sh | map against spike-in, filter BAM to keep unmapped reads, convert back to FASTQs
3_map_reference.sh | map FASTQ to bacterial genome + ΦX174, filter to MAPQ>=20, convert to SAM, concatenate and index
4_call_snps.sh | freebayes call SNPs (naïve mode), filter against biases (vcffilter), tabulate (vcf2tsv from vcflib)
--------------------
sickle trim creates some singletons - a lot? maybe it would be simpler to discard because they get propagated through subsequent steps
step 4 not yet automated
