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

## Notes on Reworking
OK, I have decided to drop sickle - it doesn't appear to be maintained but [cutadapt](https://github.com/marcelm/cutadapt) is. I've also considered [trim_galore](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md), but that is still 0.n, [not super maintained](https://github.com/FelixKrueger/TrimGalore) and reduces clarity of procedure.

Illumina adapter sequences are [linked here](https://support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html) (the link therein to the detailed adapter sequences document is broken - [here's a copy](http://www.ag.unr.edu/genomics/documents/Illumina_Adapter_Sequences.pdf)) and [explained diagrammatically here](https://support.illumina.com/bulletins/2016/04/adapter-trimming-why-are-adapter-sequences-trimmed-from-only-the--ends-of-reads.html). However the diagram in the second link appears to be inaccurate with respect to the oritentation of the index 1 read. I have drawn my own diagram to clarify:

![Image](../master/resources/adapters.jpg?raw=true)

It looks common to take the first 12 nucleotides from the CTGTCTCTTATA[CACATCT] sequence provided by Illumina for Nextera and Nextera XT. This corresponds to the (identical) reverse complement of the 3' end of the read 1 and read 2 primers. Trim galore chooses a very aggressive overlap of 1 base with the adapter before 3' trimming - that would trim 1/4 reads on average (those ending C). I am going with cutadapt's default of 3 nucleotides (we can rely on read placement bias filter to deal with errors caused by ends of reads later). I have increased the stringency by allowing up to a maximum of 2 mismatches (for 10 or more nts of adapter). This goes down to 1 mismatch from 5-9 nts and to 0 mismatches below 4 nts based on these calculations:

```
0.2*12=2.4
0.2*10=2.0
0.2*9=1.8
0.2*5=1.0
0.2*4=0.8
```

Links:
* [Cutadapt manual page](http://manpages.ubuntu.com/manpages/xenial/man1/cutadapt.1.html)
* [Cutadapt user guide](https://cutadapt.readthedocs.io/en/stable/guide.html)
* [ASCII coding of Q scores](http://drive5.com/usearch/manual/quality_score.html)

Plans:
Build yaml dependency file for conda - cutadapt is on conda.


Check later - for VCFfilter:
https://www.biostars.org/p/110670/


During the day:
[Video python data science stack](https://youtu.be/EBgUiuFXE3E)
[Slides for above](http://clstaudt.me/wp-content/uploads/2016/07/PythonDataScienceEcosystem-Slides-slides.pdf)
