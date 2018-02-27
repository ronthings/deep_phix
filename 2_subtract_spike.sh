#!/usr/bin/env bash

usage() {
  NAME=$(basename $0)
  cat <<EOF
Usage:
  ${NAME}
It is important that you run "1_filter_reads.sh" first!
This script assumes the existence of FASTA/, TMP/ and TRIM/ folders.
- TRIM/ should contain trimmed FASTQ
- FASTA/ should contain the plasmid reference genome

EOF
}

# variables to be used in main loop
reads1=(TRIM/*_trimmed_R1.fastq.gz) # collect each forward read in array, e.g. "TRIM/A_trimmed_R1.fastq.gz"
reads1=("${reads1[@]##*/}") # [@] refers to array, greedy remove */ from left, e.g. "A_trimmed_R1.fastq.gz"
reads2=("${reads1[@]/_R1/_R2}") # substitute R2 for R1, e.g. "A_trimmed_R2.fastq.gz"

# location for log file
LOGFILE=./subspike.log

# main loop
pipeline() {

echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> START: " $0 $@

for ((i=0; i<=${#reads1[@]}-1; i++)); do
  fwdrds="${reads1[$i]}" # e.g. "A_trimmed_R1.fastq.gz"
  rvsrds="${reads2[$i]}" # e.g. "A_trimmed_R2.fastq.gz"
  id="${fwdrds%%_*}" # greedy remove _* from right e.g. "A"
  sgls="${id}_singles.fastq.gz" # e.g. "A_singles.fastq.gz"

  # mapping SE reads (no read groups), need to make {1-2686,1-250} genome
  bwa mem -t 4 FASTA/pUC18_L09136.fasta TRIM/${sgls} > TMP/${id}_greedymapped_SE.sam

  # select UNMAPPED reads (f=flag present, 4=unmapped) and cleanup
  samtools view -b -f 4 -o TMP/${id}_unmapped_SE.bam TMP/${id}_greedymapped_SE.sam
  rm TMP/${id}_greedymapped_SE.sam # remove SAM

  # mapping PE reads, can I make these map in a single-end mode??
  bwa mem -t 4 FASTA/pUC18_L09136.fasta TRIM/${fwdrds} TRIM/${rvsrds} > TMP/${id}_greedymapped_PE.sam

  # select UNMAPPED reads here too, then sort by read name (-n, because we need to emit paired end files next) and cleanup
  # explanation: | is pipe and - is reference to intermediate output)
  samtools view -O SAM -h -f 4 TMP/${id}_greedymapped_PE.sam | samtools sort -O BAM -n -o TMP/${id}_unmapped_PE.bam -
  rm TMP/${id}_greedymapped_PE.sam # remove SAM

  # convert SE and PE BAM files to FASTQ (for PE singletons will be discarded by bedtools, which is conservative), cleanup
  bedtools bamtofastq -i TMP/${id}_unmapped_SE.bam -fq TRIM/${id}_unmapped_SE.fastq
  bedtools bamtofastq -i TMP/${id}_unmapped_PE.bam -fq TRIM/${id}_unmapped_R1.fastq -fq2 TRIM/${id}_unmapped_R2.fastq
  rm TMP/${id}_unmapped_*.bam # remove BAMs
  gzip TRIM/${id}_unmapped_*.fastq # zip FASTQs for space
done

echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> DONE."
} #pipeline end

pipeline 2>&1 | tee $LOGFILE
