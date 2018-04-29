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

# determine how many CPUs on the current machine
#  let NUMCPUS=$(sysctl -n hw.ncpu)
let NUMCPUS=$(cat /proc/cpuinfo | grep processor | wc -l)

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

  # mapping PE reads to the linear spike-in genome (edge effects expected)
  bwa index FASTA/pUC18_L09136.fasta
  bwa mem -t $NUMCPUS FASTA/pUC18_L09136.fasta TRIM/${fwdrds} TRIM/${rvsrds} > TMP/${id}_greedymapped_PE.sam

  # select UNMAPPED reads (f=flag present, 4=unmapped), sort by read name (-n, because we need to emit paired end files next) and cleanup
  # explanation: | is pipe and - is reference to intermediate output)
  samtools view -O SAM -h -f 4 TMP/${id}_greedymapped_PE.sam | samtools sort -O BAM -n -o TMP/${id}_unmapped_PE.bam -
  #rm TMP/${id}_greedymapped_PE.sam # remove SAM

  # convert PE BAM files to FASTQ (for PE singletons are discarded by bedtools, which is conservative), cleanup
  bedtools bamtofastq -i TMP/${id}_unmapped_PE.bam -fq TMP/${id}_protomapped_R1.fastq -fq2 TMP/${id}_protomapped_R2.fastq
  rm TMP/${id}_unmapped_*.bam # remove BAMs

  # mapping PE reads to resected spike-in genome
  #bwa index FASTA/pUC18_L09136_resected.fasta
  #bwa mem -t $NUMCPUS FASTA/pUC18_L09136_resected.fasta TMP/${id}_protomapped_R1.fastq TMP/${id}_protomapped_R2.fastq > TMP/${id}_greedymapped_PE.sam
  #rm TMP/${id}_protomapped_*.fastq # remove intermediate FASTQs

  # select UNMAPPED reads, sort by read name and cleanup (as above)
  #samtools view -O SAM -h -f 4 TMP/${id}_greedymapped_PE.sam | samtools sort -O BAM -n -o TMP/${id}_unmapped_PE.bam -
  #rm TMP/${id}_greedymapped_PE.sam # remove SAM

  # convert PE BAM files to FASTQ and cleanup
  #bedtools bamtofastq -i TMP/${id}_unmapped_PE.bam -fq TRIM/${id}_unmapped_R1.fastq -fq2 TRIM/${id}_unmapped_R2.fastq
  #rm TMP/${id}_unmapped_*.bam # remove BAMs
  #gzip -f TRIM/${id}_unmapped_*.fastq # zip FASTQs for space (-f forces deletion of original)

done

echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> DONE."
} #pipeline end

pipeline 2>&1 | tee $LOGFILE
