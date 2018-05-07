#!/usr/bin/env bash

usage() {
  NAME=$(basename $0)
  cat <<EOF
Usage:
  ${NAME}
It is important that you run "1_filter_reads.sh" first!

EOF
}

# load config script
source bash_scripts/0_config_file.sh

# location for log file
LOGFILE=./subspike.log

# variables to be used in main loop
reads1=(TRIM/*_trimmed_R1.fastq.gz) # collect each forward read in array, e.g. "TRIM/A_trimmed_R1.fastq.gz"
reads1=("${reads1[@]##*/}") # [@] refers to array, greedy remove */ from left, e.g. "A_trimmed_R1.fastq.gz"
reads2=("${reads1[@]/_R1/_R2}") # substitute R2 for R1, e.g. "A_trimmed_R2.fastq.gz"

# main loop
pipeline() {

echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> START: " $0 $@

for ((i=0; i<=${#reads1[@]}-1; i++)); do
  fwdrds="${reads1[$i]}" # e.g. "A_trimmed_R1.fastq.gz"
  rvsrds="${reads2[$i]}" # e.g. "A_trimmed_R2.fastq.gz"
  id="${fwdrds%%_*}" # greedy remove _* from right e.g. "A"

  ## 1. MAP {TRIMMED READS} TO MAIN GENOME
  # map reads to the linear spike-in genome (edge effects expected) and send to TMP/
  bwa mem -t ${NUMCPUS} ${FASTALOC}/pUC18_L09136.fasta TRIM/${fwdrds} TRIM/${rvsrds} > TMP/${id}_tmpspike_1.sam

  # select MAPPED reads (F=flag absent, 4=unmapped), send to SPK/ for later
  samtools view -O BAM -F 4 -o SPK/${id}_spikemapped_1.bam TMP/${id}_tmpspike_1.sam

  # select UNMAPPED reads (f=flag present, 4=unmapped), sort by read name and keep BAM in TMP/
  samtools view -O SAM -h -f 4 TMP/${id}_tmpspike_1.sam | samtools sort -O BAM -n -o TMP/${id}_unmapped_1.bam -

  # convert unmapped BAM files to FASTQ (for PE singletons are discarded by bedtools, which is conservative)
  bedtools bamtofastq -i TMP/${id}_unmapped_1.bam -fq TMP/${id}_unmapped_R1.fastq -fq2 TMP/${id}_unmapped_R2.fastq

  # delete unmapped BAM and temporary SAM (so we only have unmapped FASTQ in TMP/)
  rm TMP/${id}_unmapped_1.bam
  rm TMP/${id}_tmpspike_*.sam

  ## 2. MAP {TRIMMED READS} TO RESECTED GENOME (this is for later analysis of spike-in genome)
  # map trimmed reads to resected spike-in genome
  bwa mem -t ${NUMCPUS} ${FASTALOC}/pUC18_L09136_resected.fasta TRIM/${fwdrds} TRIM/${rvsrds} > TMP/${id}_tmpspike_2.sam

  # select MAPPED reads (F=flag absent, 4=unmapped), send to SPK/ for later, and delete temporary SAM
  samtools view -O BAM -F 4 -o SPK/${id}_spikemapped_2.bam TMP/${id}_tmpspike_2.sam
  rm TMP/${id}_tmpspike_*.sam

  ## 3. MAP {UNMAPPED READS in TMP/} TO RESECTED GENOME (this is for maximal accuracy in subtracting spike-in reads)
  # map unmapped reads to resected spike-in genome
  bwa mem -t ${NUMCPUS} ${FASTALOC}/pUC18_L09136_resected.fasta TMP/${id}_unmapped_R1.fastq TMP/${id}_unmapped_R2.fastq > TMP/${id}_tmpspike_3.sam

  # remove input FASTQs
  rm TMP/${id}_unmapped_*.fastq

  # select UNMAPPED reads (f=flag present, 4=unmapped), use grep -v to remove MAPPED reads from filter list, sort by read name and cleanup
  samtools view -O SAM -h -f 4 TMP/${id}_tmpspike_3.sam | samtools sort -O BAM -n -o TMP/${id}_unmapped_2.bam -
  rm TMP/${id}_tmpspike_*.sam # delete temporary SAM

  # convert unmapped BAM files to FASTQ, send to UMP/ and cleanup
  bedtools bamtofastq -i TMP/${id}_unmapped_2.bam -fq UMP/${id}_unmapped_R1.fastq -fq2 UMP/${id}_unmapped_R2.fastq
  rm TMP/${id}_unmapped_2.bam # remove BAM

  # compress FASTQ files in UMP/
  gzip -f UMP/${id}_unmapped_*.fastq # zip FASTQs for space (-f forces deletion of original)

done

echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> DONE."
} #pipeline end

pipeline 2>&1 | tee $LOGFILE


## dead darlings (grep is too slow)
# create filter list (stored in TMP) of mapped reads to remove
#grep -v "^@" SPK/${id}_spikemapped_1.sam | cut -f1 | sort -n | uniq > TMP/${id}_removelist.txt
# select UNMAPPED reads (f=flag present, 4=unmapped), use grep -v to remove MAPPED reads from filter list, sort by read name and cleanup
#samtools view -O SAM -h -f 4 TMP/${id}_tmpspike_2.sam | grep -vf TMP/${id}_removelist.txt | samtools sort -O BAM -n -o
