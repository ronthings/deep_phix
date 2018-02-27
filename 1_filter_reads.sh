#!/usr/bin/env bash

usage() {
  NAME=$(basename $0)
  cat <<EOF
Usage:
  ${NAME}
The script assumes a FASTQ/ folder exists and is populated with FASTQ files.
It will create the TMP/ and TRIM/ folders if they do not already exist.

EOF
}

# make the required output directories
for d in "TMP" "TRIM" ; do
  if [ ! -d $d ] ; then
    mkdir $d
  fi
done

# variables to be used in main loop
reads1=(FASTQ/*R1*.fastq.gz) # collect each forward read in array, e.g. "FASTQ/A_S1_L001_R1_001.fastq.gz"
reads1=("${reads1[@]##*/}") # [@] refers to array, greedy remove */ from left, e.g. "A_S1_L001_R1_001.fastq.gz"
reads2=("${reads1[@]/_R1/_R2}") # substitute R2 for R1, e.g. "A_S1_L001_R2_001.fastq.gz"

# location for log file
LOGFILE=./filter.log

# main loop
pipeline() {

echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> START: " $0 $@

for ((i=0; i<=${#reads1[@]}-1; i++)); do # i from zero to one minus length of array
  fwdrds="${reads1[$i]}" # e.g. "A_S1_L001_R1_001.fastq.gz"
  rvsrds="${reads2[$i]}" # e.g. "A_S1_L001_R2_001.fastq.gz"
  id="${fwdrds%%_*}" # greedy remove _ from right e.g. "A"

  # cutadapt processes PE reads simultaneously: first filter against Q10, 5' (g/G) & 3' (a/A) search, discard trimmed reads
  cutadapt -q 10,10 -g GGGCTCGG -a GACGCTGC -G GCAGCGTC -A CCGAGCCC --discard-trimmed \
  -o TMP/${id}_filtered_R1.fastq.gz -p TMP/${id}_filtered_R2.fastq.gz \
  FASTQ/${fwdrds} FASTQ/${rvsrds}

  # sickle: adaptively trim filtered reads to Q30, discard reads <20 bases, write paired + singletons to TRIM/
  sickle pe -t sanger -l 20 -n -g -f TMP/${id}_filtered_R1.fastq.gz -r TMP/${id}_filtered_R2.fastq.gz \
  -o TRIM/${id}_trimmed_R1.fastq.gz -p TRIM/${id}_trimmed_R2.fastq.gz -s TRIM/${id}_singles.fastq.gz

  # cleanup filtered reads
  rm TMP/${id}_filtered_R*.fastq.gz
done

echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> DONE."
} #pipeline end

pipeline 2>&1 | tee $LOGFILE
