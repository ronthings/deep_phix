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

# determine how many CPUs on the current machine
let NUMCPUS=$(sysctl -n hw.ncpu)
# let NUMCPUS=$(cat /proc/cpuinfo | grep processor | wc -l)-1

# adapter sequence for Nextera (XT)
adapter="CTGTCTCTTATA"

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

  cutadapt --quality-base=33 --quality-cutoff 30,30 \
  -a ${adapter} -A ${adapter} --error-rate=0.2 --overlap=3 \
  --trim-n --pair-filter=any --minimum-length=20 --cores=$NUMCPUS \
  -o TRIM/${id}_trimmed_R1.fastq.gz -p TRIM/${id}_trimmed_R2.fastq.gz \
  FASTQ/${fwdrds} FASTQ/${rvsrds}

done

echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> DONE."
} #pipeline end

pipeline 2>&1 | tee $LOGFILE
