#!/usr/bin/env bash

usage() {
  NAME=$(basename $0)
  cat <<EOF
Usage:
  ${NAME}
You must define global variables first (via "0_config_and_run.sh")

EOF
}

# location for log file
LOGFILE=./1_filter.log

# adapter sequence for Nextera (XT)
adapter="CTGTCTCTTATA"

# variables to be used in main loop
reads1=(${FASTQLOC}/*R1*.fastq.gz) # collect each forward read in array, e.g. "~/FASTQ/A_S1_L001_R1_001.fastq.gz"
reads1=("${reads1[@]##*/}") # [@] refers to array, greedy remove */ from left, e.g. "A_S1_L001_R1_001.fastq.gz"
reads2=("${reads1[@]/_R1/_R2}") # substitute R2 for R1, e.g. "A_S1_L001_R2_001.fastq.gz"

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
  ${FASTQLOC}/${fwdrds} ${FASTQLOC}/${rvsrds}

done

echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> DONE."
} #pipeline end

pipeline 2>&1 | tee $LOGFILE
