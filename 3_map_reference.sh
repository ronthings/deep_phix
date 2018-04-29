#!/usr/bin/env bash

usage() {
  NAME=$(basename $0)
  cat <<EOF
Usage:
  ${NAME}
It is important that you run "2_subtract_spike.sh" first!
This script assumes the existence of FASTA/, TMP/ and TRIM/ folders.
- TRIM/ should contain unmapped (spike-subtracted) FASTQ
- FASTA/ should contain the {host+phage} reference genome

EOF
}

# determine how many CPUs on the current machine
#  let NUMCPUS=$(sysctl -n hw.ncpu)
let NUMCPUS=$(cat /proc/cpuinfo | grep processor | wc -l)

# create the MAP directory if it doesn't exist
if [ ! -d "MAP" ] ; then
  mkdir MAP
fi

# variables to be used in main loop
reads1=(UMP/*_unmapped_R1.fastq.gz) # collect each forward read in array, e.g. "TRIM/A_unmapped_R1.fastq.gz"
reads1=("${reads1[@]##*/}") # [@] refers to array, greedy remove */ from left, e.g. "A_unmapped_R1.fastq.gz"
reads2=("${reads1[@]/_R1/_R2}") # substitute R2 for R1, e.g. "A_unmapped_R2.fastq.gz"

# location for log file
LOGFILE=./map.log

# main loop
pipeline() {

echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> START: " $0 $@

for ((i=0; i<=${#reads1[@]}-1; i++)); do
  fwdrds="${reads1[$i]}" # e.g. "A_unmapped_R1.fastq.gz"
  rvsrds="${reads2[$i]}" # e.g. "A_unmapped_R2.fastq.gz"
  id="${fwdrds%%_*}" # greedy remove _* from right e.g. "A"

  # map (unmapped) reads against (concatenated) host and phix genomes
  # -R = adding read group ID/sample to header
  bwa index FASTA/reference.fasta
  bwa mem -t ${NUMCPUS} -R '@RG\tID:Oye\tSM:'"$id" FASTA/reference.fasta UMP/${fwdrds} UMP/${rvsrds} > TMP/${id}_refmapped.sam

  # SAM>BAM, filter for mapped reads and MAPQ>=20 (1 in 100), pipe to sort by ref position (=default, don't use -n option); cleanup
  samtools view -bS -F 4 -q 20 TMP/${id}_refmapped.sam | samtools sort -@ 3 -o MAP/${id}_refmapped.bam -
  rm TMP/${id}_refmapped.sam # delete SAM

  # index
  samtools index MAP/${id}_refmapped.bam

  # select only those reads mapping to phix and index again
  samtools view -b MAP/${id}_refmapped.bam AF176034.1 -o MAP/${id}_phixmapped.bam
  samtools index MAP/${id}_phixmapped.bam
done

echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> DONE."
} #pipeline end

pipeline 2>&1 | tee $LOGFILE
