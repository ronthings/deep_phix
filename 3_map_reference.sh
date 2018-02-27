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

# create the MAP directory if it doesn't exist
if [ ! -d "MAP" ] ; then
  mkdir MAP
fi

# variables to be used in main loop
reads1=(TRIM/*_unmapped_R1.fastq.gz) # collect each forward read in array, e.g. "TRIM/A_unmapped_R1.fastq.gz"
reads1=("${reads1[@]##*/}") # [@] refers to array, greedy remove */ from left, e.g. "A_unmapped_R1.fastq.gz"
reads2=("${reads1[@]/_R1/_R2}") # substitute R2 for R1, e.g. "A_unmapped_R2.fastq.gz"
readsS=("${reads1[@]/_R1/_SE}") # substitute SE for R1, e.g. "A_unmapped_SE.fastq.gz"

# location for log file
LOGFILE=./map.log

# main loop
pipeline() {

echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> START: " $0 $@

for ((i=0; i<=${#reads1[@]}-1; i++)); do
  fwdrds="${reads1[$i]}" # e.g. "A_unmapped_R1.fastq.gz"
  rvsrds="${reads2[$i]}" # e.g. "A_unmapped_R2.fastq.gz"
  id="${fwdrds%%_*}" # greedy remove _* from right e.g. "A"

  # mapping (unmapped) reads against (concatenated) host and phix genomes
  # -R = adding read group ID/sample to header
  bwa mem -t 4 -R '@RG\tID:Oye\tSM:'"$id" FASTA/reference.fasta TRIM/${readsS[$i]} > TMP/${id}_refmapped_SE.sam
  bwa mem -t 4 -R '@RG\tID:Oye\tSM:'"$id" FASTA/reference.fasta TRIM/${reads1[$i]} TRIM/${reads2[$i]} > TMP/${id}_refmapped_PE.sam

  # convert SAMs to BAMs, filtering for mapped reads with MAPQ>=20 (1 in 100) and cleanup
  samtools view -bS -F 4 -q 20 -o TMP/${id}_refmapped_SE.bam TMP/${id}_refmapped_SE.sam
  samtools view -bS -F 4 -q 20 -o TMP/${id}_refmapped_PE.bam TMP/${id}_refmapped_PE.sam
  rm TMP/${id}_refmapped_*.sam # remove SAMs

  # concatenate SE and PE BAMs and then pipe to sort by ref position (=default, don't use -n option)
  samtools cat TMP/${id}_refmapped_SE.bam TMP/${id}_refmapped_PE.bam | samtools sort -@ 3 -o MAP/${id}_refmapped_merged.bam -
  rm TMP/${id}_refmapped_*.bam # remove unmerged BAMs - can use wildcard (*) because output > other directory

  # index
  samtools index MAP/${id}_refmapped_merged.bam

  # select only those reads mapping to phix and index again
  samtools view -b MAP/${id}_refmapped_merged.bam AF176034.1 -o MAP/${id}_phixmapped.bam
  samtools index MAP/${id}_phixmapped.bam
done

echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> DONE."
} #pipeline end

pipeline 2>&1 | tee $LOGFILE
