#!/usr/bin/env bash

usage() {
  NAME=$(basename $0)
  cat <<EOF
Usage:
  ${NAME}
It is important that you run "2_subtract_spike.sh" first!

EOF
}

# location for log file
LOGFILE=./3_map.log

# variables to be used in main loop
reads1=(UMP/*_unmapped_R1.fastq.gz) # collect each forward read in array, e.g. "UMP/A_unmapped_R1.fastq.gz"
reads1=("${reads1[@]##*/}") # [@] refers to array, greedy remove */ from left, e.g. "A_unmapped_R1.fastq.gz"
reads2=("${reads1[@]/_R1/_R2}") # substitute R2 for R1, e.g. "A_unmapped_R2.fastq.gz"

# main loop
pipeline() {

echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> START: " $0 $@

for ((i=0; i<=${#reads1[@]}-1; i++)); do
  fwdrds="${reads1[$i]}" # e.g. "A_unmapped_R1.fastq.gz"
  rvsrds="${reads2[$i]}" # e.g. "A_unmapped_R2.fastq.gz"
  id="${fwdrds%%_*}" # greedy remove _* from right e.g. "A"

  # choose a reference genome
  echo ${id} processing...
  ref=$(
  awk -F"," -v id=$id '$1 == id { print $3 }' ${FASTALOC}/ref_decoder.csv
  )
  echo ${ref} selected as reference

  ## MAP TO MAIN GENOME
  # map (unmapped) reads against (concatenated) host and phix genomes
  # -R = adding read group ID/sample to header
  bwa mem -t ${NUMCPUS} -R '@RG\tID:HSTSWTCH\tSM:'"$id" ${FASTALOC}/${ref}.fasta UMP/${fwdrds} UMP/${rvsrds} > TMP/${id}_refmapped_1.sam

  # SAM>BAM, filter for mapped reads and MAPQ>=20 (1 in 100), pipe to sort by ref position (=default, don't use -n option); cleanup
  samtools view -bS -F 4 -q 20 TMP/${id}_refmapped_1.sam | samtools sort -@ 3 -o MAP/${id}_refmapped_1.bam -
  rm TMP/${id}_refmapped_*.sam # delete SAM

  # let's index
  samtools index MAP/${id}_refmapped_1.bam

  # select only those reads mapping to phix and index again
  samtools view -b MAP/${id}_refmapped_1.bam AF176034.1 -o MAP/${id}_phixmapped_1.bam
  samtools index MAP/${id}_phixmapped_1.bam

  ## MAP TO RESECTED GENOME
  # map (unmapped) reads against (concatenated) host and RESECTED phix genomes
  bwa mem -t ${NUMCPUS} -R '@RG\tID:HSTSWTCH\tSM:'"$id" ${FASTALOC}/${ref}_resected.fasta UMP/${fwdrds} UMP/${rvsrds} > TMP/${id}_refmapped_2.sam

  # SAM>BAM with filter and sort - but push to BAM to TMP (because we won't keep this one)
  samtools view -bS -F 4 -q 20 TMP/${id}_refmapped_2.sam | samtools sort -@ 3 -o TMP/${id}_refmapped_2.bam -
  rm TMP/${id}_refmapped_*.sam # delete SAM

  # let's index
  samtools index TMP/${id}_refmapped_2.bam

  # fetch BAM from TMP and select only those reads mapping to resected phix, index again and cleanup
  samtools view -b TMP/${id}_refmapped_2.bam RESTART_2694_RESECTED_AF176034.1 -o MAP/${id}_phixmapped_2.bam
  samtools index MAP/${id}_phixmapped_2.bam
  rm TMP/${id}_refmapped_*.* # remove BAM

done

echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> DONE."
} #pipeline end

pipeline 2>&1 | tee $LOGFILE
