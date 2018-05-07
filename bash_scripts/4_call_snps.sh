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
LOGFILE=./4_call.log

# variables to be used in main loop
reads1=(MAP/*_phixmapped_1.bam) # collect each forward read in array, e.g. "MAP/A_phixmapped_1.bam"
reads1=("${reads1[@]##*/}") # [@] refers to array, greedy remove */ from left, e.g. "A_phixmapped_1.bam"
reads2=("${reads1[@]/_1/_2}") # substitute R2 for R1, e.g. "A_phixmapped_2.bam"

# main loop
pipeline() {

echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> START: " $0 $@

for ((i=0; i<=${#reads1[@]}-1; i++)); do
  fwdrds="${reads1[$i]}" # e.g. "A_phixmapped_1.bam"
  rvsrds="${reads2[$i]}" # e.g. "A_phixmapped_2.bam"
  id="${fwdrds%%_*}" # greedy remove _* from right e.g. "A"

  # choose a reference genome
  echo ${id} processing...
  #ref=$(
  #awk -F"," -v id=$id '$1 == id { print $3 }' ${FASTALOC}/ref_decoder.csv
  #)
  #echo ${ref} selected as reference
  ref=phix_AF176034

  ## PILEUP AGAINST MAIN GENOME
  # naive variant call - should try version preserving indels
  freebayes --fasta-reference ${FASTALOC}/${ref}.fasta --pooled-continuous \
  --min-alternate-fraction 0.01 --min-alternate-count 1 \
  --min-mapping-quality 20 --min-base-quality 30 \
  --no-indels --no-mnps --no-complex MAP/${id}_phixmapped_1.bam > VCF/${id}_phix_1.vcf

  ## PILEUP AGAINST RESECTED GENOME
  # naive variant call
  freebayes --fasta-reference ${FASTALOC}/${ref}_resected.fasta --pooled-continuous \
  --min-alternate-fraction 0.01 --min-alternate-count 1 \
  --min-mapping-quality 20 --min-base-quality 30 \
  --no-indels --no-mnps --no-complex MAP/${id}_phixmapped_2.bam > VCF/${id}_phix_2.vcf

  # recount the positions in the resected VCF file
  python python_scripts/2_recount_resected.py VCF/${id}_phix_2.vcf # will create ${id}_phix_2_reindexed.vcf

  # merge VCFs and select data with highest coverage
  python python_scripts/3_fuse_vcfs.py VCF/${id}_phix_1.vcf # will create ${id}_phix_merged.vcf

  # create table with genome-level data (based on phix_coord; echo statement is passing default ARF character to script)
  echo "*" | python vcf-codon-table/vcf_parser.py ${FASTALOC}/${ref}.fasta vcf-codon-table/phix_coord.txt VCF/${id}_phix_merged.vcf VCF/${id}_phix_unfiltered_table.tsv

  # filter against various biases using vcflib's vcffilter
  vcffilter -f "SRP > 20" -f "SAP > 20" -f "EPP > 20" -f "QUAL > 30" -f "DP > 30" VCF/${id}_phix_merged.vcf > VCF/${id}_phix_filtered.vcf

  # create table with genome-level data (based on phix_coord)
  echo "*" | python vcf-codon-table/vcf_parser.py ${FASTALOC}/${ref}.fasta vcf-codon-table/phix_coord.txt VCF/${id}_phix_filtered.vcf VCF/${id}_phix_filtered_table.tsv

done

echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> DONE."
} #pipeline end

pipeline 2>&1 | tee $LOGFILE
