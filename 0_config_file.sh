#!/usr/bin/env bash
## Configuration file for deep_seq workflow
## Edit this script to point to file locations.

# determine how many CPUs on the current machine
#  let NUMCPUS=$(sysctl -n hw.ncpu) # for MacOS
let NUMCPUS=$(cat /proc/cpuinfo | grep processor | wc -l) # for Ubuntu

# make the required output directories
# MAP: ref mapped reads; SPK: spike mapped files; TMP: temporary files; TRIM: trimmed reads; UMP: unmapped reads; VCF: phiX VCFs
for d in "MAP" "SPK" "TMP" "TRIM" "UMP" "VCF"; do
  if [ ! -d $d ] ; then
    mkdir $d
  fi
done

## locations of key files
# FASTA: file containing reference genomes (incl. ref_decoder.csv)
# FASTQ: file containing FASTQ reads
BASEDIR="/media/deepdata/HostSwitch"
FASTALOC="${BASEDIR}/FASTA"
FASTQLOC="${BASEDIR}/FASTQ"

## index genome files
bwa index ${FASTALOC}/pUC18_L09136.fasta
bwa index ${FASTALOC}/pUC18_L09136_resected.fasta
bwa index ${FASTALOC}/EC_reference.fasta
bwa index ${FASTALOC}/EC_reference_resected.fasta
bwa index ${FASTALOC}/ST_reference_FIRSTPASS.fasta
bwa index ${FASTALOC}/ST_reference_FIRSTPASS_resected.fasta
