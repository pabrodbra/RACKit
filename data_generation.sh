#!/bin/bash
# Synth Data Generation Script by Pablo Rodriguez

# datageneration.sh <multifasta> <num_reads>
MULTIFASTA=$1
N_READS=$2


# RACKIT binary
SH_DIR=$(dirname "$0")
BIN="${SH_DIR}/bin/"
# MEGAHIT Path
MEGAHIT=${MEGAHIT}
GRINDER=grinder
# Output Directory Path
OUTPUT="/home/pablorod/results/rac_test/"
# Output folders
RESULTS="${OUTPUT}inputs/"
READS_FQ="${RESULTS}reads.fa"
READS_FA="${RESULTS}reads.fa"
CONTIGS="${RESULTS}contigs.fa"

mkdir $RESULTS
#grinder -rf test_bacterial.fa -tr 10000 -rd 101 uniform 10 -md poly4 3e-3 3.3e-8 -mr 80 20 -am uniform -ql 30 10 -fq 1
# Fasta
${GRINDER} -rf $1 -tr $2 -rd 101 uniform 10 -md poly4 3e-3 3.3e-8 -mr 80 20 -am uniform
# FastQ
${GRINDER} -rf $1 -tr $2 -rd 101 uniform 10 -md poly4 3e-3 3.3e-8 -mr 80 20 -am uniform -ql 30 10 -fq 1
mv "${SH_DIR}/grinder-reads.fastq" ${READS_FQ}
rm "${SH_DIR}/grinder*"
# Fastq to Fasta
sed -n '1~4s/^@/>/p;2~4p' ${READS_FQ} > ${READS_FA}

# MEGAHIT - Assembly
${MEGAHIT} --kmin-1pass --k-min 27 --k-step 10 --k-max 87 -r ${READS_FA} -o CONTIGS #-t 12