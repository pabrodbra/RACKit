#!/bin/bash
# Synth Data Generation Script by Pablo Rodriguez

# datageneration.sh <multifasta> <num_reads>
MULTIFASTA=$1
N_READS=$2

if [ "$#" -ne 2 ]; then
	echo "***USAGE*** data_generation.sh <multifasta> <num_reads>"
	exit 1
fi

### ---------
# Tools
MEGAHIT=${MEGAHIT}
GRINDER=grinder
# Output Directory Path
OUTPUT="/home/pablorod/results/rac_test/"
INTERMEDIATE_FILES="${OUTPUT}intermediateFiles/"
ASM="${INTERMEDIATE_FILES}data_generation/"
# Output folders
INPUTS="${OUTPUT}inputs/"
# Output names
READS_FQ="${INPUTS}reads.fastq"
READS_FA="${INPUTS}reads.fasta"
CONTIGS="${INPUTS}contigs.fa"

### Folder checks
# Intermediate files
mkdir -p $INTERMEDIATE_FILES
# Inputs
if [ -d $INPUTS ]; then rm -Rf $INPUTS; fi
mkdir -p $INPUTS
# Assembly
if [ -d $ASM ]; then rm -Rf $ASM; fi

### Execution
# Fasta ${GRINDER} -rf $1 -tr $2 -rd 101 uniform 10 -md poly4 3e-3 3.3e-8 -mr 80 20 -am uniform
# FastQ
${GRINDER} -rf $1 -tr $2 -rd 101 uniform 10 -md poly4 3e-3 3.3e-8 -mr 80 20 -am uniform -od ${ASM} -ql 30 10 -fq 1
mv "${ASM}grinder-reads.fastq" ${READS_FQ}
# Fastq to Fasta
sed -n '1~4s/^@/>/p;2~4p' ${READS_FQ} > ${READS_FA}

# MEGAHIT - Assembly
${MEGAHIT} --kmin-1pass --k-min 27 --k-step 10 --k-max 87 -r ${READS_FA} -o ${ASM} #-t 12
mv "${ASM}final.contigs.fa" $CONTIGS
