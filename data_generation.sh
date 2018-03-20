#!/bin/bash
# Synth Data Generation Script by Pablo Rodriguez

MULTIFASTA=$1
N_READS=$2
OUTPUT_DIRECTORY=$3

if [ "$#" -ne 3 ]; then
	echo "***USAGE*** data_generation.sh <multifasta> <num_reads> <output_dir>"
	exit 1
fi

### ---------
# Tools
MEGAHIT=${MEGAHIT}
GRINDER=grinder
# Output Directory Path
OUTPUT=${OUTPUT_DIRECTORY}
INTERMEDIATE_FILES="${OUTPUT}intermediateFiles/"
ASM="${INTERMEDIATE_FILES}data_generation/"
ASSEMBLY="${ASM}asm/"
# Output folders
INPUTS="${OUTPUT}inputs/"
# Output names
READS_FQ="${INPUTS}reads.fastq"
READS_FA="${INPUTS}reads.fasta"
CONTIGS="${INPUTS}contigs.fa"
# Log
LOG="${OUTPUT}log.txt"

### Folder checks
# Intermediate files
mkdir -p $INTERMEDIATE_FILES
# Inputs
if [ -d $INPUTS ]; then rm -Rf $INPUTS; fi
mkdir -p $INPUTS
# Assembly
if [ -d $ASM ]; then rm -Rf $ASM; fi
mkdir -p $ASM
if [ -d $ASSEMBLY ]; then rm -Rf $ASSEMBLY; fi

### Execution
echo "### Now executing: grinder | Read Generation" &>> ${LOG}
# Fasta ${GRINDER} -rf $1 -tr $2 -rd 101 uniform 10 -md poly4 3e-3 3.3e-8 -mr 80 20 -am uniform
# FastQ
${GRINDER} -rf $1 -tr $2 -rd 101 uniform 10 -md poly4 3e-3 3.3e-8 -mr 80 20 -am uniform -od ${ASM} -ql 30 10 -fq 1
mv "${ASM}grinder-reads.fastq" ${READS_FQ} &>> ${LOG}
# Fastq to Fasta
sed -n '1~4s/^@/>/p;2~4p' ${READS_FQ} > ${READS_FA} &>> ${LOG}
# FORMAT (sed -E 's/^>([0-9^\s]*)[^=]*=([^\ ]*).*/>\1-\2/' file)
sed -i -E 's/^>([0-9^\s]*)[^=]*=([^\ ]*).*/>\1-\2/' ${READS_FA} &>> ${LOG}
echo "##! (1/2) FASTQ and FASTA reads generated" &>> ${LOG}

# MEGAHIT - Assembly
echo "### Now executing: megahit | Contig assembly" &>> ${LOG}
${MEGAHIT} --kmin-1pass --k-min 27 --k-step 10 --k-max 87 -r ${READS_FA} -o ${ASSEMBLY} &>> ${LOG}#-t 12
mv "${ASSEMBLY}final.contigs.fa" $CONTIGS &>> ${LOG}
echo "##! (2/2) Contigs assembled..." &>> ${LOG}

