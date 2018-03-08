#!/bin/bash

FASTQ=$1
FASTA="${FASTQ/.fastq/}.fa"

cat $FASTQ | paste - - - - | sed 's/^@/>/g'| cut -f1-2 | tr '\t' '\n' > $FASTA

echo "${FASTQ} --> ${FASTA}"
