#!/bin/bash
# RACKit Workflow by Pablo Rodríguez Brazzarola

#############
### Programs
#############
# RACKIT binary
SH_DIR=$(dirname "$0")
BIN="${SH_DIR}/bin/"
# MAKEBLASTDB Path
MAKEBLASTDB=makeblastdb
# BLASTN Path
BLASTN=blastn
# MEGAN Path
MEGAN=0
# QNUCLPARSERBLAST Path
QNUCLPARSERBLAST="${BIN}qnuclparserblast"
# FILTERPARSEDBLAST Path
FILTERPARSEDBLAST="${BIN}filterParsedBlast"
# TAXOMAKER Path
TAXOMAKER="${BIN}taxomaker"
# MERGEFULLFASTA Path
MERGEFULLFASTA="${BIN}mergMultiFasta"
# UNISEQDBCOVERAGE Path
UNISEQDBCOVERAGE="${BIN}uniseqDBCoverage"
# R Path
RPATH=0

#############
### Inputs
#############
# Reads fasta Path
#READS="/home/pablorod/data/metagenome/RAC_test/ssd_reads_v2.fa"
#CONTIGS="/home/pablorod/data/metagenome/RAC_test/final.contigs.fa"
READS=$1
CONTIGS=$2
# Database fasta Path
#DB="/home/pablorod/data/metagenome/RAC_test/refsoilDB.fa"
DB=$3
# Blast Params
BLASTN_EVALUE=0.0001

#############
### Outputs
#############
# Output Directory Path
#OUTPUT="/home/pablorod/results/rac_test/"
OUTPUT=$4
# Output folders
INTERMEDIATE_FILES="${OUTPUT}intermediateFiles/"
RESULTS="${OUTPUT}results/"

# Log
LOG="${OUTPUT}log.txt"

# BLASTDBs
REFDB_NAME="REFDB"
CONTIGDB_NAME="CONTIGDB"
REFDB_PATH="${INTERMEDIATE_FILES}${REFDB_NAME}"
CONTIGDB_PATH="${INTERMEDIATE_FILES}${CONTIGDB_NAME}"

if [ "$#" -ne 4 ]; then
	echo "****USAGE**** rackit.sh <reads> <contigs> <reference_db> <output>

##########################
### --- Execution --- ####
##########################
mkdir -p $INTERMEDIATE_FILES
mkdir -p $RESULTS

### Make blast databases
# STATUS: WORKING (WORKING/TESTING)
# OPTIONAL: YES (YES/NO)
: '
${MAKEBLASTDB} -in ${DB} -dbtype nucl -title ${REFDB_NAME} \
    -max_file_sz 2GB -out ${REFDB_PATH} &>> ${LOG}

# OPTIONAL: NO (YES/NO)
${MAKEBLASTDB} -in ${CONTIGS} -dbtype nucl -title ${CONTIGDB_NAME} \
    -max_file_sz 2GB -out ${CONTIGDB_PATH} &>> ${LOG}

'
### Blast reads and contigs against DBs
# STATUS: TESTING (WORKING/TESTING)
# *OPTIONAL: NO (YES/NO)
READS_REFDB_BLAST_PATH="${INTERMEDIATE_FILES}READS_VS_REFDB.blast"
CONTIGS_REFDB_BLAST_PATH="${INTERMEDIATE_FILES}CONTIGS_VS_REFDB.blast"
READS_CONTIGS_BLAST_PATH="${INTERMEDIATE_FILES}READS_VS_CONTIGS.blast"
: '
# Reads vs REFDB
${BLASTN} -task megablast -db ${REFDB_PATH} -evalue ${BLASTN_EVALUE} -ungapped \
    -query ${READS} -out ${READS_REFDB_BLAST_PATH} &>> ${LOG}
# Contigs vs REFDB
${BLASTN} -task megablast -db ${REFDB_PATH} -evalue ${BLASTN_EVALUE} -ungapped \
    -query ${CONTIGS} -out ${CONTIGS_REFDB_BLAST_PATH} &>> ${LOG}
# Reads vs CONTIGDB
${BLASTN} -task megablast -db ${CONTIGDB_PATH} -evalue ${BLASTN_EVALUE} -ungapped \
    -query ${READS} -out ${READS_CONTIGS_BLAST_PATH} &>> ${LOG}
'
### Parse blast results + Filter
# STATUS: TESTING (WORKING/TESTING)
# OPTIONAL: NO (YES/NO)
PARSED_READS_REFDB_BLAST_PATH="${INTERMEDIATE_FILES}READS_VS_REFDB.pblast"
PARSED_CONTIGS_REFDB_BLAST_PATH="${INTERMEDIATE_FILES}CONTIGS_VS_REFDB.pblast"
PARSED_READS_CONTIGS_BLAST_PATH="${INTERMEDIATE_FILES}READS_VS_CONTIGS.pblast"
FILTER_PARSED_READS_CONTIGS="${INTERMEDIATE_FILES}FILTER_READS_VS_CONTIGS.pblast"

# QNUCLPARSERBLAST
${QNUCLPARSERBLAST} ${READS_REFDB_BLAST_PATH} ${PARSED_READS_REFDB_BLAST_PATH} &>> ${LOG}
${QNUCLPARSERBLAST} ${CONTIGS_REFDB_BLAST_PATH} ${PARSED_CONTIGS_REFDB_BLAST_PATH} &>> ${LOG}
${QNUCLPARSERBLAST} ${READS_CONTIGS_BLAST_PATH} ${PARSED_READS_CONTIGS_BLAST_PATH} &>> ${LOG}

# FIX PARSEDBLAST
FIX_READS="{INTERMEDEIATE_FILES}fixed_reads.pblast"
FIX_CONTIGS="{INTERMEDEIATE_FILES}fixed_contigs.pblast"
FIX_READS_CONTIGS="${INTERMEDIATE_FILES}fixed_reads_contigs.pblast"
sed -z 's/\n>/;>/g' ${PARSED_READS_REFDB_BLAST_PATH} | sed -z 's/\n/ /g' | sed -z 's/;/\n/g' | sed -z 's/\t/;/g' > ${FIX_READS}
sed -z 's/\n>/;>/g' ${PARSED_CONTIGS_REFDB_BLAST_PATH} | sed -z 's/\n/ /g' | sed -z 's/;/\n/g' | sed -z 's/\t/;/g' > ${FIX_CONTIGS}
sed -z 's/\n>/;>/g' ${PARSED_READS_CONTIGS_BLAST_PATH} | sed -z 's/\n/ /g' | sed -z 's/;/\n/g' | sed -z 's/\t/;/g' > ${FIX_READS_CONTIGS}

FILTER_READS="${INTERMEDIATE_FILES}filtered_reads_contigs.pblast"
${FILTERPARSEDBLAST} ${FIX_READS_CONTIGS} ${FILTER_READS} 0 # or 1

### MEGAN
# STATUS: TESTING (WORKING/TESTING)
# OPTIONAL: NO (YES/NO)

### REVCO
# STATUS: TESTING (WORKING/TESTING)
# OPTIONAL: NO (YES/NO)

### DB Coverage
# TAXOMAKER
# STATUS: TESTING (WORKING/TESTING)
# OPTIONAL: YES (YES/NO)
UNI_DB="${INTERMEDIATE_FILES}UNI_REFDB.fa"
TAXO_FILE=""
${TAXOMAKER} ${DB} 0

# MERGEMULTIFAST
${MERGEFULLFASTA} ${DB} ${UNI_DB}

# UNISEQDBCOVERAGE
FIX_READS="{INTERMEDEIATE_FILES}fixed_reads.pblast"
FIX_CONTIGS="{INTERMEDEIATE_FILES}fixed_contigs.pblast"
COVERAGE_OUTPUT="${RESULTS}coverage.info"
${UNISEQDBCOVERAGE} ${FIX_READS} ${FIX_CONTIGS} ${TAXO_FILE} ${UNI_DB} ${COVERAGE_OUTPUT}

### R Results
# STATUS: TESTING (WORKING/TESTING)
# OPTIONAL: NO (YES/NO)
