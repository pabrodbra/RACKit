#!/bin/bash
# RACKit Workflow by Pablo Rodríguez Brazzarola

#############
### Programs
#############
# RACKIT binary
SH_DIR=$(dirname "$0")
BIN="${SH_DIR}/bin/"
# MEGAHIT Path
MEGAHIT=0
# MAKEBLASTDB Path
MAKEBLASTDB=makeblastdb
# BLASTN Path
BLASTN=blastn
# MEGAN Path
MEGAN=0
# QNUCLPARSERBLAST Path
QNUCLPARSERBLAST="${BIN}qnuclparserblast"
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
READS="/home/pablorod/data/metagenome/RAC_test/ssd_reads_v2.fa"
CONTIGS="/home/pablorod/data/metagenome/RAC_test/final.contigs.fa"
# Database fasta Path
DB="/home/pablorod/data/metagenome/RAC_test/refsoilDB.fa"
# Blast Params
BLASTN_EVALUE=0.0001

#############
### Outputs
#############
# Output Directory Path
OUTPUT="/home/pablorod/results/rac_test/"

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


##########################
### --- Execution --- ####
##########################

# DEBUG
: '
rm -r $INTERMEDIATE_FILES
rm -r $RESULTS

'
### Make blast databases
# STATUS: WORKING (WORKING/TESTING)
# OPTIONAL: NO (YES/NO)
: '
mkdir $INTERMEDIATE_FILES
mkdir $RESULTS

${MAKEBLASTDB} -in ${DB} -dbtype nucl -title ${REFDB_NAME} \
    -max_file_sz 2GB -out ${REFDB_PATH} &>> ${LOG}

# OPTIONAL: YES (YES/NO)
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
### Parse blast results
# STATUS: TESTING (WORKING/TESTING)
# OPTIONAL: NO (YES/NO)
PARSED_READS_REFDB_BLAST_PATH="${INTERMEDIATE_FILES}READS_VS_REFDB.pblast"
PARSED_CONTIGS_REFDB_BLAST_PATH="${INTERMEDIATE_FILES}CONTIGS_VS_REFDB.pblast"
PARSED_READS_CONTIGS_BLAST_PATH="${INTERMEDIATE_FILES}READS_VS_CONTIGS.pblast"

${QNUCLPARSERBLAST} ${READS_REFDB_BLAST_PATH} ${PARSED_READS_REFDB_BLAST_PATH} &>> ${LOG}
${QNUCLPARSERBLAST} ${CONTIGS_REFDB_BLAST_PATH} ${PARSED_CONTIGS_REFDB_BLAST_PATH} &>> ${LOG}
${QNUCLPARSERBLAST} ${READS_CONTIGS_BLAST_PATH} ${PARSED_READS_CONTIGS_BLAST_PATH} &>> ${LOG}

### MEGAN
# STATUS: TESTING (WORKING/TESTING)
# OPTIONAL: NO (YES/NO)

### REVCO
# STATUS: TESTING (WORKING/TESTING)
# OPTIONAL: NO (YES/NO)

### DB Coverage
# STATUS: TESTING (WORKING/TESTING)
# OPTIONAL: NO (YES/NO)

### R Results
# STATUS: TESTING (WORKING/TESTING)
# OPTIONAL: NO (YES/NO)
