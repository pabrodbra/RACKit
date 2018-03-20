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
# Blast Params
BLASTN_EVALUE=0.0001
# Megan Params
LCA_COVERAGE=$3
# Database fasta Path
#DB="/home/pablorod/data/metagenome/RAC_test/refsoilDB.fa"
DB=$5
DB_FORMAT="${5}.format"

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

##########################
### --- Execution --- ####
##########################
mkdir -p $INTERMEDIATE_FILES
mkdir -p $RESULTS

if [ "$#" -ne 5 ]; then
        echo "***USAGE*** rackit.sh <reads> <contigs> <lca_fold> <output> [<reference.fasta::or::ref_db_path>]"
        exit 1
fi

### Make blast databases
# STATUS: WORKING (WORKING/TESTING)
# OPTIONAL: YES (YES/NO)

if [[ $DB == *.fa* ]]; then
        echo "### Now executing: makeblastdb | Reference DB" &>> ${LOG}
        ${BIN}format ${DB} ${DB_FORMAT} &>> ${LOG}
        ${MAKEBLASTDB} -in ${DB_FORMAT} -dbtype nucl -title ${REFDB_NAME} \
            -max_file_sz 2GB -out ${REFDB_PATH} &>> ${LOG}
        echo "##! (1/20) Reference DB created" &>> ${LOG}
else
        echo "### Skipped: makeblastdb | Reference DB" &>> ${LOG}
        REFDB_PATH=$DB
        echo "##! (1/20) Reference DB alredy created" &>> ${LOG}
fi


# OPTIONAL: NO (YES/NO)
echo "### Now executing: makeblastdb | Contig DB" &>> ${LOG}
${MAKEBLASTDB} -in ${CONTIGS} -dbtype nucl -title ${CONTIGDB_NAME} \
    -max_file_sz 2GB -out ${CONTIGDB_PATH} &>> ${LOG}
echo "##! (2/20) Contig DB created" &>> ${LOG}

### Blast reads and contigs against DBs
# STATUS: TESTING (WORKING/TESTING)
# *OPTIONAL: NO (YES/NO)
READS_REFDB_BLAST_PATH="${INTERMEDIATE_FILES}READS_VS_REFDB.blast"
CONTIGS_REFDB_BLAST_PATH="${INTERMEDIATE_FILES}CONTIGS_VS_REFDB.blast"
READS_CONTIGS_BLAST_PATH="${INTERMEDIATE_FILES}READS_VS_CONTIGS.blast"

echo "### Now executing: blastn | Reads vs Reference DB" &>> ${LOG}
# Reads vs REFDB
${BLASTN} -task megablast -db ${REFDB_PATH} -evalue ${BLASTN_EVALUE} -ungapped \
    -query ${READS} -out ${READS_REFDB_BLAST_PATH} &>> ${LOG}
echo "##! (3/20) Blast Reads vs Reference DB finished" &>> ${LOG}

# Contigs vs REFDB
echo "### Now executing: blastn | Contigs vs Reference DB" &>> ${LOG}
${BLASTN} -task megablast -db ${REFDB_PATH} -evalue ${BLASTN_EVALUE} -ungapped \
    -query ${CONTIGS} -out ${CONTIGS_REFDB_BLAST_PATH} &>> ${LOG}
echo "##! (4/20) Blast Contigs vs Reference DB finished" &>> ${LOG}

# Reads vs CONTIGDB
echo "### Now executing: blastn | Reads vs Contig DB" &>> ${LOG}
${BLASTN} -task megablast -db ${CONTIGDB_PATH} -evalue ${BLASTN_EVALUE} -ungapped \
    -query ${READS} -out ${READS_CONTIGS_BLAST_PATH} &>> ${LOG}
echo "##! (5/20) Blast Reads vs Contig DB finished" &>> ${LOG}

### Parse blast results + Filter
# STATUS: TESTING (WORKING/TESTING)
# OPTIONAL: NO (YES/NO)
PARSED_READS_REFDB_BLAST_PATH="${INTERMEDIATE_FILES}READS_VS_REFDB.pblast"
PARSED_CONTIGS_REFDB_BLAST_PATH="${INTERMEDIATE_FILES}CONTIGS_VS_REFDB.pblast"
PARSED_READS_CONTIGS_BLAST_PATH="${INTERMEDIATE_FILES}READS_VS_CONTIGS.pblast"
FILTER_PARSED_READS_CONTIGS="${INTERMEDIATE_FILES}FILTER_READS_VS_CONTIGS.pblast"

# QNUCLPARSERBLAST
echo "### Now executing: qnuclparserblast | Reads vs Reference DB" &>> ${LOG}
${QNUCLPARSERBLAST} ${READS_REFDB_BLAST_PATH} ${PARSED_READS_REFDB_BLAST_PATH} &>> ${LOG}
echo "##! (6/20) Parsed Reads vs Reference DB Blast" &>> ${LOG}

echo "### Now executing: qnuclparserblast | Contigs vs Reference DB" &>> ${LOG}
${QNUCLPARSERBLAST} ${CONTIGS_REFDB_BLAST_PATH} ${PARSED_CONTIGS_REFDB_BLAST_PATH} &>> ${LOG}
echo "##! (7/20) Parsed Contigs vs Reference DB Blast" &>> ${LOG}

echo "### Now executing: qnuclparserblast | Reads vs Contig DB" &>> ${LOG}
${QNUCLPARSERBLAST} ${READS_CONTIGS_BLAST_PATH} ${PARSED_READS_CONTIGS_BLAST_PATH} &>> ${LOG}
echo "##! (8/20) Parsed Reads vs Contig DB Blast" &>> ${LOG}

# FIX PARSEDBLAST
FIX_READS="${INTERMEDIATE_FILES}fixed_reads.pblast"
FIX_CONTIGS="${INTERMEDIATE_FILES}fixed_contigs.pblast"
FIX_READS_CONTIGS="${INTERMEDIATE_FILES}fixed_reads_contigs.pblast"
FILTER_READS="${INTERMEDIATE_FILES}filtered_reads_contigs.pblast"

echo "### Now executing: Fix Parsed Blast | Reads vs Reference DB" &>> ${LOG}
sed -z 's/\n>/;>/g' ${PARSED_READS_REFDB_BLAST_PATH} | sed -z 's/\n/ /g' | sed -z 's/;/\n/g' | sed -z 's/\t/;/g' > ${FIX_READS}
echo "##! (9/20)  Fixed Reads vs Reference DB Parsed Blast" &>> ${LOG}

echo "### Now executing: Fix Parsed Blast | Contigs vs Reference DB" &>> ${LOG}
sed -z 's/\n>/;>/g' ${PARSED_CONTIGS_REFDB_BLAST_PATH} | sed -z 's/\n/ /g' | sed -z 's/;/\n/g' | sed -z 's/\t/;/g' > ${FIX_CONTIGS}
echo "##! (10/20)  Fixed Contigs vs Reference DB Parsed Blast " &>> ${LOG}

echo "### Now executing: Fix Parsed Blast |  Reads vs Contig DB" &>> ${LOG}
sed -z 's/\n>/;>/g' ${PARSED_READS_CONTIGS_BLAST_PATH} | sed -z 's/\n/ /g' | sed -z 's/;/\n/g' | sed -z 's/\t/;/g' > ${FIX_READS_CONTIGS}
echo "##! (11/20)  Fixed Reads vs Contig DB Parsed Blast " &>> ${LOG}

echo "### Now executing: Filter Parsed Blast | Reads vs Contigs" &>> ${LOG}
${FILTERPARSEDBLAST} ${FIX_READS_CONTIGS} ${FILTER_READS} 0 &>> ${LOG}
echo "##! (12/20)  Filtered Reads vs Contig DB Parsed Blast " &>> ${LOG}

### MEGAN (13-14)
# STATUS: TESTING (WORKING/TESTING)
# OPTIONAL: NO (YES/NO)

### REVCO (15)
# STATUS: TESTING (WORKING/TESTING)
# OPTIONAL: NO (YES/NO)

### DB Coverage (16)
# TAXOMAKER
# STATUS: TESTING (WORKING/TESTING)
# OPTIONAL: YES (YES/NO)
: '
UNI_DB="${INTERMEDIATE_FILES}UNI_REFDB.fa"
TAXO_FILE=""
${TAXOMAKER} ${DB} 0

# MERGEMULTIFAST (17)
${MERGEFULLFASTA} ${DB} ${UNI_DB}

# UNISEQDBCOVERAGE (18)
FIX_READS="{INTERMEDEIATE_FILES}fixed_reads.pblast"
FIX_CONTIGS="{INTERMEDEIATE_FILES}fixed_contigs.pblast"
COVERAGE_OUTPUT="${RESULTS}coverage.info"
${UNISEQDBCOVERAGE} ${FIX_READS} ${FIX_CONTIGS} ${TAXO_FILE} ${UNI_DB} ${COVERAGE_OUTPUT}

### Inconsistency Solver (19)
# STATUS: TESTING (WORKING/TESTING)
# OPTIONAL: NO (YES/NO)


### R Results (20)
# STATUS: TESTING (WORKING/TESTING)
# OPTIONAL: NO (YES/NO)
'
