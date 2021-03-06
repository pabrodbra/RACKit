#!/bin/bash
# RACKit Workflow by Pablo Rodríguez Brazzarola


# Output Directory Path
#OUTPUT="/home/pablorod/results/rac_test/"
OUTPUT=$4
# Output folders
INTERMEDIATE_FILES="${OUTPUT}intermediateFiles/"
RESULTS="${OUTPUT}results/"
REPORT="${OUTPUT}report/"
INPUT="${OUTPUT}inputs/"
mkdir -p $INTERMEDIATE_FILES
mkdir -p $RESULTS
mkdir -p $REPORT

#############
### Programs
#############
# RACKIT binary
SH_DIR=$(dirname "$0")
BIN="${SH_DIR}/bin/"
# MAKEBLASTDB Path
MAKEBLASTDB=/home/pablorod/software/ncbi-blast-2.7.1+/bin/makeblastdb
# BLASTN Path
BLASTN=/home/pablorod/software/ncbi-blast-2.7.1+/bin/blastn
# MEGAN Path
MEGAN=$MEGAN
# QNUCLPARSERBLAST Path
QNUCLPARSERBLAST="${BIN}qnuclparserblast"
# FILTERPARSEDBLAST Path
FILTERPARSEDBLAST="${BIN}filterParsedBlast"
# TAXOMAKER Path
TAXOMAKER="${BIN}taxomaker"
# MERGEFULLFASTA Path
MERGEFULLFASTA="${BIN}mergeMultiFasta"
# UNISEQDBCOVERAGE Path
UNISEQDBCOVERAGE="${BIN}uniseqDBCoverage"
# RACKIT R Path
RACKIT_R_DIR="Rscript ${SH_DIR}/src/r/"
# RACKIT.py
RACKIT_PY="python3 ${SH_DIR}/src/python/rackit.py"

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
ACCESSION_TO_TAXA=/home/pablorod/data/ncbi_taxonomy/nucl_acc2tax-Mar2018.abin
# Database fasta Path
#DB="/home/pablorod/data/metagenome/RAC_test/refsoilDB.fa"
DB=$5
DB_FORMAT="${DB}" #.format"

#############
### Outputs
#############

# Reference Database
DB="${INTERMEDIATE_FILES}REFDB.fa"
cp $5 $DB
DB_FORMAT="${DB}.format"

# Grinder Ranks
GRINDER_RANKS="${INTERMEDIATE_FILES}data_generation/grinder-ranks.txt"

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
if [ "$#" -ne 5 ]; then
        echo "***USAGE*** rackit.sh <reads> <contigs> <lca_fold> <output> [<reference.fasta::or::ref_db_path>]"
        exit 1
fi

### Make blast databases
# STATUS: WORKING (WORKING/TESTING)
# OPTIONAL: YES (YES/NO)

echo "--------------- RACKIT WORKFLOW by Pablo Rodriguez -----------------" &>> ${LOG}
echo "--- Start time: $(date)" &>> ${LOG}
echo "--- Inputs:" &>> ${LOG}
echo "------ Reads Fasta: ${READS}" &>> ${LOG}
echo "------ Contigs Fasta: ${CONTIGS}" &>> ${LOG}
echo "------ LCA Coverage Threshold: ${LCA_COVERAGE}" &>> ${LOG}
echo "------ Reference DB fasta: ${DB}" &>> ${LOG}
echo "------ Output Directory: ${OUTPUT}" &>> ${LOG}


if [[ $DB == *.fa* ]]; then
        echo "###  (1/20) Now executing: makeblastdb | Reference DB" &>> ${LOG}
        ${BIN}format ${DB} ${DB_FORMAT} &>> ${LOG}
        ${MAKEBLASTDB} -in ${DB_FORMAT} -dbtype nucl -title ${REFDB_NAME} \
            -max_file_sz 2GB -out ${REFDB_PATH} &>> ${LOG}
        echo "##!  (1/20) Reference DB created" &>> ${LOG}
else
        echo "###  (1/20) Skipped: makeblastdb | Reference DB" &>> ${LOG}
        REFDB_PATH=$DB
        echo "##!  (1/20) Reference DB alredy created" &>> ${LOG}
fi


# OPTIONAL: NO (YES/NO)
echo "### (2/20)  Now executing: makeblastdb | Contig DB" &>> ${LOG}
${MAKEBLASTDB} -in ${CONTIGS} -dbtype nucl -title ${CONTIGDB_NAME} \
    -max_file_sz 2GB -out ${CONTIGDB_PATH} &>> ${LOG}
echo "##! (2/20)  Contig DB created" &>> ${LOG}

### Blast reads and contigs against DBs
# STATUS: TESTING (WORKING/TESTING)
# *OPTIONAL: NO (YES/NO)
READS_REFDB_BLAST_PATH="${INTERMEDIATE_FILES}READS_VS_REFDB.blast"
CONTIGS_REFDB_BLAST_PATH="${INTERMEDIATE_FILES}CONTIGS_VS_REFDB.blast"
READS_CONTIGS_BLAST_PATH="${INTERMEDIATE_FILES}READS_VS_CONTIGS.blast"

# Reads vs REFDB
echo "### (3/20)  Now executing: blastn | Reads vs Reference DB" &>> ${LOG}
${BLASTN} -task megablast -db ${REFDB_PATH} -evalue ${BLASTN_EVALUE} -ungapped \
    -query ${READS} -out ${READS_REFDB_BLAST_PATH} &>> ${LOG}
echo "##! (3/20)  Blast Reads vs Reference DB finished" &>> ${LOG}

# Contigs vs REFDB
echo "### (4/20)  Now executing: blastn | Contigs vs Reference DB" &>> ${LOG}
${BLASTN} -task megablast -db ${REFDB_PATH} -evalue ${BLASTN_EVALUE} -ungapped \
    -query ${CONTIGS} -out ${CONTIGS_REFDB_BLAST_PATH} &>> ${LOG}
echo "##! (4/20)  Blast Contigs vs Reference DB finished" &>> ${LOG}

# Reads vs CONTIGDB
echo "### (5/20)  Now executing: blastn | Reads vs Contig DB" &>> ${LOG}
${BLASTN} -task megablast -db ${CONTIGDB_PATH} -evalue ${BLASTN_EVALUE} -ungapped \
    -query ${READS} -out ${READS_CONTIGS_BLAST_PATH} &>> ${LOG}
echo "##! (5/20)  Blast Reads vs Contig DB finished" &>> ${LOG}

### Parse blast results + Filter
# STATUS: TESTING (WORKING/TESTING)
# OPTIONAL: NO (YES/NO)
PARSED_READS_REFDB_BLAST_PATH="${INTERMEDIATE_FILES}READS_VS_REFDB.pblast"
PARSED_CONTIGS_REFDB_BLAST_PATH="${INTERMEDIATE_FILES}CONTIGS_VS_REFDB.pblast"
PARSED_READS_CONTIGS_BLAST_PATH="${INTERMEDIATE_FILES}READS_VS_CONTIGS.pblast"
FILTER_PARSED_READS_CONTIGS="${INTERMEDIATE_FILES}FILTER_READS_VS_CONTIGS.pblast"

# QNUCLPARSERBLAST
echo "### (6/20)  Now executing: qnuclparserblast | Reads vs Reference DB" &>> ${LOG}
${QNUCLPARSERBLAST} ${READS_REFDB_BLAST_PATH} ${PARSED_READS_REFDB_BLAST_PATH} &>> ${LOG}
echo "##! (6/20)  Parsed Reads vs Reference DB Blast" &>> ${LOG}

echo "### (7/20)  Now executing: qnuclparserblast | Contigs vs Reference DB" &>> ${LOG}
${QNUCLPARSERBLAST} ${CONTIGS_REFDB_BLAST_PATH} ${PARSED_CONTIGS_REFDB_BLAST_PATH} &>> ${LOG}
echo "##! (7/20) Parsed Contigs vs Reference DB Blast" &>> ${LOG}

echo "### (8/20)  Now executing: qnuclparserblast | Reads vs Contig DB" &>> ${LOG}
${QNUCLPARSERBLAST} ${READS_CONTIGS_BLAST_PATH} ${PARSED_READS_CONTIGS_BLAST_PATH} &>> ${LOG}
echo "##! (8/20)  Parsed Reads vs Contig DB Blast" &>> ${LOG}

# FIX PARSEDBLAST
FIX_READS="${INTERMEDIATE_FILES}fixed_reads.pblast"
FIX_CONTIGS="${INTERMEDIATE_FILES}fixed_contigs.pblast"
FIX_READS_CONTIGS="${INTERMEDIATE_FILES}fixed_reads_contigs.pblast"
FILTER_READS="${INTERMEDIATE_FILES}filtered_reads_contigs.pblast"
FILTER_UNIQUE_READS="${INTERMEDIATE_FILES}unique_filtered_reads.pblast"
FILTER_UNIQUE_CONTIGS="${INTERMEDIATE_FILES}unique_filtered_contigs.pblast"

echo "### (9/20)  Now executing: Fix Parsed Blast | Reads vs Reference DB" &>> ${LOG}
sed -z 's/\n>/;>/g' ${PARSED_READS_REFDB_BLAST_PATH} | sed -z 's/\n/ /g' | sed -z 's/;/\n/g' | sed -z 's/\t/;/g' > ${FIX_READS}
echo "##! (9/20)  Fixed Reads vs Reference DB Parsed Blast" &>> ${LOG}

echo "### (10/20)  Now executing: Fix Parsed Blast | Contigs vs Reference DB" &>> ${LOG}
sed -z 's/\n>/;>/g' ${PARSED_CONTIGS_REFDB_BLAST_PATH} | sed -z 's/\n/ /g' | sed -z 's/;/\n/g' | sed -z 's/\t/;/g' > ${FIX_CONTIGS}
echo "##! (10/20)  Fixed Contigs vs Reference DB Parsed Blast " &>> ${LOG}

echo "### (11/20)  Now executing: Fix Parsed Blast |  Reads vs Contig DB" &>> ${LOG}
sed -z 's/\n>/;>/g' ${PARSED_READS_CONTIGS_BLAST_PATH} | sed -z 's/\n/ /g' | sed -z 's/;/\n/g' | sed -z 's/\t/;/g' > ${FIX_READS_CONTIGS}
echo "##! (11/20)  Fixed Reads vs Contig DB Parsed Blast " &>> ${LOG}

echo "### (12/20)  Now executing: Filter Parsed Blast | All Fixed Parsed Blast Results" &>> ${LOG}
${FILTERPARSEDBLAST} ${FIX_READS_CONTIGS} ${FILTER_READS} 0 &>> ${LOG}
echo "##! (12/20)  Filtered all Fixed Parsed Blast Results" &>> ${LOG}

### MEGAN (13-14)
# STATUS: TESTING (WORKING/TESTING)
# OPTIONAL: NO (YES/NO)
READ_LCA="${INTERMEDIATE_FILES}lca_reads"
CONTIG_LCA="${INTERMEDIATE_FILES}lca_contigs"
READ_TAXON_PATH="${INTERMEDIATE_FILES}taxon_path_reads"
CONTIG_TAXON_PATH="${INTERMEDIATE_FILES}taxon_path_contigs"

# Format to taxon path --> sed -z 's/;[a-z]__/;/g' mod.test | sed -z 's/; [0-9]*;/;/g' | sed 's/;/,"/' | sed 's/$/"/' > test.tp
echo "### (13/20)  Now executing: Megan - Blast2LCA | Reads" &>> ${LOG}
${MEGAN} -i ${READS_REFDB_BLAST_PATH} -f BlastText -m BlastN -o ${READ_LCA} -mid ${LCA_COVERAGE} -v true -a2t ${ACCESSION_TO_TAXA} &>> ${LOG}
sed -z 's/;[a-z]__/;/g' ${READ_LCA} | sed -z 's/; [0-9]*;/;/g' | sed 's/;/,"/' | sed 's/$/"/' > ${READ_TAXON_PATH}
echo "##! (13/20)  MEGAN Reads LCA calculated" &>> ${LOG}

echo "### (14/20)  Now executing: Megan - Blast2LCA | Contigs" &>> ${LOG}
${MEGAN} -i ${CONTIGS_REFDB_BLAST_PATH} -f BlastText -m BlastN -o ${CONTIG_LCA} -mid ${LCA_COVERAGE} -v true -a2t ${ACCESSION_TO_TAXA} &>> ${LOG}
sed -z 's/;[a-z]__/;/g' ${CONTIG_LCA} | sed -z 's/; [0-9]*;/;/g' | sed 's/;/,"/' | sed 's/$/"/' > ${CONTIG_TAXON_PATH}
echo "##! (14/20)  MEGAN Contigs LCA calculated" &>> ${LOG}

### REVCO (15)
# STATUS: TESTING (WORKING/TESTING)
# OPTIONAL: NO (YES/NO)
echo "### (15/20) Now executing: RACKIT Python ToolKit" &>> ${LOG}
${RACKIT_PY} ${FIX_READS_CONTIGS} ${READ_TAXON_PATH} ${CONTIG_TAXON_PATH} 10 ${RESULTS} ${FIX_READS} ${FIX_CONTIGS} ${GRINDER_RANKS} &>> ${LOG}
echo "##! (15/20)  RACKIT Python ToolKit finished successfully" &>> ${LOG}

### DB Coverage (16)
# TAXOMAKER
# STATUS: TESTING (WORKING/TESTING)
# OPTIONAL: YES (YES/NO)
UNI_DB="${INTERMEDIATE_FILES}UNI_REFDB.fa"
TAXO_FILE="${DB}.taxo"
echo "### (16/20) Now executing: Taxomaker" &>> ${LOG}
${TAXOMAKER} ${DB} 0 &>> ${LOG}
echo "##! (16/20)  Taxomaker finished successfully" &>> ${LOG}

# MERGEMULTIFAST (17)
echo "### (17/20) Now executing: MergeMultiFasta" &>> ${LOG}
${MERGEFULLFASTA} ${DB} ${UNI_DB} &>> ${LOG}
echo "##! (17/20)  MergeMultiFasta finished successfully" &>> ${LOG}

# UNISEQDBCOVERAGE (18)
COVERAGE_OUTPUT="${RESULTS}coverage"

N_READS=$(grep -c '>' ${READS})
N_CONTIGS=$(grep -c '>' ${CONTIGS})

echo "### (18/20) Now executing: UniseqDBCoverage" &>> ${LOG}
${UNISEQDBCOVERAGE} ${FIX_READS} ${FIX_CONTIGS} ${TAXO_FILE} ${UNI_DB} ${COVERAGE_OUTPUT} ${N_READS} ${N_CONTIGS} &>> ${LOG}
echo "##! (18/20)  UniseqDBCoverage finished successfully" &>> ${LOG}

### Create original distro
ORIGINAL_DISTRO="${RESULTS}refdb_distro.csv"

echo "### (19/20) Now executing: Reference DB Species Distribution" &>> ${LOG}
echo "Specie,Percentage" > $ORIGINAL_DISTRO
awk 'NR>1' intermediateFiles/data_generation/grinder-ranks.txt | while read line; do
    stringarray=($line)
    temp_accesion=${stringarray[1]}
    temp_percentage=${stringarray[2]}
    SPECIE_NAME=$(grep $temp_accesion ${DB} | cut -d ' ' -f2,3)
    echo "${SPECIE_NAME},${temp_percentage}" >> $ORIGINAL_DISTRO
done
echo "##! (19/20)  Reference DB Species Distribution finished successfully" &>> ${LOG}

### R Results (20)
# STATUS: TESTING (WORKING/TESTING)
# OPTIONAL: NO (YES/NO)
echo "### (20/20) Now executing: RACC Results R Script" &>> ${LOG}
"${RACKIT_R_DIR}RACC_script.R" $RESULTS $REPORT $RACKIT_R_DIR
echo "##! (20/20) RACC Results R Script finished" &>> ${LOG}

echo "---- Finish time: $(date)" &>>$LOG
