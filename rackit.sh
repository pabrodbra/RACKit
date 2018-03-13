#!/bin/bash
# RACKit Workflow by Pablo Rodríguez Brazzarola

### Programs
# MAKEBLASTDB Path
MAKEBLASTDB=0
# BLASTN Path
BLASTN=0
# MEGAHIT Path
MEGAHIT=0
# MEGAN Path
MEGAN=0
# QNUCLPARSERBLAST Path
QNUCLPARSERBLAST=0
# TAXOMAKER Path
TAXOMAKER=0
# MERGEFULLFASTA Path
MERGEFULLFASTA=0
# UNISEQDBCOVERAGE Path
UNISEQDBCOVERAGE=0
# R Path
RPATH=0

### Inputs
# Reads fasta Path
READS=0
# Database fasta Path
DB=0

### Outputs
# Output Directory Path
OUTPUT_PATH=0
INTERMEDIATE_FILES="{OUTPUT_PATH}intermediateFiles/"

### --- Execution --- ####