# RACKit

Sort out tools/folders
Setup requirements.txt and Makefile(check)

Start testing workflow with example

------

Flow:
Data_Generation:
    a. obtain samples
    b. generate reads with selected distribution
    c. generate noise
    d. format + cat all reads
    e. assemble contigs (megahit) 

Setup:
    1. generate blastdb (makeblastdb)
    2. generate contigs_blastdb (makeblastdb)
Blastn:
    3. READS vs DB
    4. CONTIGS vs DB
    5. READS vs CONTIGS_DB
Parse_Results:
    6. parse all blastn results (qnuclparserblast + format_parsed_blast.sh)
    7. generate RMA (megan taxonomy) for (READS && CONTIGS) vs DB (MEGANv6 -- figure out cmd prompt)
REVCO:
    8. REVCO execution (revco.py)
UniSeqDB:
    9. index DB (taxomaker)
    10. create DB unisequence (mergeMuliFasta)
    11. Fill coverage of DB unisequence (uniseqDBCoverage)
R Scripts:
    12. FIX THEM