# RACKit

Sort out tools/folders
Setup requirements.txt and Makefile(check)

Start testing workflow with example

------

Flow:
Read_Generation:
    a. obtain samples
    b. generate reads with selected distribution
    c. generate noise
    d. format + cat all reads

Setup:
    1. generate blastdb (makeblastdb)
    2. assemble contigs (megahit)
    3. generate contigs_blastdb (makeblastdb)
Blastn:
    4. READS vs DB
    5. CONTIGS vs DB
    6. READS vs CONTIGS_DB
Parse_Results:
    7. parse all blastn results (qnuclparserblast + format_parsed_blast.sh)
    8. generate RMA (megan taxonomy) for (READS && CONTIGS) vs DB (MEGANv6 -- figure out cmd prompt)
REVCO:
    9. REVCO execution (revco.py)
UniSeqDB:
    10. index DB (taxomaker)
    11. create DB unisequence (mergeMuliFasta)
    12. Fill coverage of DB unisequence (uniseqDBCoverage)
R Scripts:
    13. FIX THEM