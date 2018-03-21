# RACKit

Sort out tools/folders
Setup requirements.txt and Makefile(check)

Download Grinder (Version 0.5.4)
* Download compressed file from 'https://sourceforge.net/projects/biogrinder/'
* 'tar -xvzf community_images.tar.gz'
* 'rm Grinder-0.5.4.tar.gz'
* 'perl Makefile.PL'
* 'make'
* 'make install'

Download MegaHIT
* 'git clone https://github.com/voutcn/megahit.git'
* 'cd megahit'
* 'make'

Download MeganV6 (Version 6.10.13)
* If not OpenJFX not installed: 'sudo apt-get install openjfx'
* 'wget http://ab.inf.uni-tuebingen.de/data/software/megan6/download/MEGAN_Community_unix_6_10_13.sh'
* 'chmod +x MEGAN_Community_unix_6_10_13.sh'
* './MEGAN_Community_unix_6_10_13.sh'
* 'rm MEGAN_Community_unix_6_10_13.sh'
* 'wget http://ab.inf.uni-tuebingen.de/data/software/megan6/download/nucl_acc2tax-Mar2018.abin.zip'
* 'zip -d nucl_acc2tax-Mar2018.abin.zip'

Download Blast (Version 2.7.1)
* 'wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.7.1+-x64-linux.tar.gz'
* 'tar -xzvf ncbi-blast-2.7.1+-x64-linux.tar.gz'
* 'rm ncbi-blast-2.7.1+-x64-linux.tar.gz'

Setup Enviroment
* 'export PATH=$PATH:/home/pablorod/software/ncbi-blast-2.7.1+/bin'
* 'export RACKIT=/home/pablorod/software/RACC/rackit.sh'
* 'export MEGAHIT=/home/pablorod/software/megahit/megahit'
* 'export MEGAN=/home/pablorod/software/meganv6/tools/blast2lca'


Start testing workflow with example (Export - How to do program paths properly)
Remake REVCO Statistics
Test REVCO py
Test UniseqDBCoverage
Test R

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
