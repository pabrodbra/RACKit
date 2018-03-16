# RACKit

Sort out tools/folders
Setup requirements.txt and Makefile(check)

Download MetaSim
* 
* 'wget http://ab.inf.uni-tuebingen.de/data/software/metasim/download/V0_9_5/MetaSim_unix_0_9_5.sh'
* 'chmod +x MetaSim_unix_0_9_5.sh'
* './MetaSim_unix_0_9_5.sh'
* 'rm MetaSim_unix_0_9_5.sh'

OR

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