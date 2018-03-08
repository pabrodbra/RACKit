rm(list=ls())
source("C:/Users/Blinsky.Blinsk/Documents/- UMA/BitLab/Metagenome/Scripts/R/RACC_functions.R")

# Obtain Synthetic Dataset Results

### SYNTHETIC DATASET
setwd("C:/Users/Blinsky.Blinsk/Documents/- UMA/BitLab/Metagenome/Synthetic Dataset/Results")
#pablorod@papaya:~/results/REVCO/semi_synth_dataset$ grep ">" -c /home/pablorod/data/metagenomics/REVCO/FAST/format_gastro_reads.fasta
#521334
#pablorod@papaya:~/results/REVCO/semi_synth_dataset$ grep ">" -c /home/pablorod/data/metagenomics/REVCO/FAST/final.contigs.fa
#52953

# Coverage
RACC.uniseqDB.coverage.plot(coverage.data = "Coverage.info", output.name = "SD_coverage.png")

# All
RACC.inconsistency.full(inconsistency.data = "species_synth.inconsistency", n.reads = 521334, n.contigs = 52953, output.name = "SPECIES", filter.flag = FALSE, binary.flag = TRUE, null.tax.flag=FALSE)
RACC.inconsistency.full(inconsistency.data = "species_synth.inconsistency", n.reads = 521334, n.contigs = 52953, output.name = "SPECIES", filter.flag = TRUE, binary.flag = TRUE, null.tax.flag=FALSE)

RACC.inconsistency.full(inconsistency.data = "family_synth.inconsistency", n.reads = 521334, n.contigs = 52953, output.name = "FAMILY", filter.flag = FALSE, binary.flag = TRUE, null.tax.flag=FALSE)
RACC.inconsistency.full(inconsistency.data = "family_synth.inconsistency", n.reads = 521334, n.contigs = 52953, output.name = "FAMILY", filter.flag = TRUE, binary.flag = TRUE, null.tax.flag=FALSE)

RACC.inconsistency.full(inconsistency.data = "class_synth.inconsistency", n.reads = 521334, n.contigs = 52953, output.name = "CLASS", filter.flag = FALSE, binary.flag = TRUE, null.tax.flag=FALSE)
RACC.inconsistency.full(inconsistency.data = "class_synth.inconsistency", n.reads = 521334, n.contigs = 52953, output.name = "CLASS", filter.flag = TRUE, binary.flag = TRUE, null.tax.flag=FALSE)

# SEMI-SYNTHETIC DATASET (NOISE DISTRIBUTION 25 25 25 25)

setwd("C:/Users/Blinsky.Blinsk/Documents/- UMA/BitLab/Metagenome/Synthetic Dataset/SemiSynth/results/")
#pablorod@papaya:~/results/REVCO/semi_synth_dataset$ grep ">" -c /home/pablorod/data/metagenomics/semisynth/EXECUTION/asm/final.contigs.fa
#56004
#pablorod@papaya:~/results/REVCO/semi_synth_dataset$ grep ">" -c /home/pablorod/data/metagenomics/semisynth/FASTA/reads/ssd_reads.fa
#499991

RACC.uniseqDB.coverage.plot(coverage.data = "Coverage.info", output.name = "SD_coverage.png")

# All
RACC.inconsistency.full(inconsistency.data = "species_synth.inconsistency", n.reads = 499991, n.contigs = 56004, output.name = "SPECIES", filter.flag = FALSE, binary.flag = TRUE, null.tax.flag=FALSE)
RACC.inconsistency.full(inconsistency.data = "species_synth.inconsistency", n.reads = 499991, n.contigs = 56004, output.name = "SPECIES", filter.flag = TRUE, binary.flag = TRUE, null.tax.flag=FALSE)

RACC.inconsistency.full(inconsistency.data = "family_synth.inconsistency", n.reads = 499991, n.contigs = 56004, output.name = "FAMILY", filter.flag = FALSE, binary.flag = TRUE, null.tax.flag=FALSE)
RACC.inconsistency.full(inconsistency.data = "family_synth.inconsistency", n.reads = 499991, n.contigs = 56004, output.name = "FAMILY", filter.flag = TRUE, binary.flag = TRUE, null.tax.flag=FALSE)

RACC.inconsistency.full(inconsistency.data = "class_synth.inconsistency", n.reads = 499991, n.contigs = 56004, output.name = "CLASS", filter.flag = FALSE, binary.flag = TRUE, null.tax.flag=FALSE)
RACC.inconsistency.full(inconsistency.data = "class_synth.inconsistency", n.reads = 499991, n.contigs = 56004, output.name = "CLASS", filter.flag = TRUE, binary.flag = TRUE, null.tax.flag=FALSE)

# ---------------------

# Without None + No_hits + Not_assigned
RACC.inconsistency.full(inconsistency.data = "species_synth.inconsistency", n.reads = 499991, n.contigs = 56004, output.name = "SPECIES_NR", filter.flag = FALSE, binary.flag = TRUE, null.tax.flag=TRUE)
RACC.inconsistency.full(inconsistency.data = "species_synth.inconsistency", n.reads = 499991, n.contigs = 56004, output.name = "SPECIES_NR", filter.flag = TRUE, binary.flag = TRUE, null.tax.flag=TRUE)

RACC.inconsistency.full(inconsistency.data = "family_synth.inconsistency", n.reads = 499991, n.contigs = 56004, output.name = "FAMILY_NR", filter.flag = FALSE, binary.flag = TRUE, null.tax.flag=TRUE)
RACC.inconsistency.full(inconsistency.data = "family_synth.inconsistency", n.reads = 499991, n.contigs = 56004, output.name = "FAMILY_NR", filter.flag = TRUE, binary.flag = TRUE, null.tax.flag=TRUE)

RACC.inconsistency.full(inconsistency.data = "class_synth.inconsistency", n.reads = 499991, n.contigs = 56004, output.name = "CLASS_NR", filter.flag = FALSE, binary.flag = TRUE, null.tax.flag=TRUE)
RACC.inconsistency.full(inconsistency.data = "class_synth.inconsistency", n.reads = 499991, n.contigs = 56004, output.name = "CLASS_NR", filter.flag = TRUE, binary.flag = TRUE, null.tax.flag=TRUE)

# SEMI-SYNTHETIC DATASET V2 (NOISE DISTRIBUTION 20 30 30 20)
setwd("C:/Users/Blinsky.Blinsk/Documents/- UMA/BitLab/Metagenome/Synthetic Dataset/SemiSynth/v2/results/")
#pablorod@papaya:~/results/REVCO/semi_synth_dataset$ grep ">" -c /home/pablorod/data/metagenomics/semisynth/EXECUTION/asm/final.contigs.fa
#56004
#pablorod@papaya:~/results/REVCO/semi_synth_dataset$ grep ">" -c /home/pablorod/data/metagenomics/semisynth/FASTA/reads/ssd_reads.fa
#499991

RACC.uniseqDB.coverage.plot(coverage.data = "Coverage.info", output.name = "SD_coverage.png")

# All
RACC.inconsistency.full(inconsistency.data = "species_synth.inconsistency", n.reads = 499991, n.contigs = 56004, output.name = "SPECIES", filter.flag = FALSE, binary.flag = TRUE, null.tax.flag=FALSE)
RACC.inconsistency.full(inconsistency.data = "species_synth.inconsistency", n.reads = 499991, n.contigs = 56004, output.name = "SPECIES", filter.flag = TRUE, binary.flag = TRUE, null.tax.flag=FALSE)

RACC.inconsistency.full(inconsistency.data = "family_synth.inconsistency", n.reads = 499991, n.contigs = 56004, output.name = "FAMILY", filter.flag = FALSE, binary.flag = TRUE, null.tax.flag=FALSE)
RACC.inconsistency.full(inconsistency.data = "family_synth.inconsistency", n.reads = 499991, n.contigs = 56004, output.name = "FAMILY", filter.flag = TRUE, binary.flag = TRUE, null.tax.flag=FALSE)

RACC.inconsistency.full(inconsistency.data = "class_synth.inconsistency", n.reads = 499991, n.contigs = 56004, output.name = "CLASS", filter.flag = FALSE, binary.flag = TRUE, null.tax.flag=FALSE)
RACC.inconsistency.full(inconsistency.data = "class_synth.inconsistency", n.reads = 499991, n.contigs = 56004, output.name = "CLASS", filter.flag = TRUE, binary.flag = TRUE, null.tax.flag=FALSE)
