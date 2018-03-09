# REVCO Plot Generator
rm(list=ls())
setwd("C:/Users/Blinsky.Blinsk/Documents/- UMA/BitLab/Metagenome/Synthetic Dataset/Results")

#inconsistency_path <- "synth.inconsistencyseq"
inconsistency_path <- "assigned.inconsistency"

inconsistency_info <- read.csv(inconsistency_path, header=TRUE)

#Get taxons
read.taxons <- as.vector(unique(inconsistency_info$ReadTaxonomy))
contig.taxons <- as.vector(unique(inconsistency_info$ContigTaxonomy))


# Count inconsistencies
read.inconsistencies <- sapply(1:length(read.taxons), function(x) 
  sum(inconsistency_info$ReadTaxonomy == read.taxons[[x]]) )

contig.inconsistencies <- sapply(1:length(contig.taxons), function(x) 
  sum(inconsistency_info$ContigTaxonomy == contig.taxons[[x]]) )

names(read.inconsistencies) <- read.taxons
names(contig.inconsistencies) <- contig.taxons

read.inconsistencies <- sort(read.inconsistencies, decreasing = TRUE)
contig.inconsistencies <- sort(contig.inconsistencies, decreasing = TRUE)

#### ----------------
save(list = ls(all.names = TRUE), file="MARKDsums.RData")

#################################
load("INCsums.RData")
load("ASGsums.RData")
load("MARKDsums.RData")

colorfunc <- colorRampPalette(c("dark blue", "light blue", "light green", "dark green", "yellow", "orange", "dark red"))

###################

png(filename="all-reads.png", width = 1024, height = 1024)
barplot(read.inconsistencies, main = "Reads Inconsistencies per Taxonomy", ylab = "Inconsistencies", xlab = "Taxonomy", 
        col = colorfunc(length(read.inconsistencies)), beside = TRUE, cex.names = 0.6, log="y")
legend("topright", names(read.inconsistencies), cex = 1.5, fill = colorfunc(length(read.inconsistencies)), ncol = 2)
dev.off()

png(filename="all-contigs.png", width = 1024, height = 1024)
barplot(contig.inconsistencies, main = "Contigs Inconsistencies per Taxonomy", ylab = "Inconsistencies", xlab = "Taxonomy", 
        col = colorfunc(length(contig.inconsistencies)), beside = TRUE, cex.names = 0.6, log="y")
legend("topright", names(contig.inconsistencies), cex = 1.5, fill = colorfunc(length(contig.inconsistencies)), ncol = 2)
dev.off()

###################
### 2DO ###

#pablorod@papaya:~/data/metagenomics/REVCO$ grep ">" format_gastro_reads.fasta | wc -l
#521334
#pablorod@papaya:~/data/metagenomics/REVCO$ grep ">" contig.simple | wc -l       
#829331

#pablorod@papaya:~/software/RACC/REVCO$ wc -l ContigID-SpeciesIndex.csv
#720386 ContigID-SpeciesIndex.csv
#pablorod@papaya:~/software/RACC/REVCO$ wc -l ReadID-SpeciesIndex.csv
#455416 ReadID-SpeciesIndex.csv

# Reads Assigned
tot.colors <- c("yellow", "orange")
tot.names <- c("Consistant", "Inconsistant")
n.reads.total <- 455416
n.reads.inconsistent <- length(unique(inconsistency_info$ReadID))
reads.label = c(signif((n.reads.total-n.reads.inconsistent)*100/n.reads.total), signif(n.reads.inconsistent*100/n.reads.total))

png(filename="total-reads.png", width = 1024, height = 1024)
pie(c(n.reads.total, n.reads.inconsistent), labels = paste("Percentage: ",reads.label, "\nAmount: ", c(n.reads.total, n.reads.inconsistent)), col = tot.colors
    , main = "Reads Taxonomy Status after classification")
legend("topleft", tot.names, cex = 1.5, fill = tot.colors)
dev.off()

# Contigs Assigned
n.contigs.total <- 720386
n.contigs.inconsistent <- length(unique(inconsistency_info$ContigID))
contigs.label <-  c(signif((n.contigs.total-n.contigs.inconsistent)*100/n.contigs.total), signif(n.contigs.inconsistent*100/n.contigs.total))
png(filename="total-contigs.png", width = 1024, height = 1024)
pie(c(n.contigs.total, n.contigs.inconsistent), labels = paste("Percentage: ", contigs.label, "\nAmount: ", c(n.contigs.total, n.contigs.inconsistent)), col = tot.colors
    , main = "Contigs Taxonomy Status after classification")
legend("topleft", tot.names, cex = 1.5, fill = tot.colors)
dev.off()

################### OR #######################

reads.inconsistent.plus <- inconsistency_info[which(inconsistency_info$SemiConsistant == '+'),]$ReadID
reads.inconsistent.minus <- inconsistency_info[which(inconsistency_info$SemiConsistant == '-'),]$ReadID
n.reads.inconsistent.plus <- length(unique(reads.inconsistent.plus))
n.reads.inconsistent.minus <- length(unique(reads.inconsistent.minus))

contigs.inconsistent.plus <- inconsistency_info[which(inconsistency_info$SemiConsistant == '+'),]$ContigID
contigs.inconsistent.minus <- inconsistency_info[which(inconsistency_info$SemiConsistant == '-'),]$ContigID
n.contigs.inconsistent.plus <- length(unique(contigs.inconsistent.plus))
n.contigs.inconsistent.minus <- length(unique(contigs.inconsistent.minus))

# Reads Assigned
tot.colors <- c("yellow", "orange", "red")
tot.names <- c("Consistant", "Semi-Inconsistant", "Inconsistant")
n.reads.total <- 455416
reads.label = c(signif((n.reads.total-n.reads.inconsistent)*100/n.reads.total), signif(n.reads.inconsistent.plus*100/n.reads.total)
                , signif(n.reads.inconsistent.minus*100/n.reads.total))

png(filename="marked-total-reads.png", width = 1024, height = 1024)
pie(c(n.reads.total, n.reads.inconsistent.plus, n.reads.inconsistent.minus), labels = paste("Percentage: ",reads.label, "\nAmount: ", c(n.reads.total, n.reads.inconsistent.plus, n.reads.inconsistent.minus)), col = tot.colors
    , main = "Reads Taxonomy Status after classification")
legend("topleft", tot.names, cex = 1.5, fill = tot.colors)
dev.off()

# Contigs Assigned
n.contigs.total <- 720386
n.contigs.inconsistent <- length(unique(inconsistency_info$ContigID))
contigs.label <-  c(signif((n.contigs.total-n.contigs.inconsistent)*100/n.contigs.total), signif(n.contigs.inconsistent.plus*100/n.contigs.total)
                    , signif(n.contigs.inconsistent.minus*100/n.reads.total))

png(filename="marked-total-contigs.png", width = 1024, height = 1024)
pie(c(n.contigs.total, n.contigs.inconsistent.plus, n.contigs.inconsistent.minus), labels = paste("Percentage: ", contigs.label,
    "\nAmount: ", c(n.contigs.total, n.contigs.inconsistent.plus, n.contigs.inconsistent.minus)), col = tot.colors
    , main = "Contigs Taxonomy Status after classification")
legend("topleft", tot.names, cex = 1.5, fill = tot.colors)
dev.off()


##########################
paste(names(read.inconsistencies), "\n", read.inconsistencies)

pie(read.inconsistencies, col = colorfunc(length(read.inconsistencies)),  
    labels = "", main = "Reads Inconsistencies per Taxonomy")
legend("topright", names(read.inconsistencies), cex = 0.3, fill = colorfunc(length(read.inconsistencies)), ncol = 2)

pie(contig.inconsistencies, col = colorfunc(length(contig.inconsistencies)), 
    labels = "", main = "Contigs Inconsistencies per Taxonomy")
legend("topright", names(contig.inconsistencies), cex = 0.3, fill = colorfunc(length(contig.inconsistencies)), ncol = 2)


############ NONE REMOVED ############
full.contig.inconsistencies
pie(assigned.read.inconsistencies, col = colorfunc(length(assigned.read.inconsistencies)), 
    labels = "", main = "Removed None Reads - Reads Inconsistencies per Taxonomy")
#legend("bottomleft", names(assigned.read.inconsistencies), cex = 0.5, fill = colorfunc(length(full.read.inconsistencies)), ncol = 2)

pie(assigned.contig.inconsistencies, col = colorfunc(length(assigned.contig.inconsistencies)), 
    labels = "", main = "Removed None Reads - Contigs Inconsistencies per Taxonomy")
#legend("bottomleft", names(assigned.contig.inconsistencies), cex = 0.5, fill = colorfunc(length(full.read.inconsistencies)), ncol = 2)


assigned.contig.inconsistencies

###################