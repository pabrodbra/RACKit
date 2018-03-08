# REVCO Plot Generator - CONST, SEMICONST, INCONST

rm(list=ls())
setwd("C:/Users/Blinsky.Blinsk/Documents/- UMA/BitLab/Metagenome/Synthetic Dataset/Results")

inconsistency_path <- "synth_marked.inconsistency"
assigned_inconsistency_path <- "marked_assigned.inconsistency"

inconsistency_path <- "species_synth.inconsistency"
inconsistency_info <- read.csv(inconsistency_path, header=TRUE)
assigned.inconsistency_info <- read.csv(assigned_inconsistency_path, header=TRUE)

#Get taxons
read.taxons <- as.vector(unique(inconsistency_info$ReadTaxonomy))
contig.taxons <- as.vector(unique(inconsistency_info$ContigTaxonomy))

assigned.read.taxons <- as.vector(unique(assigned.inconsistency_info$ReadTaxonomy))
assigned.contig.taxons <- as.vector(unique(assigned.inconsistency_info$ContigTaxonomy))

# Count inconsistencies
read.inconsistencies <- sapply(1:length(read.taxons), function(x) 
  sum(inconsistency_info$ReadTaxonomy == read.taxons[[x]]) )

contig.inconsistencies <- sapply(1:length(contig.taxons), function(x) 
  sum(inconsistency_info$ContigTaxonomy == contig.taxons[[x]]) )

names(read.inconsistencies) <- read.taxons
names(contig.inconsistencies) <- contig.taxons

read.inconsistencies <- sort(read.inconsistencies, decreasing = TRUE)
contig.inconsistencies <- sort(contig.inconsistencies, decreasing = TRUE)

#-
assigned.read.inconsistencies <- sapply(1:length(assigned.read.taxons), function(x) 
  sum(assigned.inconsistency_info$ReadTaxonomy == assigned.read.taxons[[x]]) )

assigned.contig.inconsistencies <- sapply(1:length(assigned.contig.taxons), function(x) 
  sum(assigned.inconsistency_info$ContigTaxonomy == assigned.contig.taxons[[x]]) )

names(assigned.read.inconsistencies) <- assigned.read.taxons
names(assigned.contig.inconsistencies) <- assigned.contig.taxons

assigned.read.inconsistencies <- sort(assigned.read.inconsistencies, decreasing = TRUE)
assigned.contig.inconsistencies <- sort(assigned.contig.inconsistencies, decreasing = TRUE)
#### ----------------
save(list = ls(all.names = TRUE), file="MARKDsums.RData")

load("MARKDsums.RData")

colorfunc <- colorRampPalette(c("dark blue", "light blue", "light green", "dark green", "yellow", "orange", "dark red"))

###################
png(filename="all-reads.png", width = 1024, height = 1024)
barplot(read.inconsistencies, main = "Reads Inconsistencies per Taxonomy", ylab = "Inconsistencies", xlab = "Taxonomy", 
        col = colorfunc(length(read.inconsistencies)), beside = TRUE, cex.names = 0.6, log="y")
legend("topright", names(read.inconsistencies)[1:10], cex = 1.5, fill = colorfunc(length(read.inconsistencies)), ncol = 2)
dev.off()

png(filename="all-contigs.png", width = 1024, height = 1024)
barplot(contig.inconsistencies, main = "Contigs Inconsistencies per Taxonomy", ylab = "Inconsistencies", xlab = "Taxonomy", 
        col = colorfunc(length(contig.inconsistencies)), beside = TRUE, cex.names = 0.6, log="y")
legend("topright", names(contig.inconsistencies)[1:10], cex = 1.5, fill = colorfunc(length(contig.inconsistencies)), ncol = 2)
dev.off()

png(filename="assigned-reads.png", width = 1024, height = 1024)
barplot(assigned.read.inconsistencies, main = "Assigned Reads Inconsistencies per Taxonomy", ylab = "Inconsistencies", xlab = "Taxonomy", 
        col = colorfunc(length(assigned.read.inconsistencies)), beside = TRUE, cex.names = 0.6, log="y")
legend("topright", names(assigned.read.inconsistencies)[1:10], cex = 1.5, fill = colorfunc(length(assigned.read.inconsistencies)), ncol = 2)
dev.off()
png(filename="assigned-contigs.png", width = 1024, height = 1024)
barplot(assigned.contig.inconsistencies, main = "Assigned Contigs Inconsistencies per Taxonomy", ylab = "Inconsistencies", xlab = "Taxonomy", 
        col = colorfunc(length(assigned.contig.inconsistencies)), beside = TRUE, cex.names = 0.6, log="y")
legend("topright", names(assigned.contig.inconsistencies)[1:10], cex = 1.5, fill = colorfunc(length(assigned.contig.inconsistencies)), ncol = 2)
dev.off()

# filtered
filtered.read.inconsistencies <- read.inconsistencies[which(read.inconsistencies > sum(read.inconsistencies)/100)]
filtered.contigs.inconsistencies <-  contig.inconsistencies[which(contig.inconsistencies > sum(contig.inconsistencies)/100)]
assigned.filtered.reads.inconsistencies <- assigned.read.inconsistencies[which(assigned.read.inconsistencies > sum(assigned.read.inconsistencies)/100)]
assigned.filtered.contigs.inconsistencies <- assigned.contig.inconsistencies[which(assigned.contig.inconsistencies > sum(assigned.contig.inconsistencies)/100)]

png(filename="filtered_all-reads.png", width = 1024, height = 1024)
barplot(filtered.read.inconsistencies, main = "Reads Inconsistencies per Taxonomy", ylab = "Inconsistencies", xlab = "Taxonomy", 
        col = colorfunc(length(filtered.read.inconsistencies)), beside = TRUE, cex.names = 0.6, log="y")
legend("topright", names(filtered.read.inconsistencies), cex = 1.5, fill = colorfunc(length(filtered.read.inconsistencies)), ncol = 2)
dev.off()

png(filename="filtered_all-contigs.png", width = 1024, height = 1024)
barplot(filtered.contigs.inconsistencies, main = "Contigs Inconsistencies per Taxonomy", ylab = "Inconsistencies", xlab = "Taxonomy", 
        col = colorfunc(length(filtered.contigs.inconsistencies)), beside = TRUE, cex.names = 0.6, log="y")
legend("topright", names(filtered.contigs.inconsistencies), cex = 1.5, fill = colorfunc(length(filtered.contigs.inconsistencies)), ncol = 2)
dev.off()

png(filename="filtered_assigned-reads.png", width = 1024, height = 1024)
barplot(assigned.filtered.reads.inconsistencies, main = "Assigned Reads Inconsistencies per Taxonomy", ylab = "Inconsistencies", xlab = "Taxonomy", 
        col = colorfunc(length(assigned.filtered.reads.inconsistencies)), beside = TRUE, cex.names = 0.6, log="y")
legend("topright", names(assigned.filtered.reads.inconsistencies), cex = 1.5, fill = colorfunc(length(assigned.filtered.reads.inconsistencies)), ncol = 2)
dev.off()

png(filename="filtered_assigned-contigs.png", width = 1024, height = 1024)
barplot(assigned.filtered.contigs.inconsistencies, main = "Assigned Contigs Inconsistencies per Taxonomy", ylab = "Inconsistencies", xlab = "Taxonomy", 
        col = colorfunc(length(assigned.filtered.contigs.inconsistencies)), beside = TRUE, cex.names = 0.6, log="y")
legend("topright", names(assigned.filtered.contigs.inconsistencies), cex = 1.5, fill = colorfunc(length(assigned.filtered.contigs.inconsistencies)), ncol = 2)
dev.off()
###################

#pablorod@papaya:~/data/metagenomics/REVCO$ grep ">" format_gastro_reads.fasta | wc -l
#521334
#pablorod@papaya:~/data/metagenomics/REVCO$ grep ">" contig.simple | wc -l       
#829331

#pablorod@papaya:~/software/RACC/REVCO$ wc -l ContigID-SpeciesIndex.csv
#720386 ContigID-SpeciesIndex.csv
#pablorod@papaya:~/software/RACC/REVCO$ wc -l ReadID-SpeciesIndex.csv
#455416 ReadID-SpeciesIndex.csv

tot.colors <- c("yellow", "orange", "red")
tot.names <- c("Consistent", "Semi-Inconsistent", "Inconsistent")

## All
reads.inconsistent.plus <- inconsistency_info[which(inconsistency_info$SemiConsistant == '+'),]$ReadID
reads.inconsistent.minus <- inconsistency_info[which(inconsistency_info$SemiConsistant == '-'),]$ReadID
n.reads.inconsistent.plus <- length(unique(reads.inconsistent.plus))
n.reads.inconsistent.minus <- length(unique(reads.inconsistent.minus))

contigs.inconsistent.plus <- inconsistency_info[which(inconsistency_info$SemiConsistant == '+'),]$ContigID
contigs.inconsistent.minus <- inconsistency_info[which(inconsistency_info$SemiConsistant == '-'),]$ContigID
n.contigs.inconsistent.plus <- length(unique(contigs.inconsistent.plus))
n.contigs.inconsistent.minus <- length(unique(contigs.inconsistent.minus))

# Reads
n.reads.total <- 455416
n.reads.consistant <- n.reads.total-n.reads.inconsistent.plus-n.reads.inconsistent.minus
reads.label = c(signif(n.reads.consistant*100/n.reads.total), signif(n.reads.inconsistent.plus*100/n.reads.total)
                , signif(n.reads.inconsistent.minus*100/n.reads.total))

png(filename="marked-total-reads.png", width = 1024, height = 1024)
pie(c(n.reads.consistant, n.reads.inconsistent.plus, n.reads.inconsistent.minus), labels = paste("Percentage: ",reads.label, 
    "\nAmount: ", c(n.reads.consistant, n.reads.inconsistent.plus, n.reads.inconsistent.minus)), col = tot.colors
    , main = "Reads Taxonomy Status after classification")
legend("topleft", tot.names, cex = 1.5, fill = tot.colors)
dev.off()

# Contigs
n.contigs.total <- 720386
n.contigs.inconsistent <- length(unique(inconsistency_info$ContigID))
n.contigs.consistant <- n.contigs.total-n.contigs.inconsistent.plus-n.contigs.inconsistent.minus
contigs.label <-  c(signif(n.contigs.consistant*100/n.contigs.total), signif(n.contigs.inconsistent.plus*100/n.contigs.total)
                    , signif(n.contigs.inconsistent.minus*100/n.reads.total))

png(filename="marked-total-contigs.png", width = 1024, height = 1024)
pie(c(n.contigs.consistant, n.contigs.inconsistent.plus, n.contigs.inconsistent.minus), labels = paste("Percentage: ", contigs.label,
    "\nAmount: ", c(n.contigs.consistant, n.contigs.inconsistent.plus, n.contigs.inconsistent.minus)), col = tot.colors, 
    main = "Contigs Taxonomy Status after classification")
legend("topleft", tot.names, cex = 1.5, fill = tot.colors)
dev.off()

## Assigned

assigned.reads.inconsistent.plus <- assigned.inconsistency_info[which(assigned.inconsistency_info$SemiConsistant == '+'),]$ReadID
assigned.reads.inconsistent.minus <- assigned.inconsistency_info[which(assigned.inconsistency_info$SemiConsistant == '-'),]$ReadID
assigned.n.reads.inconsistent.plus <- length(unique(assigned.reads.inconsistent.plus))
assigned.n.reads.inconsistent.minus <- length(unique(assigned.reads.inconsistent.minus))

assigned.contigs.inconsistent.plus <- assigned.inconsistency_info[which(assigned.inconsistency_info$SemiConsistant == '+'),]$ContigID
assigned.contigs.inconsistent.minus <- assigned.inconsistency_info[which(assigned.inconsistency_info$SemiConsistant == '-'),]$ContigID
assigned.n.contigs.inconsistent.plus <- length(unique(assigned.contigs.inconsistent.plus))
assigned.n.contigs.inconsistent.minus <- length(unique(assigned.contigs.inconsistent.minus))

# Reads
assigned.n.reads.consistant <- n.reads.total-assigned.n.reads.inconsistent.plus-assigned.n.reads.inconsistent.minus
assigned.reads.label = c(signif(assigned.n.reads.consistant*100/n.reads.total), signif(assigned.n.reads.inconsistent.plus*100/n.reads.total)
                         , signif(assigned.n.reads.inconsistent.minus*100/n.reads.total))

png(filename="assigned_marked-total-reads.png", width = 1024, height = 1024)
pie(c(assigned.n.reads.consistant, assigned.n.reads.inconsistent.plus, assigned.n.reads.inconsistent.minus), labels = paste("Percentage: ",assigned.reads.label, 
                                                                                                                            "\nAmount: ", c(assigned.n.reads.consistant, assigned.n.reads.inconsistent.plus, assigned.n.reads.inconsistent.minus)), col = tot.colors
    , main = "Assigned Reads Taxonomy Status after classification")
legend("topleft", tot.names, cex = 1.5, fill = tot.colors)
dev.off()

# Contigs
assigned.n.contigs.consistant <- n.contigs.total-assigned.n.contigs.inconsistent.plus-assigned.n.contigs.inconsistent.minus
assigned.contigs.label <-  c(signif(assigned.n.contigs.consistant*100/n.contigs.total), signif(assigned.n.contigs.inconsistent.plus*100/n.contigs.total)
                             , signif(assigned.n.contigs.inconsistent.minus*100/n.reads.total))

png(filename="assigned_marked-total-contigs.png", width = 1024, height = 1024)
pie(c(assigned.n.contigs.consistant, assigned.n.contigs.inconsistent.plus, assigned.n.contigs.inconsistent.minus), labels = paste("Percentage: ", assigned.contigs.label,
                                                                                                                                  "\nAmount: ", c(assigned.n.contigs.consistant, assigned.n.contigs.inconsistent.plus, assigned.n.contigs.inconsistent.minus)), col = tot.colors
    , main = "Assigned Contigs Taxonomy Status after classification")
legend("topleft", tot.names, cex = 1.5, fill = tot.colors)
dev.off()

# Mark Nones

tot.colors <- c("yellow", "blue","orange", "red")
tot.names <- c("Consistent", "Unassigned", "Semi-Inconsistent", "Inconsistent")
none.removed.reads <- n.reads.inconsistent.minus+n.reads.inconsistent.plus - (assigned.n.reads.inconsistent.plus+assigned.n.reads.inconsistent.minus)
none.removed.contigs <- n.contigs.inconsistent.minus+n.contigs.inconsistent.plus - (assigned.n.contigs.inconsistent.plus+assigned.n.contigs.inconsistent.minus)

# Reads
assigned.n.reads.consistant <- n.reads.total-assigned.n.reads.inconsistent.plus-assigned.n.reads.inconsistent.minus-none.removed.reads
assigned.reads.label = c(signif(assigned.n.reads.consistant*100/n.reads.total), signif(none.removed.reads*100/n.reads.total), signif(assigned.n.reads.inconsistent.plus*100/n.reads.total)
                         , signif(assigned.n.reads.inconsistent.minus*100/n.reads.total))

png(filename="nr_assigned_marked-total-reads.png", width = 1024, height = 1024)
pie(c(assigned.n.reads.consistant, none.removed.reads, assigned.n.reads.inconsistent.plus, assigned.n.reads.inconsistent.minus), labels = paste("Percentage: ",assigned.reads.label, 
    "\nAmount: ", c(assigned.n.reads.consistant, none.removed.reads, assigned.n.reads.inconsistent.plus, assigned.n.reads.inconsistent.minus)), col = tot.colors
    , main = "Assigned Reads Taxonomy Status after classification")
legend("topleft", tot.names, cex = 1.5, fill = tot.colors)
dev.off()

# Contigs
assigned.n.contigs.consistant <- n.contigs.total-assigned.n.contigs.inconsistent.plus-assigned.n.contigs.inconsistent.minus-none.removed.contigs
assigned.contigs.label <-  c(signif(assigned.n.contigs.consistant*100/n.contigs.total), signif(none.removed.contigs*100/n.contigs.total), signif(assigned.n.contigs.inconsistent.plus*100/n.contigs.total)
                             , signif(assigned.n.contigs.inconsistent.minus*100/n.reads.total))

png(filename="nr_assigned_marked-total-contigs.png", width = 1024, height = 1024)
pie(c(assigned.n.contigs.consistant, none.removed.contigs, assigned.n.contigs.inconsistent.plus, assigned.n.contigs.inconsistent.minus), labels = paste("Percentage: ", assigned.contigs.label,
    "\nAmount: ", c(assigned.n.contigs.consistant, none.removed.contigs, assigned.n.contigs.inconsistent.plus, assigned.n.contigs.inconsistent.minus)), col = tot.colors
    , main = "Assigned Contigs Taxonomy Status after classification")
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