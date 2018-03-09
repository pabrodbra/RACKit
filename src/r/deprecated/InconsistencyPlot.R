rm(list=ls())
setwd("C:/Users/Blinsky.Blinsk/Documents/- UMA/BitLab/Metagenome/REVCO")

# Load data
forward.full <- read.csv("forward-REVCO.inconsistency", header=TRUE)
reverse.full <- read.csv("reverse-REVCO.inconsistency", header=TRUE)
forward.assigned <- read.csv("forward-REVCO-assigned.inconsistency", header=TRUE)
reverse.assigned <- read.csv("reverse-REVCO-assigned.inconsistency", header=TRUE)

#Get taxons
read.full.list.taxons <- as.vector(unique(forward.full$ReadTaxonomy))
read.assigned.list.taxons <- as.vector(unique(forward.assigned$ReadTaxonomy))


contig.full.list.taxons <- as.vector(unique(reverse.full$ContigTaxonomy))
contig.assigned.list.taxons <- as.vector(unique(forward.assigned$ContigTaxonomy))

# Count inconsistencies
full.read.inconsistencies <- sapply(1:length(read.full.list.taxons), function(x) 
        sum(forward.full$ReadTaxonomy == read.full.list.taxons[[x]]) + sum(reverse.full$ReadTaxonomy == read.full.list.taxons[[x]]))

full.contig.inconsistencies <- sapply(1:length(contig.full.list.taxons), function(x) 
  sum(forward.full$ContigTaxonomy == contig.full.list.taxons[[x]]) + sum(reverse.full$ContigTaxonomy == contig.full.list.taxons[[x]]))

assigned.read.inconsistencies <- sapply(1:length(read.assigned.list.taxons), function(x) 
  sum(forward.assigned$ReadTaxonomy == read.assigned.list.taxons[[x]]) + sum(reverse.assigned$ReadTaxonomy == read.assigned.list.taxons[[x]])) 

assigned.contig.inconsistencies <- sapply(1:length(contig.assigned.list.taxons), function(x) 
  sum(forward.assigned$ContigTaxonomy == contig.assigned.list.taxons[[x]]) + sum(reverse.assigned$ContigTaxonomy == contig.assigned.list.taxons[[x]])) 

names(full.read.inconsistencies) <- read.full.list.taxons
names(full.contig.inconsistencies) <- contig.full.list.taxons
names(assigned.read.inconsistencies) <- read.assigned.list.taxons
names(assigned.contig.inconsistencies) <- contig.assigned.list.taxons

full.read.inconsistencies <- sort(full.read.inconsistencies, decreasing = TRUE)
full.contig.inconsistencies <- sort(full.contig.inconsistencies, decreasing = TRUE)
assigned.read.inconsistencies <- sort(assigned.read.inconsistencies, decreasing = TRUE)
assigned.contig.inconsistencies <- sort(assigned.contig.inconsistencies, decreasing = TRUE)

save(list = ls(all.names = TRUE), file="sums.RData")
#################################
load("sums.RData")
# colorfunc <- colorRampPalette(sample(colors(), length(assigned.contig.inconsistencies)))
colorfunc <- colorRampPalette(c("dark blue", "light blue", "light green", "dark green", "yellow", "orange", "dark red"))

png(filename="all-reads.png", width = 960, height = 960)
barplot(full.read.inconsistencies, main = "Reads Inconsistencies per Taxonomy", ylab = "Inconsistencies", xlab = "Taxonomy", 
        col = colorfunc(length(full.read.inconsistencies)), beside = TRUE, cex.names = 0.6)
legend("topright", names(full.read.inconsistencies), cex = 1.5, fill = colorfunc(length(full.read.inconsistencies)), ncol = 2)
dev.off()

png(filename="all-contigs.png", width = 960, height = 960)
barplot(full.contig.inconsistencies, main = "Contigs Inconsistencies per Taxonomy", ylab = "Inconsistencies", xlab = "Taxonomy", 
        col = colorfunc(length(full.contig.inconsistencies)), beside = TRUE, cex.names = 0.6)
legend("topright", names(full.contig.inconsistencies), cex = 1.5, fill = colorfunc(length(full.contig.inconsistencies)), ncol = 2)
dev.off()

png(filename="nn-reads.png", width = 960, height = 960)
barplot(assigned.read.inconsistencies, main = "Removed None Reads - Reads Inconsistencies per Taxonomy", ylab = "Inconsistencies", xlab = "Taxonomy", 
        col = colorfunc(length(assigned.read.inconsistencies)), beside = TRUE, cex.names = 0.6)
legend("topright", names(assigned.read.inconsistencies), cex = 1.5, fill = colorfunc(length(assigned.read.inconsistencies)), ncol = 2)
dev.off()

png(filename="nn-contigs.png", width = 960, height = 960)
barplot(assigned.contig.inconsistencies, main = "Removed None Reads - Contigs Inconsistencies per Taxonomy", ylab = "Inconsistencies", xlab = "Taxonomy", 
        col = colorfunc(length(assigned.contig.inconsistencies)), beside = TRUE, cex.names = 0.6)
legend("topright", names(assigned.contig.inconsistencies), cex = 1.5, fill = colorfunc(length(assigned.contig.inconsistencies)), ncol = 2)
dev.off()

#########################

#pablorod@mango:~/metagenome$ grep ">" reverse_paired.fasta | wc -l
#4177249
#pablorod@mango:~/metagenome$ grep ">" forward_paired.fasta | wc -l
#4177249
#pablorod@mango:~/metagenome$ grep ">" contigDB/asm.contig | wc -l
#168703

# Reads Assigned
tot.colors <- c("yellow", "orange")
tot.names <- c("Consistant", "Inconsistant")
n.reads.total <- 4177249 + 4177249
n.reads.inconsistent <- length(unique(forward.full$ReadID)) + length(unique(reverse.full$ReadID))
reads.label = c(signif((n.reads.total-n.reads.inconsistent)*100/n.reads.total), signif(n.reads.inconsistent*100/n.reads.total))

png(filename="total-reads.png", width = 960, height = 960)
pie(c(n.reads.total, n.reads.inconsistent), labels = paste("Percentage: ",reads.label, "\nAmount: ", c(n.reads.total, n.reads.inconsistent)), col = tot.colors
    , main = "Reads Taxonomy Status after classification")
legend("topleft", tot.names, cex = 1.5, fill = tot.colors)
dev.off()

# Contigs Assigned
n.contigs.total <- 168703
n.contigs.inconsistent <- length(unique(forward.full$ContigID)) + length(unique(reverse.full$ContigID))
contigs.label <-  c(signif((n.contigs.total-n.contigs.inconsistent)*100/n.contigs.total), signif(n.contigs.inconsistent*100/n.contigs.total))
png(filename="total-contigs.png", width = 960, height = 960)
pie(c(n.contigs.total, n.contigs.inconsistent), labels = paste("Percentage: ", contigs.label, "\nAmount: ", c(n.contigs.total, n.contigs.inconsistent)), col = tot.colors
    , main = "Contigs Taxonomy Status after classification")
legend("topleft", tot.names, cex = 1.5, fill = tot.colors)
dev.off()

##########################
paste(names(full.read.inconsistencies), "\n", full.read.inconsistencies)

pie(full.read.inconsistencies, col = colorfunc(length(full.read.inconsistencies)),  
    labels = "", main = "Reads Inconsistencies per Taxonomy")
legend("bottomleft", names(full.read.inconsistencies), cex = 0.5, fill = colorfunc(length(full.read.inconsistencies)), ncol = 2)

pie(full.contig.inconsistencies, col = colorfunc(length(full.contig.inconsistencies)), 
    labels = "", main = "Contigs Inconsistencies per Taxonomy")
legend("bottomleft", names(full.contig.inconsistencies), cex = 0.5, fill = colorfunc(length(full.read.inconsistencies)), ncol = 2)
full.contig.inconsistencies
pie(assigned.read.inconsistencies, col = colorfunc(length(assigned.read.inconsistencies)), 
    labels = "", main = "Removed None Reads - Reads Inconsistencies per Taxonomy")
legend("bottomleft", names(assigned.read.inconsistencies), cex = 0.5, fill = colorfunc(length(full.read.inconsistencies)), ncol = 2)

pie(assigned.contig.inconsistencies, col = colorfunc(length(assigned.contig.inconsistencies)), 
    labels = "", main = "Removed None Reads - Contigs Inconsistencies per Taxonomy")
legend("bottomleft", names(assigned.contig.inconsistencies), cex = 0.5, fill = colorfunc(length(full.read.inconsistencies)), ncol = 2)


assigned.contig.inconsistencies
