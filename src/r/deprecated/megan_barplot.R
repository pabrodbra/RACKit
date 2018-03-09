rm(list = ls())

### Prepare Dataset

reads.csv <- read.csv("C:/Users/Blinsky.Blinsk/Documents/- UMA/BitLab/Metagenome/MEGAN/3 - Family/BLAST_paired-ex family.txt", header = FALSE)
reads.csv <- reads.csv[order(reads.csv$V1),]
#contigs.csv <-  read.csv("C:/Users/Blinsky.Blinsk/Documents/- UMA/BitLab/Metagenome/MEGAN/3 - Family/BLAST_SdN_contig-ex family.txt", header = FALSE)
contigs.csv <-  read.csv("C:/Users/Blinsky.Blinsk/Documents/- UMA/BitLab/Metagenome/Scripts/Family-SumTax.csv-normalized", header = FALSE)

#reads.csv <- read.csv("C:/Users/Blinsky.Blinsk/Documents/- UMA/BitLab/Metagenome/MEGAN/2 - Genus/BLAST_paired-ex.txt", header = FALSE)
#contigs1.csv <- read.csv("C:/Users/Blinsky.Blinsk/Documents/- UMA/BitLab/Metagenome/MEGAN/2 - Genus/BLAST_SdN_contig.txt", header = FALSE)
#contigs.csv <-  read.csv("C:/Users/Blinsky.Blinsk/Documents/- UMA/BitLab/Metagenome/Scripts/Genus-SumTax.csv-normalized", header = FALSE)

common.names <- intersect(contigs.csv[,1], reads.csv[,1])

reads.data <- data.frame(reads.csv[,2])
rownames(reads.data) <- reads.csv[,1]
colnames(reads.data) <- "Percentages"

contigs.data <- data.frame(contigs.csv[,2])
rownames(contigs.data) <- contigs.csv[,1]
colnames(reads.data) <- "Percentages"

logical.reads.common <- rownames(reads.data) %in% common.names
logical.contigs.common <- rownames(contigs.data) %in% common.names

reads.common <- reads.data[logical.reads.common,]
contigs.common <- contigs.data[logical.contigs.common,]
names(contigs.common) <- common.names

common.data <- rbind(reads.common, contigs.common)
rownames(common.data) <- c("Reads", "Contigs")

### Plotting
common.max_threshold <- max(common.data)

t1<- sort(common.data, decreasing=TRUE)
common.data
t1
png(filename="classified-taxonomy.png", width = 960, height = 960)
barplot(common.data, beside = TRUE, col = c("blue", "red"), cex.names = 0.4, ylim = c(0, common.max_threshold), ylab="Percentage of Sequences", 
        xlab = "Name of Taxonomy", main = "Percentage of Classified Taxonomy")
legend("topright", rownames(common.data), cex = 1.5, fill = c("blue", "red"), ncol = 2)
dev.off()

png(filename="classified-taxonomy-difference.png", width = 960, height = 960)
barplot(reads.common-contigs.common, beside = TRUE, cex.names = 0.4, ylab="Percentage (Reads-Contigs)", 
        xlab = "Name of Taxonomy", main = "Percentage Difference of Classified Taxonomy", col = colorfunc(length(common.data)))
legend("topright", colnames(common.data), cex = 1.5, fill = colorfunc(length(common.data)), ncol = 2)
dev.off()
 