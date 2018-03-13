# RACC_inconsistency_solver
# Pablo Rodriguez Brazzarola

rm(list=ls())
dev.off()

library(ggplot2)
library(reshape2)

###################################
### ---- FULL-SYNTH DATASET --- ###
###################################
setwd("C:/Users/Blinsky.Blinsk/Documents/- UMA/BitLab/Metagenome/Use Cases/Synthetic Dataset/Results")

fsd.solved.data <- read.csv("species_fsd_solved_inconsistencies.csv",header = TRUE, sep = ',')

solved.data <- fsd.solved.data
rank.names <- c("Species","Genus","Family","Order","Class","Phylum","Kingdom","Domain", "Organism")

# Where are inconsistencies solved
rank.solved.inconsistency <- sort(unique(solved.data$RankSolved), decreasing=TRUE)
count.rank.solved.inconsistency <- sapply(1:length(rank.solved.inconsistency), function(x)
  sum(solved.data$RankSolved == rank.solved.inconsistency[[x]]))
names(count.rank.solved.inconsistency) <- rank.names
count.rank.solved.inconsistency

# Summed where are inconsistencies solved
sum.solved.ranks <- sapply(1:length(count.rank.solved.inconsistency), function(x)
  sum(count.rank.solved.inconsistency[1:x]))
names(sum.solved.ranks) <- rank.names
sum.solved.ranks

# Type of inconsistencies
types.possible <- c("WR","WC","H")
count.types <- sapply(1:length(types.possible), function(x)
  sum(solved.data$Type == types.possible[[x]]))
names(count.types) <- types.possible
count.types

# Where are inconsistencies solved for each type
solved.inconsistencies.WR <- sapply(1:length(rank.solved.inconsistency), function(x)
  sum(solved.data$RankSolved[which(solved.data$Type == "WR")] == rank.solved.inconsistency[[x]]))
solved.inconsistencies.WC <- sapply(1:length(rank.solved.inconsistency), function(x)
  sum(solved.data$RankSolved[which(solved.data$Type == "WC")] == rank.solved.inconsistency[[x]]))
solved.inconsistencies.H <- sapply(1:length(rank.solved.inconsistency), function(x)
  sum(solved.data$RankSolved[which(solved.data$Type == "H")] == rank.solved.inconsistency[[x]]))


solved.inconsistencies.per.type <- rbind(solved.inconsistencies.WR, solved.inconsistencies.WC, solved.inconsistencies.H)
colnames(solved.inconsistencies.per.type) <- rank.names
solved.inconsistencies.per.type

# summed where inconsistencies are solved per type
sum.solved.ranks.WR <- sapply(1:length(count.rank.solved.inconsistency), function(x)
  sum(solved.inconsistencies.per.type[1,1:x]))
sum.solved.ranks.WC <- sapply(1:length(count.rank.solved.inconsistency), function(x)
  sum(solved.inconsistencies.per.type[2,1:x]))
sum.solved.ranks.H <- sapply(1:length(count.rank.solved.inconsistency), function(x)
  sum(solved.inconsistencies.per.type[3,1:x]))

sum.solved.ranks.per.type <- rbind(sum.solved.ranks.WR,sum.solved.ranks.WC,sum.solved.ranks.H)
colnames(sum.solved.ranks.per.type) <- rank.names
sum.solved.ranks.per.type# <- sum.solved.ranks.per.type[,2:ncol(sum.solved.ranks.per.type)]

normalize.sum.solved.ranks.per.type <- sum.solved.ranks.per.type*100/c(sum.solved.ranks.per.type[,ncol(sum.solved.ranks.per.type)])  #as.matrix(sum.solved.ranks.per.type)

melt.sum.solved.ranks.per.type <-melt(normalize.sum.solved.ranks.per.type,id.vars=names(normalize.sum.solved.ranks.per.type))

ggplot(melt.sum.solved.ranks.per.type,aes(x=factor(Var2),y=value,fill=factor(Var1)))+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_brewer(palette = "Set1", name="Inconsistency Type",
                    labels=c("Weak Reads", "Weak Contigs", "Hard"))+
  xlab("Taxonomic Rank")+
  ylab("Percentage of Inconsistencies Solved")+
  ggtitle("Taxonomic rank where Inconsistencies are Solved (FSD)", subtitle = "Percentage for each Type")+
  theme(axis.title = element_text(size=22),
        axis.text = element_text(size=18),
        plot.title =element_text(size=25,face="bold"),
        plot.subtitle = element_text(size=18),
        legend.text=element_text(size=15),
        legend.title =element_text(size=18))
  ###################################
### ---- SEMI-SYNTH DATASET --- ###
###################################
setwd("C:/Users/Blinsky.Blinsk/Documents/- UMA/BitLab/Metagenome/Use Cases/SemiSynth Dataset/results/")

ssd.solved.data <- read.csv("species_ssd_solved_inconsistencies.csv",header = TRUE, sep = ',')

solved.data <- ssd.solved.data
rank.names <- c("Species","Genus","Family","Order","Class","Phylum","Kingdom","Domain")

# Where are inconsistencies solved
rank.solved.inconsistency <- sort(unique(solved.data$RankSolved), decreasing=TRUE)
count.rank.solved.inconsistency <- sapply(1:length(rank.solved.inconsistency), function(x)
  sum(solved.data$RankSolved == rank.solved.inconsistency[[x]]))
names(count.rank.solved.inconsistency) <- rank.names
count.rank.solved.inconsistency

# Summed where are inconsistencies solved
sum.solved.ranks <- sapply(1:length(count.rank.solved.inconsistency), function(x)
  sum(count.rank.solved.inconsistency[1:x]))
names(sum.solved.ranks) <- rank.names
sum.solved.ranks

# Type of inconsistencies
types.possible <- c("WR","WC","H")
count.types <- sapply(1:length(types.possible), function(x)
  sum(solved.data$Type == types.possible[[x]]))
names(count.types) <- types.possible
count.types

# Where are inconsistencies solved for each type
solved.inconsistencies.WR <- sapply(1:length(rank.solved.inconsistency), function(x)
  sum(solved.data$RankSolved[which(solved.data$Type == "WR")] == rank.solved.inconsistency[[x]]))
solved.inconsistencies.WC <- sapply(1:length(rank.solved.inconsistency), function(x)
  sum(solved.data$RankSolved[which(solved.data$Type == "WC")] == rank.solved.inconsistency[[x]]))
solved.inconsistencies.H <- sapply(1:length(rank.solved.inconsistency), function(x)
  sum(solved.data$RankSolved[which(solved.data$Type == "H")] == rank.solved.inconsistency[[x]]))


solved.inconsistencies.per.type <- rbind(solved.inconsistencies.WR, solved.inconsistencies.WC, solved.inconsistencies.H)
colnames(solved.inconsistencies.per.type) <- rank.names
solved.inconsistencies.per.type

# summed where inconsistencies are solved per type
sum.solved.ranks.WR <- sapply(1:length(count.rank.solved.inconsistency), function(x)
  sum(solved.inconsistencies.per.type[1,1:x]))
sum.solved.ranks.WC <- sapply(1:length(count.rank.solved.inconsistency), function(x)
  sum(solved.inconsistencies.per.type[2,1:x]))
sum.solved.ranks.H <- sapply(1:length(count.rank.solved.inconsistency), function(x)
  sum(solved.inconsistencies.per.type[3,1:x]))
  
sum.solved.ranks.per.type <- rbind(sum.solved.ranks.WR,sum.solved.ranks.WC,sum.solved.ranks.H)
colnames(sum.solved.ranks.per.type) <- rank.names
sum.solved.ranks.per.type# <- sum.solved.ranks.per.type[,2:ncol(sum.solved.ranks.per.type)]

normalize.sum.solved.ranks.per.type <- sum.solved.ranks.per.type*100/c(sum.solved.ranks.per.type[,ncol(sum.solved.ranks.per.type)])  #as.matrix(sum.solved.ranks.per.type)

melt.sum.solved.ranks.per.type <-melt(normalize.sum.solved.ranks.per.type,id.vars=names(normalize.sum.solved.ranks.per.type))


ggplot(melt.sum.solved.ranks.per.type,aes(x=factor(Var2),y=value,fill=factor(Var1)))+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_brewer(palette = "Set1", name="Inconsistency Type",
                    labels=c("Weak Reads", "Weak Contigs", "Hard"))+
  xlab("Taxonomic Rank")+
  ylab("Percentage of Inconsistencies Solved")+
  ggtitle("Taxonomic rank where Inconsistencies are Solved (SSD)", subtitle = "Percentage for each Type")+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=14),
        plot.title =element_text(size=25,face="bold"),
        plot.subtitle = element_text(size=18),
        legend.text=element_text(size=15),
        legend.title =element_text(size=18))
###################################
#-------------------
read.taxon.ranks <- sort(unique(solved.data$ReadTopRank), decreasing = TRUE)
contig.taxon.ranks <- sort(unique(solved.data$ContigTopRank), decreasing = TRUE)

read.rank.counter <-  sapply(1:length(read.taxon.ranks), function(x) 
  sum(solved.data$ReadTopRank == read.taxon.ranks[[x]]) )
contig.rank.counter <-  sapply(1:length(contig.taxon.ranks), function(x) 
  sum(solved.data$ContigTopRank == contig.taxon.ranks[[x]]) )

barplot(rbind(read.rank.counter,contig.rank.counter), beside=TRUE, log = "y", main = "Taxonomic Rank where Inconsistencies are Solved")

