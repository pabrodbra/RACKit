# Made by Pablo Rodríguez Brazzarola
# RACC Results Processing Functions
# Libraries
library(ggplot2)
library(reshape2)

# Constants
rank.names <- c("Species","Genus","Family","Order","Class","Phylum","Kingdom","Domain", "Organism")
inconsistency.types <- c("WR","WC","H")
plot.width <- 1024
plot.height <- 

#1 BEAR pa file -> orig.distro.per.taxa
#2 MEGAN Reads Distro per Taxa -> reads.distro.per.taxa
#3 MEGAN Contigs Distro per Taxa -> contigs.distro.per.taxa

### Original Distro
# I: 1
# O: Plot

### Distro Comparison
# I: 1,2,3
# O: Plot

### RMSE
# I: 1,2,3
# O: Table CSV

### Inconsistencies Found
# I: 
# O: Table CSV

### Inconsistency Resolution
# I: Solved_Inconsistencies_CSV
# O: Plot
inconsistency_resolution <- function(solved_inconsistencies_csv, output.path){
    # Load data
    solved.data <-read.csv2(solved_inconsistencies_csv, header = TRUE, sep = ',')

    # Sort by RankSolved
    rank.solved.inconsistency <- sort(unique(solved.data$RankSolved), decreasing=TRUE)

    # Where are Inconsistencies Solved
    count.rank.solved.inconsistency <- sapply(1:length(rank.solved.inconsistency), function(x)
        sum(solved.data$RankSolved == rank.solved.inconsistency[[x]]))
    names(count.rank.solved.inconsistency) <- rank.names

    # Summed where are Inconsistencies Solved
    sum.solved.ranks <- sapply(1:length(count.rank.solved.inconsistency), function(x)
        sum(count.rank.solved.inconsistency[1:x]))
    names(sum.solved.ranks) <- rank.names

    # Type of Inconsistencies
    count.types <- sapply(1:length(inconsistency.types), function(x)
    sum(solved.data$Type == inconsistency.types[[x]]))
    names(count.types) <- inconsistency.types

    # Where are Inconsistencies Solved for each Inconsistency Type
    solved.inconsistencies.WR <- sapply(1:length(rank.solved.inconsistency), function(x)
        sum(solved.data$RankSolved[which(solved.data$Type == "WR")] == rank.solved.inconsistency[[x]]))
    solved.inconsistencies.WC <- sapply(1:length(rank.solved.inconsistency), function(x)
        sum(solved.data$RankSolved[which(solved.data$Type == "WC")] == rank.solved.inconsistency[[x]]))
    solved.inconsistencies.H <- sapply(1:length(rank.solved.inconsistency), function(x)
        sum(solved.data$RankSolved[which(solved.data$Type == "H")] == rank.solved.inconsistency[[x]]))

    solved.inconsistencies.per.type <- rbind(solved.inconsistencies.WR, solved.inconsistencies.WC, solved.inconsistencies.H)
    colnames(solved.inconsistencies.per.type) <- rank.names

    # Sum where Inconsistencies are Solved per Inconsistency Type
    sum.solved.ranks.WR <- sapply(1:length(count.rank.solved.inconsistency), function(x)
    sum(solved.inconsistencies.per.type[1,1:x]))
    sum.solved.ranks.WC <- sapply(1:length(count.rank.solved.inconsistency), function(x)
    sum(solved.inconsistencies.per.type[2,1:x]))
    sum.solved.ranks.H <- sapply(1:length(count.rank.solved.inconsistency), function(x)
    sum(solved.inconsistencies.per.type[3,1:x]))

    sum.solved.ranks.per.type <- rbind(sum.solved.ranks.WR,sum.solved.ranks.WC,sum.solved.ranks.H)
    colnames(sum.solved.ranks.per.type) <- rank.names

    # Normalize
    normalize.sum.solved.ranks.per.type <- sum.solved.ranks.per.type*100/c(sum.solved.ranks.per.type[,ncol(sum.solved.ranks.per.type)])  #as.matrix(sum.solved.ranks.per.type)
    # Melt
    melt.sum.solved.ranks.per.type <-melt(normalize.sum.solved.ranks.per.type,id.vars=names(normalize.sum.solved.ranks.per.type))
    # Plot
    png(filename = output.path, width = plot.width, height = plot.height)
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
    dev.off()

    output.return <- list(
        "count.rank.solved.inconsistency" = count.rank.solved.inconsistency,
        "sum.solved.ranks" = sum.solved.ranks,
        "count.types" = count.types,
        "solved.inconsistencies.per.type" = solved.inconsistencies.per.type,
        "sum.solved.ranks.per.type" = sum.solved.ranks.per.type,
        "normalize.sum.solved.ranks.per.type" = normalize.sum.solved.ranks.per.type,
    )
    return(output.return)
}


### DB Coverage Comparison
# I: 
# O: Table CSV

### Statistical Measurements (TP/FP/TN/FN)
# I:
# O: Table CSV