# Made by Pablo Rodríguez Brazzarola
# RACC Results Processing Functions
# Libraries
list.of.packages <- c("ggplot2", "reshape2", "knitr", "rmarkdown")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

library(ggplot2)
library(reshape2)
library(knitr)

# Constants
rank.names <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
inconsistency.types <- c("WR","WC","H")
plot.width <- 1024
plot.height <- 1024
plot.width.in <- 12
plot.height.in <- 12
CEX.MAIN <- 2
CEX.LAB <- 2.5
CEX.AXIS <- 2
CEX.LEGEND <- 1.6
GG.AXIS.TITLE <- 22
GG.AXIS.TEXT <- 18
GG.PLOT.MAIN <- 25
GG.PLOT.SUB <- 18
GG.LEGEND.TEXT <- 15
GG.LEGEND.TITLE <- 18

#1 GRINDER Rank -> orig.distro.per.taxa
#2 MEGAN Reads Distro per Taxa -> reads.distro.per.taxa
#3 MEGAN Contigs Distro per Taxa -> contigs.distro.per.taxa
### Helpers
is.nan.data.frame <- function(x){
  do.call(cbind, lapply(x, is.nan))
}

list.to.dataframe <- function(l, set.column.name, set.row.names = ""){
  df <- data.frame(matrix(unlist(l), nrow=length(l), byrow=T))
  colnames(df) <- set.column.name
  ifelse(set.row.names=="",rownames(df) <- names(l), rownames(df) <- set.row.names)
  
  return(df)
}

join_list_comparisons <- function(vector.list.comparisons){
  ret <- setNames(data.frame(matrix(ncol = length(names(vector.list.comparisons)), nrow = 0)), 
                  names(vector.list.comparisons))
  for (list.comp in vector.list.comparisons){
    ret <- rbind(ret, list.comp, stringsAsFactors = FALSE)
  }
  return(ret)
}

### Prepare descendant original, reads and contigs distro data + absolute difference relative to original
# I: original.distro.file, read.distro.file, contig.distro.file
# O: Normalized data
preprocess_distribution_data <- function(original.distro.path,
                                         reads.distro.path,
                                         contigs.distro.path){
  original.distro.per.specie <- read.csv2(original.distro.path, header=TRUE, sep = '\t')
  reads.distro.per.specie <- read.csv2(reads.distro.path, header=FALSE, sep = ',')
  contigs.distro.per.specie <- read.csv2(contigs.distro.path, header=FALSE, sep = ',')
  
  original.distro.per.specie[,2] <- as.numeric(original.distro.per.specie[,2])
  # reads.distro.per.specie[,2] <- as.numeric(reads.distro.per.specie[,2])
  # contigs.distro.per.specie[,2] <- as.numeric(contigs.distro.per.specie[,2])
  
  # Make all species be in all sets
  which.reads <- which(!(reads.distro.per.specie$V1 %in% contigs.distro.per.specie$V1)==TRUE)
  which.contigs <- which(!(contigs.distro.per.specie$V1 %in% reads.distro.per.specie$V1)==TRUE)
  which.reads.info <- reads.distro.per.specie[which.reads,]
  which.contigs.info <- contigs.distro.per.specie[which.contigs,]
  which.reads.null <- cbind(as.vector(which.contigs.info$V1), rep(0, length(which.contigs.info$V1)))
  which.contigs.null <- cbind(as.vector(which.reads.info$V1), rep(0, length(which.reads.info$V1)))
  
  if( !identical(which.reads, integer(0)) )
    reads.distro.per.specie <- rbind(reads.distro.per.specie[-which.reads,],which.reads.info,which.reads.null)
  if( !identical(which.contigs, integer(0)) )
    contigs.distro.per.specie <- rbind(contigs.distro.per.specie[-which.contigs,], which.contigs.null, which.contigs.info)
  
  # Original distro descendant
  orig.desc <- c()
  for(specie in reads.distro.per.specie[,1]){
    orig.desc <- c(orig.desc, sum(original.distro.per.specie[grep(specie, original.distro.per.specie[,1]),2])) 
  }
  
  # Prepare Values
  orig.desc <- as.numeric(orig.desc)/sum(orig.desc)
  reads.desc <- as.numeric( as.numeric(reads.distro.per.specie[,2]) )*100/sum( as.numeric(reads.distro.per.specie[,2]) )
  contigs.desc <- as.numeric( as.numeric(contigs.distro.per.specie[,2]) )*100/sum( as.numeric(contigs.distro.per.specie[,2]) )
  names(orig.desc) <- names(reads.desc) <- names(contigs.desc) <- reads.distro.per.specie[,1]
  
  # Distribution info
  distro.info <- rbind(reads.desc, orig.desc, contigs.desc)
  distro.info.2 <- rbind(abs(orig.desc-reads.desc), abs(orig.desc-contigs.desc))
  
  ret <- list(
    "reads.descendant.distribution" = distro.info[1,],
    "orig.descendant.distribution" = distro.info[2,],
    "contigs.descendant.distribution" = distro.info[3,],
    "reads.absolute.difference" = distro.info.2[1,],
    "contigs.absolute.difference" = distro.info.2[2,]
  )
  return (ret)
}

### Original Distro
# I: original.distro.file
# O: Plot
original_distribution_plot <- function(orig.desc,
                                       output.path="IMG-original_ditribution_plot.png",
                                       PNG=TRUE, title="", max.species=50){
  # Distribution lines
  sort.data <- data.frame(orig=orig.desc)
  sort.data <- orig.desc[order(sort.data$orig, decreasing = TRUE)]
  sort.data <- sort.data[sort.data!=0]
  if(length(orig.desc) > max.species)
    sort.data <- sort.data[seq_len(max.species)]
  
  
  tmp.df <- data.frame(value = sort.data, specie = names(sort.data), group = rep("orig", length(sort.data)))
  dat <- melt(tmp.df, value.name = "value", varnames = c("specie", "group"))
  dat$group <- factor(dat$group, levels = "orig")
  dat$specie = factor(dat$specie, levels = dat$specie[order(-dat$value[dat$group=="orig"])])
  
  # Plot
  if(PNG)
    png(filename = output.path, width = plot.width, height = plot.height)
  t <- sapply(1:length(names(orig.desc)), function(x) "")
  par(mar=c(6,6,5,1.5))
  bks <-c(0.1,0.5,1,5,10,25,50,100)
  
  p <- ggplot(dat, aes(x = specie, y = value, group = group)) + 
    geom_line(aes(color=group), size=1.5) + #expand_limits(y=0) + 
    scale_y_continuous(labels = function (x) paste0(x,"%"), breaks = bks) + 
    coord_trans(y="log10", limy = c(0.1,100)) + 
    scale_colour_brewer(palette = "Set1", name="", labels = "Original")+
    xlab("Species")+
    ylab("Percentage of sequences")+
    ggtitle("Relative abundance per Specie", subtitle = title)+
    theme(axis.title = element_text(size=GG.AXIS.TITLE),
          axis.text = element_text(size=GG.AXIS.TEXT),
          axis.text.x = element_blank(),
          plot.title =element_text(size=GG.PLOT.MAIN,face="bold"),
          plot.subtitle = element_text(size=GG.PLOT.SUB),
          legend.text=element_text(size=GG.LEGEND.TEXT),
          legend.title =element_text(size=GG.LEGEND.TITLE))
  plot(p)
  if(PNG)
    dev.off()
}

### Distro Comparison
# I: original.distro.file, read.distro.file, contig.distro.file
# O: Plot
distribution_comparison_plot <- function(orig.desc,
                                         reads.desc,
                                         contigs.desc,
                                         output.path="IMG-compared_ditribution_plot.png",
                                         PNG=TRUE, title = "", max.species=50){
  # Distribution lines
  sort.data <- data.frame(orig=orig.desc*100,reads=reads.desc,contig=contigs.desc)
  sort.data <- sort.data[order(sort.data$orig, decreasing = TRUE),]
  row.sub <- apply(sort.data, 1, function(x) all(x!=0))
  sort.data[row.sub,]
  
  if(length(orig.desc) > max.species)
    sort.data <- sort.data[seq_len(max.species),]
  
  tmp.df <- data.frame(value = c(sort.data$orig, sort.data$reads, sort.data$contig),
                       specie = rep(rownames(sort.data), 3),
                       group = c( rep("orig", nrow(sort.data)), rep("reads", nrow(sort.data)), rep("contig", nrow(sort.data)) ))
  dat <- melt(tmp.df, value.name = "value", varnames = c("specie", "group"))
  dat$group = factor(dat$group, levels=c("orig","reads","contig"))
  dat$specie = factor(dat$specie, levels = dat$specie[order(-dat$value[dat$group=="orig"])])
  
  # Plot
  if(PNG)
    png(filename = output.path, width = plot.width, height = plot.height)
  t <- sapply(1:length(rownames(sort.data)), function(x) "")
  par(mar=c(6,6,5,1.5))
  bks <-c(0.1,0.5,1,5,10,25,50,100)
  
  p <- ggplot(dat, aes(x = specie, y = value, group = group)) + 
    geom_line(aes(color=group), size=1.5) + 
    scale_y_continuous(labels = function (x) paste0(x,"%"), breaks = bks) + 
    coord_trans(y="log10", limy = c(0.1,100)) + 
    scale_colour_brewer(palette = "Set1", name="", labels = c("Original", "Reads", "Contigs"))+#labs(color="")+
    xlab("Species")+
    ylab("Percentage of sequences")+
    ggtitle("Relative abundance per Specie", subtitle = title)+
    theme(axis.title = element_text(size=GG.AXIS.TITLE),
          axis.text = element_text(size=GG.AXIS.TEXT),
          axis.text.x = element_blank(),
          plot.title =element_text(size=GG.PLOT.MAIN,face="bold"),
          plot.subtitle = element_text(size=GG.PLOT.SUB),
          legend.text=element_text(size=GG.LEGEND.TEXT),
          legend.title =element_text(size=GG.LEGEND.TITLE))
  plot(p)
  if(PNG)
    dev.off()
}

### RMSE
# I: 1,2,3
# O: Table CSV
rmse_comparison <- function(orig.desc,
                            reads.desc,
                            contigs.desc){
  reads.distro.rmse <- (sum((orig.desc-reads.desc)^2)/length(orig.desc) )^(1/2)
  contigs.distro.rmse <- (sum((orig.desc-contigs.desc)^2)/length(orig.desc) )^(1/2)
  
  ret <- list(
    "reads.rmse" = reads.distro.rmse,
    "contigs.rmse" = contigs.distro.rmse
  )
  return(ret)
}

### Inconsistency Resolution
# I: Solved_Inconsistencies_CSV
# O: Plot
inconsistency_resolution <- function(solved_inconsistencies_csv, 
                                     output.path="IMG-inconsistency_solver.png",
                                     img.title = "Taxonomic rank where Inconsistencies are Solved",
                                     PNG=TRUE){
    # Load data
    solved.data <-read.csv2(solved_inconsistencies_csv, header = TRUE, sep = ',')

    # Sort by RankSolved
    rank.solved.inconsistency <- sort(unique(solved.data$RankSolved), decreasing=TRUE)

    # Where are Inconsistencies Solved
    count.rank.solved.inconsistency <- sapply(1:length(rank.solved.inconsistency), function(x)
        sum(solved.data$RankSolved == rank.solved.inconsistency[[x]]))
    names(count.rank.solved.inconsistency) <- rank.names[rank.solved.inconsistency]

    # Summed where are Inconsistencies Solved
    sum.solved.ranks <- sapply(1:length(count.rank.solved.inconsistency), function(x)
        sum(count.rank.solved.inconsistency[1:x]))
    names(sum.solved.ranks) <- rank.names[rank.solved.inconsistency]

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
    
    if(ncol(solved.inconsistencies.per.type) == 8)
      colnames(solved.inconsistencies.per.type) <- c(rank.names[rank.solved.inconsistency], "Organism")
    else
      colnames(solved.inconsistencies.per.type) <- rank.names[rank.solved.inconsistency]
    
    # Sum where Inconsistencies are Solved per Inconsistency Type
    sum.solved.ranks.WR <- sapply(1:length(count.rank.solved.inconsistency), function(x)
    sum(solved.inconsistencies.per.type[1,1:x]))
    sum.solved.ranks.WC <- sapply(1:length(count.rank.solved.inconsistency), function(x)
    sum(solved.inconsistencies.per.type[2,1:x]))
    sum.solved.ranks.H <- sapply(1:length(count.rank.solved.inconsistency), function(x)
    sum(solved.inconsistencies.per.type[3,1:x]))

    sum.solved.ranks.per.type <- rbind(sum.solved.ranks.WR,sum.solved.ranks.WC,sum.solved.ranks.H)
    
    if(ncol(sum.solved.ranks.per.type) == 8)
      colnames(sum.solved.ranks.per.type) <- c(rank.names[rank.solved.inconsistency], "Organism")
    else
      colnames(sum.solved.ranks.per.type) <- rank.names[rank.solved.inconsistency]

    # Normalize
    normalize.sum.solved.ranks.per.type <- sum.solved.ranks.per.type*100/c(sum.solved.ranks.per.type[,ncol(sum.solved.ranks.per.type)])  #as.matrix(sum.solved.ranks.per.type)
    normalize.sum.solved.ranks.per.type[is.nan.data.frame(normalize.sum.solved.ranks.per.type)] <- 0
    # Melt
    melt.sum.solved.ranks.per.type <-melt(normalize.sum.solved.ranks.per.type,id.vars=names(normalize.sum.solved.ranks.per.type))
    # Plot
    #png(filename = output.path, width = plot.width, height = plot.height)
    p <- ggplot(melt.sum.solved.ranks.per.type,aes(x=factor(Var2),y=value,fill=factor(Var1)))+
        geom_bar(stat="identity",position="dodge")+
        scale_fill_brewer(palette = "Set1", name="Inconsistency Type",
                        labels=c("Weak Reads", "Weak Contigs", "Hard"))+
        xlab("Taxonomic Rank")+
        ylab("Percentage of Inconsistencies Solved")+
        ggtitle(img.title, subtitle = "Percentage for each Type")+
        theme(axis.title = element_text(size=GG.AXIS.TITLE),
            axis.text = element_text(size=GG.AXIS.TEXT),
            plot.title =element_text(size=GG.PLOT.MAIN,face="bold"),
            plot.subtitle = element_text(size=GG.PLOT.SUB),
            legend.text=element_text(size=GG.LEGEND.TEXT),
            legend.title =element_text(size=GG.LEGEND.TITLE))
    plot(p)
    if(PNG)
      ggsave(filename = output.path, width = plot.height.in, height = plot.height.in, units = 'in')

    output.return <- list(
        "count.rank.solved.inconsistency" = count.rank.solved.inconsistency,
        "sum.solved.ranks" = sum.solved.ranks,
        "count.types" = count.types,
        "solved.inconsistencies.per.type" = solved.inconsistencies.per.type,
        "sum.solved.ranks.per.type" = sum.solved.ranks.per.type,
        "normalize.sum.solved.ranks.per.type" = normalize.sum.solved.ranks.per.type
    )
    return(output.return)
}


### DB Coverage Comparison
# I: coverage.info
# O: Table CSV
coverage_comparison <- function(coverage.info.path){
  # Load data
  coverage.info <- read.csv2(coverage.info.path, header=TRUE, sep = '\t')
  #coverage.info <- coverage.info[c("ReadHits", "ContigHits", "BothYes",  "BothNo",)]
  coverage.total <- sum(as.numeric(coverage.info))
  
  coverage.percentage <- c(signif(coverage.info[1]*100/coverage.total), signif(coverage.info[2]*100/coverage.total)
                      , signif(coverage.info[3]*100/coverage.total), signif(coverage.info[4]*100/coverage.total))
  
  ret <- list(
    "read.only.cov" = coverage.percentage$ReadHits,
    "contig.only.cov" = coverage.percentage$ContigHits,
    "read.total.cov" = coverage.percentage$ReadHits + coverage.percentage$BothYes,
    "contig.total.cov" = coverage.percentage$ContigHits + coverage.percentage$BothYes,
    "shared.cov" = coverage.percentage$BothYes,
    "not.cov" = coverage.percentage$BothNo
  )
  
  return (ret)
}

#--------------TODO

### Statistical Measurements (TP/FP/TN/FN)
# I: rac-statistical_measurements.csv
# O: Table CSV


### Inconsistencies Found
# I: inc-solver.out
# O: Table CSV


### Multi RACKIT Results analysis