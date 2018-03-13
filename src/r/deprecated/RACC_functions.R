# RACC functions
# By: Pablo Rodriguez Brazzarola

RACC.uniseqDB.coverage.plot <- function(coverage.data, output.name="DB_coverage.png"){
  coverage.info <- read.csv(coverage.data, header=TRUE, sep = '\t')
  coverage.info <- coverage.info[c("ReadOnly", "BothYes", "ContigOnly", "BothNo")]
  coverage.total <- sum(as.numeric(coverage.info))
  
  pie.colors <- c("dodgerblue4", "dodgerblue1","deepskyblue", "olivedrab2")
  pie.names <- c("Read Only", "Both", "Contig Only", "None")
  
  coverage.info["ReadOnly"] <- coverage.info["ReadOnly"]-coverage.info["BothYes"]
  coverage.info["ContigOnly"] <- coverage.info["ContigOnly"]-coverage.info["BothYes"]
  
  coverage.label <- c(signif(coverage.info[1]*100/coverage.total), signif(coverage.info[2]*100/coverage.total)
                      , signif(coverage.info[3]*100/coverage.total), signif(coverage.info[4]*100/coverage.total))
  
  png(filename = output.name, width = 1024, height = 1024)
  pie(as.numeric(coverage.info), labels = paste("Percentage: ", coverage.label, "\nNucleotides (Mbp): ", signif(coverage.info/1000000)), col = pie.colors
      , main = paste("Percentage of Database Covered (Total Mbp = ", signif(sum(coverage.info/1000000)), ")", sep=""), init.angle = 45)
  legend("bottomleft", pie.names, cex = 1.5, fill = pie.colors)
  dev.off()
}

RACC.inconsistency.full <- function(inconsistency.data, n.reads, n.contigs, output.name="INCplot", filter.flag=FALSE, binary.flag=TRUE, null.tax.flag = FALSE){
  inconsistency_info <- read.csv(inconsistency.data, header=TRUE)
  
  if(null.tax.flag == TRUE){
    inconsistency_info <- inconsistency_info[-which(inconsistency_info$ReadTaxonomy=="None"),]
    inconsistency_info <- inconsistency_info[-which(inconsistency_info$ContigTaxonomy=="None"),]
  }
  
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
  
  RACC.inconsistency.barplots(read.inconsistencies = read.inconsistencies, contig.inconsistencies = contig.inconsistencies, output.name = output.name, filter=filter.flag)
  RACC.inconsistency.barplots_v2(inconsistency_info = inconsistency_info, output.name = output.name, filter = filter.flag)
  RACC.inconsistency.pieplots_v2(inconsistency_info = inconsistency_info, n.reads.total = n.reads, n.contigs.total = n.contigs, output.name = output.name, binary=binary.flag)
}

RACC.inconsistency.barplots <- function(read.inconsistencies, contig.inconsistencies, output.name="INCplot", filter=FALSE){
  colorfunc <- colorRampPalette(c("light green", "dark green", "light blue", "dark blue", "purple", "dark red", "orange"))
  
  if(filter == TRUE){
    filtered.read.inconsistencies <- read.inconsistencies[which(read.inconsistencies > sum(read.inconsistencies)/100)]
    filtered.contigs.inconsistencies <-  contig.inconsistencies[which(contig.inconsistencies > sum(contig.inconsistencies)/100)]
    
    png(filename=paste(output.name,"_filtered_reads_inconsistencies.png", sep = ""), width = 1024, height = 1024)
    barplot(filtered.read.inconsistencies, main = "Reads Inconsistencies (>1%) per Taxonomy", ylab = "Inconsistencies", xlab = "Taxonomy", 
            col = colorfunc(length(filtered.read.inconsistencies)), beside = TRUE, cex.names = 0.6, log="y")
    legend("topright", names(filtered.read.inconsistencies), cex = 1.5, fill = colorfunc(length(filtered.read.inconsistencies)), ncol = 2)
    dev.off()
    
    png(filename=paste(output.name,"_filtered_contigs_inconsistencies.png", sep = ""), width = 1024, height = 1024)
    barplot(filtered.contigs.inconsistencies, main = "Contigs Inconsistencies (>1%) per Taxonomy", ylab = "Inconsistencies", xlab = "Taxonomy", 
            col = colorfunc(length(filtered.contigs.inconsistencies)), beside = TRUE, cex.names = 0.6, log="y")
    legend("topright", names(filtered.contigs.inconsistencies), cex = 1.5, fill = colorfunc(length(filtered.contigs.inconsistencies)), ncol = 2)
    dev.off()
    
  }else{
    png(filename=paste(output.name,"_reads_inconsistencies.png", sep = ""), width = 1024, height = 1024)
    barplot(read.inconsistencies, main = "Reads Inconsistencies per Taxonomy", ylab = "Inconsistencies", xlab = "Taxonomy", 
            col = colorfunc(length(read.inconsistencies)), beside = TRUE, cex.names = 0.6, log="y")
    legend("topright", names(read.inconsistencies)[1:15], cex = 1.5, fill = colorfunc(length(read.inconsistencies)), ncol = 2)
    dev.off()
    
    png(filename=paste(output.name,"_contigs_inconsistencies.png", sep = ""), width = 1024, height = 1024)
    barplot(contig.inconsistencies, main = "Contigs Inconsistencies per Taxonomy", ylab = "Inconsistencies", xlab = "Taxonomy", 
            col = colorfunc(length(contig.inconsistencies)), beside = TRUE, cex.names = 0.6, log="y")
    legend("topright", names(contig.inconsistencies)[1:15], cex = 1.5, fill = colorfunc(length(contig.inconsistencies)), ncol = 2)
    dev.off()
  }
  
}

RACC.inconsistency.barplots_v2 <- function(inconsistency_info, output.name="INCplot", filter=FALSE){
  colorfunc <- colorRampPalette(c("light green", "dark green", "light blue", "dark blue", "purple", "dark red", "orange"))
  
  # Get taxons
  inconsistency_info <- inconsistency_info[-which(inconsistency_info$ReadTaxonomy == "None"),]
  inconsistency_info <- inconsistency_info[-which(inconsistency_info$ContigTaxonomy == "None"),]
  
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
  
  if(filter==FALSE){
    #Draw Images
    png(filename=paste(output.name,"_HARD_reads_inconsistencies.png", sep = ""), width = 1024, height = 1024)
    barplot(read.inconsistencies, main = "Reads Strong Inconsistencies per Taxonomy", ylab = "Inconsistencies", xlab = "Taxonomy", 
            col = colorfunc(length(read.inconsistencies)), beside = TRUE, cex.names = 0.6, log="y")
    legend("topright", names(read.inconsistencies)[1:15], cex = 1.5, fill = colorfunc(length(read.inconsistencies)), ncol = 2)
    dev.off()
    
    png(filename=paste(output.name,"_HARD_contigs_inconsistencies.png", sep = ""), width = 1024, height = 1024)
    barplot(contig.inconsistencies, main = "Contigs Strong Inconsistencies per Taxonomy", ylab = "Inconsistencies", xlab = "Taxonomy", 
            col = colorfunc(length(contig.inconsistencies)), beside = TRUE, cex.names = 0.6, log="y")
    legend("topright", names(contig.inconsistencies)[1:15], cex = 1.5, fill = colorfunc(length(contig.inconsistencies)), ncol = 2)
    dev.off()
  }else{
    filtered.read.inconsistencies <- read.inconsistencies[which(read.inconsistencies > sum(read.inconsistencies)/100)]
    filtered.contigs.inconsistencies <-  contig.inconsistencies[which(contig.inconsistencies > sum(contig.inconsistencies)/100)]
  
    #Draw Images
    png(filename=paste(output.name,"_HARD_filtered_reads_inconsistencies.png", sep = ""), width = 1024, height = 1024)
    barplot(filtered.read.inconsistencies, main = "Reads Strong Inconsistencies (>1%) per Taxonomy", ylab = "Inconsistencies", xlab = "Taxonomy", 
            col = colorfunc(length(filtered.read.inconsistencies)), beside = TRUE, cex.names = 0.6, log="y")
    legend("topright", names(filtered.read.inconsistencies), cex = 1.5, fill = colorfunc(length(filtered.read.inconsistencies)), ncol = 2)
    dev.off()
    
    png(filename=paste(output.name,"_HARD_filtered_contigs_inconsistencies.png", sep = ""), width = 1024, height = 1024)
    barplot(filtered.contigs.inconsistencies, main = "Contigs Strong Inconsistencies (>1%) per Taxonomy", ylab = "Inconsistencies", xlab = "Taxonomy", 
            col = colorfunc(length(filtered.contigs.inconsistencies)), beside = TRUE, cex.names = 0.6, log="y")
    legend("topright", names(filtered.contigs.inconsistencies)[1:15], cex = 1.5, fill = colorfunc(length(filtered.contigs.inconsistencies)), ncol = 2)
    dev.off()
  }
}

RACC.inconsistency.pieplots <- function(inconsistency_info, n.reads.total, n.contigs.total, output.name="INCplot", binary=TRUE){
  if(binary == TRUE){
    tot.colors <- c("blue", "orange")
    tot.names <- c("Consistent", "Inconsistent")
    
    n.reads.inconsistent <- length(unique(inconsistency_info$ReadID))
    reads.label = c(signif((n.reads.total-n.reads.inconsistent)*100/n.reads.total), signif(n.reads.inconsistent*100/n.reads.total))
    
    png(filename=paste(output.name,"_total_reads.png", sep = ""), width = 1024, height = 1024)
    pie(c(n.reads.total-n.reads.inconsistent, n.reads.inconsistent), labels = paste("Percentage: ",reads.label, "\nAmount: ", c(n.reads.total-n.reads.inconsistent, n.reads.inconsistent)), col = tot.colors
        , main = "Reads Taxonomy Status after classification")
    legend("topleft", tot.names, cex = 1.5, fill = tot.colors)
    dev.off()
    
    # Contigs Assigned
    n.contigs.inconsistent <- length(unique(inconsistency_info$ContigID))
    contigs.label <-  c(signif((n.contigs.total-n.contigs.inconsistent)*100/n.contigs.total), signif(n.contigs.inconsistent*100/n.contigs.total))
    
    png(filename=paste(output.name,"_total_contigs.png", sep = ""), width = 1024, height = 1024)
    pie(c(n.contigs.total-n.contigs.inconsistent, n.contigs.inconsistent), labels = paste("Percentage: ", contigs.label, "\nAmount: ", c(n.contigs.total-n.contigs.inconsistent, n.contigs.inconsistent)), col = tot.colors
        , main = "Contigs Taxonomy Status after classification")
    legend("topleft", tot.names, cex = 1.5, fill = tot.colors)
    dev.off()
  }
  else{
    # Obtain Semiconsistent (-) and Inconsistent (+)
    tot.colors <- c("blue", "orange", "red")
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
    n.reads.consistant <- n.reads.total-n.reads.inconsistent.plus-n.reads.inconsistent.minus
    reads.label = c(signif(n.reads.consistant*100/n.reads.total), signif(n.reads.inconsistent.plus*100/n.reads.total)
                    , signif(n.reads.inconsistent.minus*100/n.reads.total))
    
    png(filename=paste(output.name,"_nb_total_reads.png"), width = 1024, height = 1024)
    pie(c(n.reads.consistant, n.reads.inconsistent.plus, n.reads.inconsistent.minus), labels = paste("Percentage: ",reads.label, 
                                                                                                     "\nAmount: ", c(n.reads.consistant, n.reads.inconsistent.plus, n.reads.inconsistent.minus)), col = tot.colors
        , main = "Reads Taxonomy Status after classification")
    legend("topleft", tot.names, cex = 1.5, fill = tot.colors)
    dev.off()
    
    # Contigs
    n.contigs.inconsistent <- length(unique(inconsistency_info$ContigID))
    n.contigs.consistant <- n.contigs.total-n.contigs.inconsistent.plus-n.contigs.inconsistent.minus
    contigs.label <-  c(signif(n.contigs.consistant*100/n.contigs.total), signif(n.contigs.inconsistent.plus*100/n.contigs.total)
                        , signif(n.contigs.inconsistent.minus*100/n.reads.total))
    
    png(filename=paste(output.name,"_nb_total_contigs.png"), width = 1024, height = 1024)
    pie(c(n.contigs.consistant, n.contigs.inconsistent.plus, n.contigs.inconsistent.minus), labels = paste("Percentage: ", contigs.label,
                                                                                                           "\nAmount: ", c(n.contigs.consistant, n.contigs.inconsistent.plus, n.contigs.inconsistent.minus)), col = tot.colors, 
        main = "Contigs Taxonomy Status after classification")
    legend("topleft", tot.names, cex = 1.5, fill = tot.colors)
    dev.off()
  }
}


RACC.inconsistency.pieplots_v2 <- function(inconsistency_info, n.reads.total, n.contigs.total, output.name="INCplot", binary=TRUE){
  if(binary == TRUE){
    tot.colors <- c("blue", "chocolate1", "chocolate4", "coral3")
    tot.names <- c("Consistent", "Weakly Inconsistent Reads", "Strongly Inconsistent", "Weakly Inconsistent Contigs")
    
    n.reads.inconsistent <- length(unique(inconsistency_info$ReadID))
    n.none.r.reads <- length(unique(inconsistency_info[which(inconsistency_info$ReadTaxonomy == "None"),]$ReadID))
    n.none.r.contigs <- length(unique(inconsistency_info[which(inconsistency_info$ContigTaxonomy == "None"),]$ReadID))
    n.none.r.total <- n.none.r.reads + n.none.r.contigs
    
    reads.label = c(signif((n.reads.total-n.reads.inconsistent)*100/n.reads.total), signif(n.none.r.reads*100/n.reads.total)
                                                                                           , signif((n.reads.inconsistent-n.none.r.total)*100/n.reads.total)
                                                                                                    , signif(n.none.r.contigs*100/n.reads.total))
    reads.data <- c((n.reads.total-n.reads.inconsistent),n.none.r.reads,(n.reads.inconsistent-n.none.r.total), n.none.r.contigs)
    
    png(filename=paste(output.name,"_total_reads.png", sep = ""), width = 1024, height = 1024)
    pie(reads.data, labels = paste("Percentage: ",reads.label, "\nAmount: ", reads.data), col = tot.colors
        , main = "Reads Taxonomy Status after classification")
    legend("topleft", paste(tot.names, " : ", reads.data, sep=""), cex = 1.5, fill = tot.colors)
    dev.off()
    
    # Contigs Assigned
    n.contigs.inconsistent <- length(unique(inconsistency_info$ContigID))
    n.none.c.reads <- length(unique(inconsistency_info[which(inconsistency_info$ReadTaxonomy == "None"),]$ContigID))
    n.none.c.contigs <- length(unique(inconsistency_info[which(inconsistency_info$ContigTaxonomy == "None"),]$ContigID))
    n.none.c.total <- n.none.c.reads + n.none.c.contigs
    
    contigs.label = c(signif((n.contigs.total-n.contigs.inconsistent)*100/n.contigs.total), signif(n.none.c.reads*100/n.contigs.total)
                    , signif((n.contigs.inconsistent-n.none.c.total)*100/n.contigs.total)
                    , signif(n.none.c.contigs*100/n.contigs.total))
    contigs.data <- c((n.contigs.total-n.contigs.inconsistent),n.none.c.reads,(n.contigs.inconsistent-n.none.c.total), n.none.c.contigs)
    
    png(filename=paste(output.name,"_total_contigs.png", sep = ""), width = 1024, height = 1024)
    pie(contigs.data, labels = paste("Percentage: ", contigs.label, "\nAmount: ", contigs.data), col = tot.colors
        , main = "Contigs Taxonomy Status after classification")
    legend("topleft", paste(tot.names, " : ", contigs.data, sep=""), cex = 1.5, fill = tot.colors)
    dev.off()
  }
  else{
    # Obtain Semiconsistent (-) and Inconsistent (+)
    tot.colors <- c("blue", "orange", "red")
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
    n.reads.consistant <- n.reads.total-n.reads.inconsistent.plus-n.reads.inconsistent.minus
    reads.label = c(signif(n.reads.consistant*100/n.reads.total), signif(n.reads.inconsistent.plus*100/n.reads.total)
                    , signif(n.reads.inconsistent.minus*100/n.reads.total))
    
    png(filename=paste(output.name,"_nb_total_reads.png"), width = 1024, height = 1024)
    pie(c(n.reads.consistant, n.reads.inconsistent.plus, n.reads.inconsistent.minus), labels = paste("Percentage: ",reads.label, 
                                                                                                     "\nAmount: ", c(n.reads.consistant, n.reads.inconsistent.plus, n.reads.inconsistent.minus)), col = tot.colors
        , main = "Reads Taxonomy Status after classification")
    legend("topleft", tot.names, cex = 1.5, fill = tot.colors)
    dev.off()
    
    # Contigs
    n.contigs.inconsistent <- length(unique(inconsistency_info$ContigID))
    n.contigs.consistant <- n.contigs.total-n.contigs.inconsistent.plus-n.contigs.inconsistent.minus
    contigs.label <-  c(signif(n.contigs.consistant*100/n.contigs.total), signif(n.contigs.inconsistent.plus*100/n.contigs.total)
                        , signif(n.contigs.inconsistent.minus*100/n.reads.total))
    
    png(filename=paste(output.name,"_nb_total_contigs.png"), width = 1024, height = 1024)
    pie(c(n.contigs.consistant, n.contigs.inconsistent.plus, n.contigs.inconsistent.minus), labels = paste("Percentage: ", contigs.label,
                                                                                                           "\nAmount: ", c(n.contigs.consistant, n.contigs.inconsistent.plus, n.contigs.inconsistent.minus)), col = tot.colors, 
        main = "Contigs Taxonomy Status after classification")
    legend("topleft", tot.names, cex = 1.5, fill = tot.colors)
    dev.off()
  }
}