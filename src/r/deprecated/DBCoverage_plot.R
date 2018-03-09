# UniseqDB Coverage

rm(list=ls())
setwd("C:/Users/Blinsky.Blinsk/Documents/- UMA/BitLab/Metagenome/Synthetic Dataset/Results")

coverage.info.path <- "Coverage.info"

coverage.info <- read.csv(coverage.info.path, header=TRUE, sep = '\t')
coverage.info <- coverage.info[c("ReadOnly", "BothYes", "ContigOnly", "BothNo")]
coverage.total <- sum(as.numeric(coverage.info))

coverage.label <- c(signif(coverage.info[1]*100/coverage.total), signif(coverage.info[2]*100/coverage.total)
                    , signif(coverage.info[3]*100/coverage.total), signif(coverage.info[4]*100/coverage.total))


tot.colors <- c("orange", "dark blue","dark green", "red")
tot.names <- c("Read Only", "Both", "Contig Only", "None")

png(filename="DB_coverage.png", width = 1024, height = 1024)
pie(as.numeric(coverage.info), labels = paste("Percentage: ", coverage.label, "\nNucleotides: ", coverage.info), col = tot.colors
    , main = "Percentage of Database Covered")
legend("bottomleft", tot.names, cex = 1.5, fill = tot.colors)
dev.off()
