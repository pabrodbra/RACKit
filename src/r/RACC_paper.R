# Made by Pablo Rodriguez Brazzarola
# RACC PAPER Results
rm(list = ls())
dev.off()

setwd("~/software/RACC/src/r")
source("RACC_functions.R")
setwd("~/results/rackit/")

###############################
### --- FSD ---
###############################

case="fsd/"

original.distro.file = "gastro_high.pa" #<- read.csv2("grinder-ranks.txt", sep = '\t', header = TRUE)
read.distro.file = paste0(case,"species-read_count.csv") # TODO
contig.distro.file = paste0(case,"species-contig_count.csv") # TODO
inconsistency.solver.output <- paste0(case,"inc_solver.out")
coverage.info.path <- paste0(case,"coverage.info")
inconsistencies.found.output <- paste0(case,"inc-finder.out")
statistical.measurements.output <- paste0(case,"rac-statistical_measurements.csv")

### --- Execution ---
#debug(preprocess_distribution_data)
processed_distro_data <- preprocess_distribution_data(original.distro.file, read.distro.file, contig.distro.file)

# t <- names(processed_distro_data$orig.descendant.distribution)
# strsplit(t[1], " ")[[1]][1]
# t2 <- lapply(seq_len(length(t)), function(x) strsplit(t[[x]], " ")[[1]][1])
# unique(t2)

dev.off()
### Original Distro
# I: original.distro.file
# O: Plot
original_distribution_plot(processed_distro_data$orig.descendant.distribution*100,
                           "IMG-original_ditribution_plot.png", PNG = FALSE, 
                           title="Gastrointestinal Tract Synthetic Dataset",
                           max.species = 20)

### Distro Comparison
# I: original.distro.file, read.distro.file, contig.distro.file
# O: Plot
distribution_comparison_plot(processed_distro_data$orig.descendant.distribution,
                             processed_distro_data$reads.descendant.distribution,
                             processed_distro_data$contigs.descendant.distribution,
                             "IMG-compared_ditribution_plot.png", PNG = FALSE, 
                             title="Gastrointestinal Tract Synthetic Dataset",
                             max.species = 20)

rmse.results <- rmse_comparison(processed_distro_data$orig.descendant.distribution,
                                processed_distro_data$reads.descendant.distribution,
                                processed_distro_data$contigs.descendant.distribution)
rmse.results

### Inconsistency Resolution
# I: inc-solver.out
# O: Plot
inconsistency.resolution.output <- paste(output.directory, "IMG-inconsistency_solver.png", sep = "")
debug(inconsistency_resolution)
solved.inconsistencies.results <- inconsistency_resolution(inconsistency.solver.output,
                                                           inconsistency.resolution.output, PNG=FALSE)
### DB Coverage Comparison
# I: coverage.info
# O: Table CSV
coverage.comparison.results <- coverage_comparison(coverage.info.path)
coverage.comparison.results

###############################
### --- FSD2 ---
###############################

case="fsd_v2/"

original.distro.file = paste0(case,"full_genomes_gastro.pa")
read.distro.file = paste0(case,"species-read_count.csv") # TODO
contig.distro.file = paste0(case,"species-contig_count.csv") # TODO
inconsistency.solver.output <- paste0(case,"inc_solver.out")
coverage.info.path <- paste0(case,"coverage.info")
inconsistencies.found.output <- paste0(case,"inc-finder.out")
statistical.measurements.output <- paste0(case,"rac-statistical_measurements.csv")

### --- Execution ---
debug(preprocess_distribution_data)
processed_distro_data <- preprocess_distribution_data(original.distro.file, read.distro.file, contig.distro.file)

### Original Distro
# I: original.distro.file
# O: Plot
undebug(original_distribution_plot)
original_distribution_plot(processed_distro_data$orig.descendant.distribution*100,
                           "IMG-original_ditribution_plot.png", PNG = FALSE, 
                           title="Gastrointestinal Tract Synthetic Dataset",
                           max.species = 20)

### Distro Comparison
# I: original.distro.file, read.distro.file, contig.distro.file
# O: Plot
undebug(distribution_comparison_plot)
distribution_comparison_plot(processed_distro_data$orig.descendant.distribution,
                             processed_distro_data$reads.descendant.distribution,
                             processed_distro_data$contigs.descendant.distribution,
                             "IMG-compared_ditribution_plot.png", PNG = FALSE, 
                             title="Gastrointestinal Tract Synthetic Dataset",
                             max.species = 20)

rmse.results <- rmse_comparison(processed_distro_data$orig.descendant.distribution,
                                processed_distro_data$reads.descendant.distribution,
                                processed_distro_data$contigs.descendant.distribution)
rmse.results

### Inconsistency Resolution
# I: inc-solver.out
# O: Plot
inconsistency.resolution.output <- paste(output.directory, "IMG-inconsistency_solver.png", sep = "")
solved.inconsistencies.results <- inconsistency_resolution(inconsistency.solver.output,
                                                           inconsistency.resolution.output, PNG=FALSE)
### DB Coverage Comparison
# I: coverage.info
# O: Table CSV
coverage.comparison.results <- coverage_comparison(coverage.info.path)
coverage.comparison.results

###############################
### --- SSD ---
###############################

case="ssd/"

original.distro.file = paste0(case,"gastro_high.pa")
read.distro.file = paste0(case,"species-read_count.csv") # TODO
contig.distro.file = paste0(case,"species-contig_count.csv") # TODO
inconsistency.solver.output <- paste0(case,"inc_solver.out")
coverage.info.path <- paste0(case,"coverage.info")
inconsistencies.found.output <- paste0(case,"inc-finder.out")
statistical.measurements.output <- paste0(case,"rac-statistical_measurements.csv")

