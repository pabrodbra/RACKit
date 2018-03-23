# Made by Pablo Rodriguez Brazzarola
# RACC Results Processing

setwd("~/Code/work/software/RACKit/src/r")
source("RACC_functions.R")
### Parameter introduction
# TODO

setwd("C:/Users/Blinsky.Blinsk/Documents/- UMA/TFG/RACKIT/results/rac_test/results")

original.distro.file = "" #<- read.csv2("grinder-ranks.txt", sep = '\t', header = TRUE)
read.distro.file = "" # TODO
contig.distro.file = "" # TODO
inconsistency.solver.output <- "inc_solver.out"
coverage.info.path <- "coverage.info"

processed_distro_data <- preprocess_distribution_data(original.distro.file, read.distro.file, contigs.distro.file)

### Original Distro
# I: original.distro.file
# O: Plot
original_distribution_plot(processed_distro_data$orig.descendant.distribution,
                           "original_ditribution_plot.png")

### Distro Comparison
# I: original.distro.file, read.distro.file, contig.distro.file
# O: Plot
distribution_comparison_plot(processed_distro_data$orig.descendant.distribution,
                           processed_distro_data$reads.descendant.distribution,
                           processed_distro_data$contigs.descendant.distribution,
                           "compared_ditribution_plot.png")
### RMSE
# I: original.distro.file, read.distro.file, contig.distro.file
# O: Table CSV
rmse.results <- rmse_comparison(processed_distro_data$orig.descendant.distribution,
                                processed_distro_data$reads.descendant.distribution,
                                processed_distro_data$contigs.descendant.distribution)

### Inconsistency Resolution
# I: inc-solver.out
# O: Plot
solved.inconsistencies.results <- inconsistency_resolution(inconsistency.solver.output,
                                                           "IMG-inconsistency_solver.png")
### DB Coverage Comparison
# I: coverage.info
# O: Table CSV
coverage.comparison.results <- coverage_comparison(coverage.info.path)

### Inconsistencies Found
# I: inc-finder.out
# O: Table CSV

### Statistical Measurements (TP/FP/TN/FN)
# I: rac-statistical_measurements.csv
# O: Table CSV