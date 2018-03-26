# RACC R script
# Made by: Pablo Rodriguez Brazzarola



args <- commandArgs(trailingOnly = TRUE)

if (length(args)!=3) {
  stop("***USAGE*** Rscript RACC_script.R <RACC_results_directory> <RACC_report_directory> <Absolute_RACKit_src_r_functions>\n", call.=FALSE)
}

relative.path <- getwd()
### --- Inputs ---
results.directory <- paste(relative.path, args[1], sep = "")
original.distro.file <-  paste(results.directory, "refdb_distro.csv", sep = "")
read.distro.file <- paste(results.directory, "species-read_count.csv", sep = "")
contig.distro.file <- paste(results.directory, "species-contig_count.csv", sep = "")
inconsistency.solver.output <- paste(results.directory, "inc_solver.out", sep = "")
coverage.info.path <- paste(results.directory, "coverage.info", sep = "")
inconsistencies.found.output <- paste(results.directory, "inc-finder.out", sep = "")
statistical.measurements.output <- paste(results.directory, "rac-statistical_measurements.csv", sep = "")

output.directory <- paste(relative.path, args[2], sep = "")

source(args[3])

### --- Execution ---
processed_distro_data <- preprocess_distribution_data(original.distro.file, read.distro.file, contig.distro.file)

### Original Distro
# I: original.distro.file
# O: Plot
original.distro.output <- paste(output.directory, "IMG-original_ditribution_plot.png", sep = "")

original_distribution_plot(processed_distro_data$orig.descendant.distribution,
                           original.distro.output)

### Distro Comparison
# I: original.distro.file, read.distro.file, contig.distro.file
# O: Plot
distro.comparison.output <- paste(output.directory, "IMG-compared_ditribution_plot.png", sep = "")
distribution_comparison_plot(processed_distro_data$orig.descendant.distribution,
                             processed_distro_data$reads.descendant.distribution,
                             processed_distro_data$contigs.descendant.distribution,
                             distro.comparison.output)
### RMSE
# I: original.distro.file, read.distro.file, contig.distro.file
# O: Table CSV
rmse.results <- rmse_comparison(processed_distro_data$orig.descendant.distribution,
                                processed_distro_data$reads.descendant.distribution,
                                processed_distro_data$contigs.descendant.distribution)
rmse.results
### Inconsistency Resolution
# I: inc-solver.out
# O: Plot
inconsistency.resolution.output <- paste(output.directory, "IMG-inconsistency_solver.png", sep = "")
solved.inconsistencies.results <- inconsistency_resolution(inconsistency.solver.output,
                                                           inconsistency.resolution.output)
### DB Coverage Comparison
# I: coverage.info
# O: Table CSV
coverage.comparison.results <- coverage_comparison(coverage.info.path)
coverage.comparison.results