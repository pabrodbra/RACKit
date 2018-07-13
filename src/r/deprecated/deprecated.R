# Deprecated
join_coverage_comparisons <- function(vector.coverage.comparison){
  ret <- setNames(data.frame(matrix(ncol = 6, nrow = 0)), 
                  c("read.only.cov", "contig.only.cov", "read.total.cov", 
                    "contig.total.cov", "shared.cov", "not.cov"))
  for (cov.comp in vector.coverage.comparison){
    ret <- rbind(ret, cov.comp, stringsAsFactors = FALSE)
  }
  return(ret)
}
# Deprecated
join_rmse_comparisons <- function(vector.rmse.comparison){
  ret <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), 
                  c("reads.rmse", "contigs.rmse"))
  for (rmse.comp in vector.rmse.comparison){
    ret <- rbind(ret, rmse.comp, stringsAsFactors = FALSE)
  }
  return(ret)
}