missing <- c()
for (i in 1:100) {
  if (!file.exists(paste0("table3_sim_", i, ".csv"))) {
    missing <- c(missing, i)
  }
}
cat("Missing files:", paste(missing, collapse=", "), "\n")
cat("Total missing:", length(missing), "\n")
