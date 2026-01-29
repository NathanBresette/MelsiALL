missing <- c()
for (i in 1:100) {
  file_name <- paste0("table3_sim_", i, ".csv")
  if (!file.exists(file_name)) {
    missing <- c(missing, i)
  }
}
cat("Missing files:", paste(missing, collapse=", "), "\n")
cat("Total missing:", length(missing), "\n")
cat("Files present:", 100 - length(missing), "\n")
