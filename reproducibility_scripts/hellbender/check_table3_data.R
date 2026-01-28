# Check Table 3 files for data quality
files <- list.files(pattern="table3_sim_.*\\.csv")
if (length(files) == 0) {
  cat("No Table 3 files found\n")
} else {
  cat("Found", length(files), "Table 3 files\n\n")
  
  # Check first file
  data <- read.csv(files[1], stringsAsFactors = FALSE)
  cat("=== COLUMN NAMES ===\n")
  cat(paste(colnames(data), collapse=", "), "\n\n")
  
  cat("=== SAMPLE ROW (first file) ===\n")
  print(data[1, ])
  cat("\n")
  
  # Check for NA values in traditional methods
  cat("=== CHECKING FOR NA VALUES IN TRADITIONAL METHODS ===\n")
  trad_cols <- c("Euclidean_F", "BrayCurtis_F", "Jaccard_F", "WeightedUniFrac_F", "UnweightedUniFrac_F")
  for (col in trad_cols) {
    if (col %in% colnames(data)) {
      na_count <- sum(is.na(data[[col]]))
      cat(col, ":", na_count, "NA values out of", nrow(data), "rows\n")
    } else {
      cat(col, ": COLUMN NOT FOUND\n")
    }
  }
  
  # Check multiple files
  cat("\n=== CHECKING MULTIPLE FILES ===\n")
  sample_files <- files[1:min(5, length(files))]
  for (f in sample_files) {
    d <- read.csv(f, stringsAsFactors = FALSE)
    has_all_trad <- all(c("Euclidean_F", "BrayCurtis_F", "Jaccard_F", "WeightedUniFrac_F", "UnweightedUniFrac_F") %in% colnames(d))
    cat(basename(f), ":", ifelse(has_all_trad, "Has all 5 traditional methods", "Missing some methods"), "\n")
  }
}
