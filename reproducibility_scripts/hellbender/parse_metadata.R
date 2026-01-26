
# Once you have the SRA metadata file, run this to extract group labels
# Expected format: CSV with columns including Run, SampleName, or disease/condition info

if (!file.exists("sra_metadata.csv")) {
  stop("Please download SRA metadata first from SRA Run Selector")
}

metadata <- read.csv("sra_metadata.csv", stringsAsFactors = FALSE)
print(colnames(metadata))

# Look for disease/condition columns
disease_cols <- grep("disease|condition|diagnosis|group|phenotype", 
                     colnames(metadata), ignore.case = TRUE, value = TRUE)
cat("Disease-related columns:", paste(disease_cols, collapse = ", "), "
")

# Match to our sample IDs
load("skiome_data_loaded.RData")
sample_ids <- rownames(counts)

# Create group labels
# This will need to be customized based on actual metadata structure

