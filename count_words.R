# Count words excluding References, tables, and figure legends
library(stringr)

# Read the manuscript
text <- readLines("manuscript/MeLSI_Research_Paper_mSystems.md", warn = FALSE)

# Find where References section starts
ref_start <- which(grepl("^## REFERENCES", text))
if (length(ref_start) == 0) {
  ref_start <- length(text) + 1
}

# Extract text up to References
main_text <- text[1:(ref_start - 1)]

# Remove tables (markdown tables start with |)
# Tables are multi-line blocks, so we need to identify table blocks
in_table <- FALSE
cleaned_lines <- c()

for (i in seq_along(main_text)) {
  line <- main_text[i]
  
  # Check if this line starts a table (starts with |)
  if (grepl("^\\|", line)) {
    in_table <- TRUE
    next  # Skip table lines
  }
  
  # Check if we're still in a table (empty line or another table line)
  if (in_table) {
    if (grepl("^\\|", line) || grepl("^\\s*$", line)) {
      next  # Still in table
    } else {
      in_table <- FALSE  # End of table
    }
  }
  
  # Skip figure legends (between \footnotesize and \normalsize)
  if (grepl("\\\\footnotesize", line)) {
    next
  }
  if (grepl("\\\\normalsize", line)) {
    next
  }
  
  # Skip lines that are part of figure legends (between footnotesize and normalsize)
  # We'll handle this by tracking state
  if (!in_table) {
    cleaned_lines <- c(cleaned_lines, line)
  }
}

# Now remove figure legends more carefully
# Find all \footnotesize and \normalsize pairs
final_lines <- c()
skip_until_normalsize <- FALSE

for (line in cleaned_lines) {
  if (grepl("\\\\footnotesize", line)) {
    skip_until_normalsize <- TRUE
    next
  }
  if (grepl("\\\\normalsize", line)) {
    skip_until_normalsize <- FALSE
    next
  }
  if (!skip_until_normalsize) {
    final_lines <- c(final_lines, line)
  }
}

# Remove markdown formatting, LaTeX commands, and count words
clean_text <- paste(final_lines, collapse = " ")

# Remove LaTeX commands (\command{...} or \command)
clean_text <- gsub("\\\\[a-zA-Z]+(\\{[^}]*\\})?", "", clean_text)

# Remove markdown formatting
clean_text <- gsub("\\*\\*([^*]+)\\*\\*", "\\1", clean_text)  # Bold
clean_text <- gsub("\\*([^*]+)\\*", "\\1", clean_text)  # Italic
clean_text <- gsub("`([^`]+)`", "\\1", clean_text)  # Code
clean_text <- gsub("\\[([^\\]]+)\\]\\([^\\)]+\\)", "\\1", clean_text)  # Links
clean_text <- gsub("!\\[([^\\]]*)\\]\\([^\\)]+\\)", "", clean_text)  # Images
clean_text <- gsub("\\$[^$]+\\$", "", clean_text)  # Math mode
clean_text <- gsub("\\$\\$[^$]+\\$\\$", "", clean_text)  # Display math

# Remove special characters but keep spaces
clean_text <- gsub("[^a-zA-Z0-9\\s]", " ", clean_text)

# Split into words and count
words <- strsplit(clean_text, "\\s+")[[1]]
words <- words[words != ""]

word_count <- length(words)

cat("Word count (excluding References, tables, and figure legends):", word_count, "\n")
cat("Target: 5,000 words\n")
cat("Difference:", word_count - 5000, "words\n")
if (word_count > 5000) {
  cat("Status: OVER by", word_count - 5000, "words\n")
} else {
  cat("Status: UNDER by", 5000 - word_count, "words\n")
}
