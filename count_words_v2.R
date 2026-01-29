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

# Remove tables and figure legends more carefully
cleaned_lines <- c()
in_table <- FALSE
in_figure_legend <- FALSE

for (i in seq_along(main_text)) {
  line <- main_text[i]
  
  # Check for figure legend start
  if (grepl("\\\\footnotesize", line)) {
    in_figure_legend <- TRUE
    next
  }
  
  # Check for figure legend end
  if (grepl("\\\\normalsize", line)) {
    in_figure_legend <- FALSE
    next
  }
  
  # Skip if in figure legend
  if (in_figure_legend) {
    next
  }
  
  # Check if this line is part of a table
  if (grepl("^\\|", line) || grepl("^\\|.*\\|", line)) {
    in_table <- TRUE
    next
  }
  
  # If we were in a table and hit a non-table line (not empty), we're out
  if (in_table && !grepl("^\\s*$", line) && !grepl("^\\|", line)) {
    in_table <- FALSE
  }
  
  # Skip if still in table
  if (in_table) {
    next
  }
  
  # Skip image references
  if (grepl("^!\\[", line)) {
    next
  }
  
  # Skip LaTeX figure environments
  if (grepl("\\\\begin\\{figure\\}", line) || grepl("\\\\end\\{figure\\}", line) || 
      grepl("\\\\centering", line) || grepl("\\\\includegraphics", line)) {
    next
  }
  
  cleaned_lines <- c(cleaned_lines, line)
}

# Combine and clean
full_text <- paste(cleaned_lines, collapse = " ")

# Remove LaTeX commands more carefully (keep content inside braces)
full_text <- gsub("\\\\[a-zA-Z]+\\{([^}]*)\\}", "\\1", full_text)  # \command{content} -> content
full_text <- gsub("\\\\[a-zA-Z]+", "", full_text)  # \command -> remove

# Remove markdown headers
full_text <- gsub("^#+\\s+", "", full_text)

# Remove markdown formatting but keep words
full_text <- gsub("\\*\\*([^*]+)\\*\\*", "\\1", full_text)  # Bold
full_text <- gsub("\\*([^*]+)\\*", "\\1", full_text)  # Italic
full_text <- gsub("`([^`]+)`", "\\1", full_text)  # Code
full_text <- gsub("\\[([^\\]]+)\\]\\([^\\)]+\\)", "\\1", full_text)  # Links

# Remove math expressions (but this might remove too much - let's be careful)
# Actually, let's keep simple math and just remove complex LaTeX math
full_text <- gsub("\\$\\$[^$]+\\$\\$", " ", full_text)  # Display math
full_text <- gsub("\\$[^$]+\\$", " ", full_text)  # Inline math

# Remove special formatting characters but keep alphanumeric and spaces
full_text <- gsub("[^a-zA-Z0-9\\s]", " ", full_text)

# Split and count words
words <- unlist(strsplit(full_text, "\\s+"))
words <- words[nchar(words) > 0]

word_count <- length(words)

cat("Word count (excluding References, tables, and figure legends):", word_count, "\n")
cat("Target: 5,000 words\n")
cat("Difference:", word_count - 5000, "words\n")
if (word_count > 5000) {
  cat("Status: OVER by", word_count - 5000, "words\n")
} else {
  cat("Status: UNDER by", 5000 - word_count, "words\n")
}

# Also show a sample of what we're counting
cat("\nFirst 200 characters of cleaned text:\n")
cat(substr(full_text, 1, 200), "\n")
