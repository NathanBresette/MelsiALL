#!/usr/bin/env Rscript
# Script to remove tables, figures, and non-text sections for word count

library(stringr)

# Read the manuscript
input_file <- "manuscript/MeLSI_Research_Paper_mSystems.md"
output_file <- "manuscript/MeLSI_Research_Paper_mSystems_wordcount.md"

text <- readLines(input_file)

# Track what we're removing
removed_sections <- c()

# Remove tables (Table 1-6)
# Pattern: **Table X.** followed by markdown table until next section or paragraph
tables_removed <- 0
in_table <- FALSE
table_start <- NULL
new_text <- c()
i <- 1

while (i <= length(text)) {
  line <- text[i]
  
  # Check if this is a table header
  if (str_detect(line, "^\\*\\*Table [0-9]+\\.")) {
    in_table <- TRUE
    table_start <- i
    # Skip the header line
    i <- i + 1
    next
  }
  
  # If we're in a table, check for end conditions
  if (in_table) {
    # End of table: empty line followed by non-table content, or section header, or figure
    if (str_detect(line, "^## |^### |^\\*\\*Figure|^!\\[|^\\\\begin\\{figure\\}")) {
      in_table <- FALSE
      # Don't include this line yet, process it in next iteration
      i <- i - 1
      next
    }
    
    # Skip markdown table lines (starting with |)
    if (str_detect(line, "^\\|")) {
      i <- i + 1
      next
    }
    
    # Skip table footnotes (lines starting with \noindent Abbreviations or \noindent that follow tables)
    if (str_detect(line, "^\\\\noindent Abbreviations|^\\\\noindent.*Abbreviations")) {
      # Check if next line is also part of table context
      i <- i + 1
      next
    }
    
    # If we hit a paragraph that's not part of table, we're done with table
    if (nchar(trimws(line)) > 0 && !str_detect(line, "^\\\\noindent|^\\s*$")) {
      in_table <- FALSE
      # Process this line normally
    } else {
      # Empty line or continuation, skip
      i <- i + 1
      next
    }
  }
  
  # Add line if not in table
  if (!in_table) {
    new_text <- c(new_text, line)
  }
  
  i <- i + 1
}

text <- new_text

# Remove figure captions and image references
# Pattern: ![](figures/...) or \begin{figure}...\end{figure} followed by Figure X caption
new_text <- c()
i <- 1
skip_next_n <- 0

while (i <= length(text)) {
  line <- text[i]
  
  # Skip image references
  if (str_detect(line, "^!\\[|^\\\\begin\\{figure\\}")) {
    # Skip until \end{figure} if LaTeX figure
    if (str_detect(line, "^\\\\begin\\{figure\\}")) {
      while (i <= length(text) && !str_detect(text[i], "^\\\\end\\{figure\\}")) {
        i <- i + 1
      }
      # Skip the \end{figure} line too
      i <- i + 1
      next
    } else {
      # Regular markdown image, skip it
      i <- i + 1
      next
    }
  }
  
  # Skip figure captions (lines starting with **Figure X.**)
  if (str_detect(line, "^\\\\noindent \\\\footnotesize|^\\\\noindent.*\\\\footnotesize")) {
    # This might be start of figure caption, skip until \normalsize
    while (i <= length(text) && !str_detect(text[i], "^\\\\normalsize")) {
      i <- i + 1
    }
    # Skip the \normalsize line too
    i <- i + 1
    next
  }
  
  new_text <- c(new_text, line)
  i <- i + 1
}

text <- new_text

# Remove specific sections
sections_to_remove <- c(
  "^## SUPPLEMENTARY MATERIAL",
  "^## FUNDING",
  "^## ACKNOWLEDGMENTS",
  "^## AUTHOR CONTRIBUTIONS",
  "^## COMPETING INTERESTS",
  "^## ORCID",
  "^## AUTHOR AFFILIATIONS",
  "^## REFERENCES",
  "^## DATA AVAILABILITY"
)

new_text <- c()
i <- 1
skip_section <- FALSE

while (i <= length(text)) {
  line <- text[i]
  
  # Check if this is a section to remove
  should_remove <- any(sapply(sections_to_remove, function(pattern) {
    str_detect(line, pattern)
  }))
  
  if (should_remove) {
    skip_section <- TRUE
    # Skip until next ## section
    while (i <= length(text)) {
      if (i < length(text) && str_detect(text[i + 1], "^## ")) {
        break
      }
      i <- i + 1
    }
    skip_section <- FALSE
    next
  }
  
  if (!skip_section) {
    new_text <- c(new_text, line)
  }
  
  i <- i + 1
}

# Write the cleaned text
writeLines(new_text, output_file)

cat("Removed tables, figures, and non-text sections.\n")
cat("Original lines:", length(readLines(input_file)), "\n")
cat("Cleaned lines:", length(new_text), "\n")
