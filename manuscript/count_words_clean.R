#!/usr/bin/env Rscript
# Word count script for cleaned manuscript

library(stringr)

# Read the cleaned manuscript
text <- readLines("MeLSI_Research_Paper_mSystems_wordcount.md")

# Remove LaTeX commands and formatting
# Remove \noindent, \footnotesize, \normalsize, etc.
text <- str_remove_all(text, "\\\\noindent|\\\\footnotesize|\\\\normalsize|\\\\textcolor\\{[^}]+\\}|\\\\underline\\{[^}]+\\}")

# Remove markdown formatting but keep words
text <- str_remove_all(text, "\\*\\*|\\*|#+\\s*|\\$\\^|\\$|\\{|\\}")

# Remove URLs but keep text
text <- str_remove_all(text, "https?://[^\\s]+")

# Remove email addresses
text <- str_remove_all(text, "[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\\.[a-zA-Z]{2,}")

# Remove DOIs
text <- str_remove_all(text, "DOI:?\\s*[0-9.]+/[^\\s]+")

# Remove figure/image references
text <- str_remove_all(text, "!\\[.*?\\]\\(.*?\\)")

# Remove LaTeX math delimiters but keep content
text <- str_replace_all(text, "\\$\\$([^$]+)\\$\\$", "\\1")
text <- str_replace_all(text, "\\$([^$]+)\\$", "\\1")

# Remove special characters that aren't part of words
text <- str_replace_all(text, "[^a-zA-Z0-9\\s]", " ")

# Split into words and count
all_words <- unlist(str_split(text, "\\s+"))
all_words <- all_words[nchar(all_words) > 0]  # Remove empty strings

# Count words
word_count <- length(all_words)

cat("Total word count:", word_count, "\n")
cat("Total characters:", sum(nchar(all_words)), "\n")

# Write word count to file
writeLines(c(
  paste("Word count:", word_count),
  paste("Character count:", sum(nchar(all_words))),
  "",
  "Note: This count excludes:",
  "- Tables",
  "- Figure captions",
  "- References",
  "- Funding, Acknowledgments, Author sections",
  "- URLs, DOIs, email addresses",
  "- LaTeX formatting commands"
), "word_count.txt")

cat("\nWord count saved to word_count.txt\n")
