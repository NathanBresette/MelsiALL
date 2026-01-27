#!/usr/bin/env Rscript
# Convert PostScript files to PNG

# Function to convert PS to PNG
convert_ps_to_png <- function(ps_file, png_file) {
  # Read PostScript file
  ps_content <- readLines(ps_file)
  
  # Create a temporary PDF first (R can handle PDF better)
  # Actually, let's use grDevices::postscript and then convert
  
  # Alternative: Use system command with ghostscript if available
  # Or use R's built-in capabilities
  
  # For now, let's try using R's postscript device capabilities
  # We'll need to render the PS file
  
  cat("Converting", ps_file, "to", png_file, "\n")
  
  # Try using system ghostscript if available
  gs_cmd <- Sys.which("gs")
  if (gs_cmd != "") {
    cmd <- paste0(
      gs_cmd, 
      " -dNOPAUSE -dBATCH -sDEVICE=png16m -r300 -dGraphicsAlphaBits=4 -dTextAlphaBits=4 ",
      "-sOutputFile=", png_file, " ", ps_file
    )
    system(cmd)
    if (file.exists(png_file)) {
      cat("Successfully converted using ghostscript\n")
      return(TRUE)
    }
  }
  
  # Fallback: Try using ImageMagick convert if available
  convert_cmd <- Sys.which("convert")
  if (convert_cmd != "") {
    cmd <- paste0(convert_cmd, " -density 300 ", ps_file, " ", png_file)
    system(cmd)
    if (file.exists(png_file)) {
      cat("Successfully converted using ImageMagick\n")
      return(TRUE)
    }
  }
  
  cat("ERROR: Could not convert. Need ghostscript or ImageMagick installed.\n")
  return(FALSE)
}

# Convert DietSwap figures
ps_dir <- "reproducibility_scripts/hellbender/figures"
png_dir <- "manuscript/figures"

# Convert PCoA
convert_ps_to_png(
  file.path(ps_dir, "dietswap_pcoa.ps"),
  file.path(png_dir, "dietswap_pcoa.png")
)

# Convert VIP combined
convert_ps_to_png(
  file.path(ps_dir, "dietswap_vip_combined.ps"),
  file.path(png_dir, "dietswap_vip_combined.png")
)

cat("\nConversion complete!\n")
