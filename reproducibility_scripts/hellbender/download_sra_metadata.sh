#!/bin/bash
# ==============================================================================
# Download SRA Metadata for Group Labels
# ==============================================================================
# This script helps download SRA metadata for SRP214545
# ==============================================================================

echo "================================================================================"
echo "Downloading SRA Metadata for SRP214545"
echo "================================================================================"
echo ""

# Method 1: Try to download via SRA Run Selector API
echo "Method 1: Attempting direct download from SRA..."
curl -s "https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP214545&go=go" > sra_page.html

if [ -f sra_page.html ] && [ -s sra_page.html ]; then
    echo "✅ Downloaded SRA page"
    # Look for download links
    grep -i "download\|metadata" sra_page.html | head -5
else
    echo "❌ Could not download SRA page"
fi

echo ""
echo "================================================================================"
echo "MANUAL DOWNLOAD INSTRUCTIONS:"
echo "================================================================================"
echo ""
echo "1. Open this URL in your browser:"
echo "   https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP214545"
echo ""
echo "2. Click the 'Run Selector' button (or it may open automatically)"
echo ""
echo "3. On the Run Selector page:"
echo "   - Click 'Download' button (top right)"
echo "   - Select 'Metadata' from dropdown"
echo "   - Save the file as 'sra_metadata.csv' in this directory"
echo ""
echo "4. The metadata file should contain columns like:"
echo "   - Run (SRR numbers)"
echo "   - SampleName"
echo "   - BioSample"
echo "   - disease/condition information"
echo ""
echo "5. Once downloaded, run: Rscript parse_metadata.R"
echo ""
