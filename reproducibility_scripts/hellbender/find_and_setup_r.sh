#!/bin/bash
# ==============================================================================
# Find and Setup R on Hellbender
# ==============================================================================

echo "================================================================================"
echo "Finding R on Hellbender"
echo "================================================================================"
echo ""

# Check if R is already available
if command -v R &> /dev/null; then
    echo "✅ R is already available:"
    R --version | head -1
    exit 0
fi

echo "R not found in PATH. Searching for R modules..."
echo ""

# Try to find R modules
echo "1. Checking module system for R..."
R_MODULES=$(module avail 2>&1 | grep -iE "^[[:space:]]*[rR]/" 2>/dev/null || echo "")

if [ -n "$R_MODULES" ]; then
    echo "   Found R modules:"
    echo "$R_MODULES" | sed 's/^/     /'
    echo ""
    echo "   Try loading one:"
    FIRST_R=$(echo "$R_MODULES" | head -1 | awk '{print $1}')
    echo "     module load $FIRST_R"
    module load $FIRST_R 2>/dev/null
    if command -v R &> /dev/null; then
        echo "   ✅ Successfully loaded $FIRST_R"
        R --version | head -1
        exit 0
    fi
else
    echo "   No R modules found in module system"
fi

echo ""
echo "2. Checking for R in common locations..."
for path in /usr/bin/R /usr/local/bin/R /opt/R /cluster/software/*/R; do
    if [ -f "$path" ]; then
        echo "   ✅ Found R at: $path"
        echo "   Add to PATH or create symlink"
        exit 0
    fi
done

echo ""
echo "================================================================================"
echo "⚠️  R Not Found"
echo "================================================================================"
echo ""
echo "R does not appear to be available on Hellbender."
echo ""
echo "Options:"
echo "  1. Contact HPC support to request R installation"
echo "  2. Install R in your home directory (may require compilation)"
echo "  3. Use conda/miniconda to install R (if available)"
echo ""
echo "To check for conda:"
echo "  which conda"
echo "  module avail conda"
echo ""
echo "For help, contact: HPC Support or check documentation at:"
echo "  https://docs.itrss.umsystem.edu/pub/hpc/hellbender"
echo ""
