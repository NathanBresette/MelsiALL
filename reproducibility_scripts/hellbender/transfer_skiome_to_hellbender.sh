#!/bin/bash
# ==============================================================================
# Transfer SKIOME Validation Files to Hellbender
# ==============================================================================

echo "================================================================================"
echo "Transferring SKIOME Validation Files to Hellbender"
echo "================================================================================"
echo ""

# Files to transfer
FILES=(
  "skiome_validation.R"
  "skiome_validation_job.sh"
  "SraRunTable.csv"
  "skiome_data_loaded.RData"
  "load_skiome_data.R"
)

REMOTE_DIR="~/melsi_simulations/hellbender/"
REMOTE_HOST="nbhtd@hellbender.rnet.missouri.edu"

echo "Files to transfer:"
for file in "${FILES[@]}"; do
  if [ -f "$file" ]; then
    echo "  ✓ $file"
  else
    echo "  ✗ $file (NOT FOUND)"
  fi
done
echo ""

echo "Transferring to: $REMOTE_HOST:$REMOTE_DIR"
echo ""

# Create remote directory if needed
ssh "$REMOTE_HOST" "mkdir -p $REMOTE_DIR" 2>/dev/null

# Transfer files
for file in "${FILES[@]}"; do
  if [ -f "$file" ]; then
    echo "Transferring $file..."
    scp "$file" "$REMOTE_HOST:$REMOTE_DIR"
    if [ $? -eq 0 ]; then
      echo "  ✓ Success"
    else
      echo "  ✗ Failed"
    fi
  fi
done

echo ""
echo "================================================================================"
echo "Transfer Complete!"
echo "================================================================================"
echo ""
echo "Next steps:"
echo "1. SSH into Hellbender: ssh $REMOTE_HOST"
echo "2. cd $REMOTE_DIR"
echo "3. sbatch skiome_validation_job.sh"
echo ""
