# SKIOME Data Loading Status

## ✅ Data Successfully Loaded

**Dataset:** PRJNA554499 (SRP214545)  
**Source:** MGnify processed data (MGYS00005299)  
**File:** `extdata/SRP214545_taxonomy_abundances_SSU_v5.0.tsv`

### Data Summary
- **Samples:** 511 individual samples (SRR numbers)
- **Taxa:** 1,856 taxonomic assignments
- **Format:** Samples × Taxa abundance matrix (transposed from MGnify format)
- **Data Type:** 16S rRNA SSU taxonomic abundances

### Sample IDs
Sample IDs are SRR numbers (e.g., SRR9678960, SRR9678970, etc.)

## ⚠️ Pending: Metadata for Group Labels

**Status:** Need to identify which samples belong to which group:
- Atopic Dermatitis
- Psoriasis  
- Healthy controls (if present)

### Options to Get Metadata:
1. **SRA Run Selector:** Download metadata table from https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP214545
2. **Publication Supplementary:** Check Nature Communications paper (DOI: 10.1038/s41467-019-12253-y) for sample metadata
3. **BioSample Records:** Query individual SRR numbers to get BioSample attributes

### Next Steps:
1. Download SRA metadata table (Run Selector → Download → Metadata)
2. Match SRR numbers to group labels
3. Update `load_skiome_data.R` to include group labels in metadata
4. Run MeLSI validation

## Files Created:
- `load_skiome_data.R` - Data loading script
- `skiome_data_loaded.RData` - Loaded data (counts + metadata structure)

## Data Location:
- **Abundance table:** `/Users/nathanbresette/Desktop/MeLSI/extdata/SRP214545_taxonomy_abundances_SSU_v5.0.tsv`
- **Loaded data:** `/Users/nathanbresette/Desktop/MeLSI/reproducibility_scripts/hellbender/skiome_data_loaded.RData`
