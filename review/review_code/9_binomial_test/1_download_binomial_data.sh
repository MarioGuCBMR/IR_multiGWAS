#!/bin/bash
##############
#INTRODUCTION#
##############

#This code downloads the data to replicate Lotta et al. 2017 Fig 3d binomial enrichment replication

# Downloads:
#   1. GARFIELD data (5.9 GB) — LD tags + MAF/TSS distance (UK10K-based)
#   2. GWAS Catalog associations (hg38, zipped)
#   3. GWAS Catalog ancestry file (for European filtering)
#   4. liftOver chain file (hg38 -> hg19)

###################################
#STEP 1: set up download directory#
###################################

DATA_DIR="/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/review/review_raw_data/binomial_data"
mkdir -p "$DATA_DIR"
cd "$DATA_DIR"

echo "============================================"
echo "Downloading binomial enrichment data"
echo "Output directory: $(pwd)"
echo "============================================"

# --- 1. GARFIELD data (5.9 GB compressed, ~83 GB full) ---
https://www.ebi.ac.uk/birney-srv/GARFIELD/package-v2/garfield-data.tar.gz

echo "  Extracting maftssd/ and tags/ only (skipping annotation/, pval/)..."
# Extract only needed subdirectories to save disk space
tar -xzf garfield-data.tar.gz garfield-data/maftssd/ 2>/dev/null || true
tar -xzf garfield-data.tar.gz garfield-data/tags/r01/ 2>/dev/null || true
tar -xzf garfield-data.tar.gz garfield-data/tags/r08/ 2>/dev/null || true
echo "  GARFIELD data extracted to garfield-data/"

# --- 2. GWAS Catalog associations (hg38) ---
echo ""
echo "[2/4] Downloading GWAS Catalog associations..."
wget -c "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated-full.zip" \
    -O gwas-catalog-associations.zip
unzip -o gwas-catalog-associations.zip
# Find the TSV file (name may vary by release)
GWAS_TSV=$(find . -maxdepth 1 -name "gwas-catalog-associations*.tsv" | head -1)
echo "  GWAS Catalog associations: $GWAS_TSV"

# --- 3. GWAS Catalog ancestry ---
echo ""
echo "[3/4] Downloading GWAS Catalog ancestry..."
wget -c "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-ancestry.tsv" \
    -O gwas-catalog-ancestry.tsv
echo "  GWAS Catalog ancestry: gwas-catalog-ancestry.tsv"

# --- 4. liftOver chain file (hg38 -> hg19) ---
echo ""
echo "[4/4] Downloading liftOver chain file (hg38 -> hg19)..."
wget -c "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"
gunzip -f hg38ToHg19.over.chain.gz
echo "  liftOver chain: hg38ToHg19.over.chain"

# --- 5. Download data to perform LD-clumping for the bg:

wget "https://github.com/MRCIEU/genetics.binaRies/raw/master/binaries/Linux/plink"
wget "http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz"
