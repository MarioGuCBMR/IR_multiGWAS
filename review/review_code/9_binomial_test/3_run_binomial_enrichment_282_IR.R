##############
#INTRODUCTION#
##############

#This code runs an R script that sets up the data for to replicate the  
#Lotta et al. 2017 Fig 3d Replication: Matched-Locus Binomial Enrichment
#Most 


# Method (Lotta et al. Online Methods):
#   - 1M permutations to estimate expected overlap proportion
#   - Binomial test per epigenome
#   - Annotations: ROADMAP 15-state model, 98 H3K27ac epigenomes
#   - Active enhancers: 6_EnhG + 7_Enh

###################
#Loading libraries#
###################

library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(svglite)
library(tidyverse)

#######################
#Setting up parameters#
#######################

input_path = "/maps/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/review/review_output/5_enrichment_analyses/2_binomial_lotta/input"
raw_data_path = "/maps/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/review/review_raw_data"
additional_path = "/maps/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/output/3_go_shifter/curated_atac_seq"

TEST_SNPS_FILE <- paste(input_path, "/all_282_ir_input.txt", sep = "")

# Data paths (from download_binomial_data.sh)
GARFIELD_DIR    <- paste(raw_data_path, "/binomial_data/garfield-data", sep = "")
GWAS_CATALOG    <- paste(raw_data_path, "/binomial_data/gwas-catalog-download-associations-alt-full.tsv", sep = "")
GWAS_ANCESTRY   <- paste(raw_data_path, "/binomial_data/gwas-catalog-ancestry.tsv", sep = "")
LIFTOVER_CHAIN  <- paste(raw_data_path,"/binomial_data/hg38ToHg19.over.chain", sep = "")

# Roadmap 15-state annotation files (one per epigenome, hg19)
ROADMAP_DIR      <- paste(raw_data_path, "/binomial_data/roadmap_15", sep = "")
EIDS_FILE        <- paste(raw_data_path, "/binomial_data/eids_98_H3K27ac.txt", sep = "")
EID_METADATA     <- paste(raw_data_path, "/binomial_data/EID_metadata.tab", sep = "")

#Let's work only with 

# Analysis parameters
N_PERM          <- 1e6       # Number of permutations (Lotta used 1M)
K_MATCH         <- 20        # Matched background loci per test locus
ENHANCER_STATES <- c("6_EnhG", "7_Enh")  # 15-state active enhancers (excl. bivalent)
P_THRESHOLD     <- 5e-8      # GWAS significance threshold
LD_PRUNE_R2     <- 0.1       # Pruning threshold for background
LD_PROXY_R2     <- 0.8       # Proxy threshold for locus definition

# Output
OUTPUT_DIR <- "/maps/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/review/review_output/5_enrichment_analyses/2_binomial_lotta/output/all_282_ir"

dir.create(OUTPUT_DIR)

##################################################
# ==== SECTION 1: Load Roadmap Annotations ======#
##################################################

#Loading ROADMAP 15-state annotations...\n")
eids <- scan(EIDS_FILE, what = character(), quiet = TRUE)
cat(sprintf("  %d epigenomes to process\n", length(eids)))

# Load EID metadata for tissue names
eid_meta <- fread(EID_METADATA)
meta_cols <- names(eid_meta)
mnemonic_col <- grep("mnemonic", meta_cols, ignore.case = TRUE, value = TRUE)[1]
std_name_col <- grep("std_name", meta_cols, ignore.case = TRUE, value = TRUE)[1]
group_col <- grep("group", meta_cols, ignore.case = TRUE, value = TRUE)[1]

# Load enhancer annotations for each epigenome
enhancer_grs <- list()
for (eid in eids) {
  bed_file <- file.path(ROADMAP_DIR, paste0("/", eid, "_15_coreMarks_dense.bed.gz"))
  if (!file.exists(bed_file)) {
    # Try without .gz
    bed_file <- file.path(ROADMAP_DIR, paste0("/", eid, "_15_coreMarks_dense.bed"))
  }
  if (!file.exists(bed_file)) {
    cat(sprintf("    WARNING: %s not found, skipping\n", eid))
    next
  }
  bed <- fread(bed_file, skip = 1, header = FALSE)
  bed = bed %>%
      dplyr::select(V1, V2, V3, V4) #we only need the first 4
  
  setnames(bed, c("chr", "start", "end", "state"))
  enh <- bed[state %in% ENHANCER_STATES]
  if (nrow(enh) > 0) {
    enhancer_grs[[eid]] <- GRanges(
      seqnames = enh$chr,
      ranges = IRanges(start = enh$start, end = enh$end)
    )
  } else {
    enhancer_grs[[eid]] <- GRanges()
  }
}

#Let's add the data from Perrin et al!

adipose <- fread(file.path(additional_path, "peaks_adipose_tissue.bed.gz"), header = FALSE)
sgbs_day0 <- fread(file.path(additional_path, "peaks_sgbs_day0.bed.gz"), header = FALSE)
sgbs_day4 <- fread(file.path(additional_path, "peaks_sgbs_day4.bed.gz"), header = FALSE)
sgbs_day14 <- fread(file.path(additional_path, "peaks_sgbs_day14.bed.gz"), header = FALSE)

enhancer_grs[["Adipose (Perrin et al)"]] <- GRanges(
  seqnames = adipose$V1,
  ranges = IRanges(start = adipose$V2, end = adipose$V3)
)

enhancer_grs[["SGBS day 0 (Perrin et al)"]] <- GRanges(
  seqnames = sgbs_day0$V1,
  ranges = IRanges(start = sgbs_day0$V2, end = sgbs_day0$V3)
)

enhancer_grs[["SGBS day 4 (Perrin et al)"]] <- GRanges(
  seqnames = sgbs_day4$V1,
  ranges = IRanges(start = sgbs_day4$V2, end = sgbs_day4$V3)
)

enhancer_grs[["SGBS day 14 (Perrin et al)"]] <- GRanges(
  seqnames = sgbs_day14$V1,
  ranges = IRanges(start = sgbs_day14$V2, end = sgbs_day14$V3)
)

#######################################
#===== SECTION 2: Permutations =======#
#######################################

#Let's first load the data needed 

test_loci_dt = fread(file.path(OUTPUT_DIR, "test_loci_features.tsv"))
matched_pools = readRDS(file.path(OUTPUT_DIR, "matched_nearest_bg_loci_2_test.RDS"))
overlap_bg = readRDS(file.path(OUTPUT_DIR, "observed_bg_overlap.RDS"))
observed_overlap = readRDS(file.path(OUTPUT_DIR, "observed_test_overlap.RDS"))

#And also add some needed settings:

n_epi <- length(enhancer_grs)
epi_names <- names(enhancer_grs)
n_test <- nrow(test_loci_dt)

# Precompute: for each test locus i and epigenome e, the overlap values of its K matched bg loci
# pool_overlaps[[e]]: matrix [n_test x K_match] of 0/1

pool_overlaps <- lapply(seq_len(n_epi), function(e) {
  mat <- matrix(0L, nrow = n_test, ncol = K_MATCH)
  for (i in seq_len(n_test)) {
    pool <- matched_pools[[i]]
    mat[i, 1:length(pool)] <- overlap_bg[pool, e]
  }
  mat
})

# Sample 1M indices from each pool
sampled_idx <- matrix(sample(1:K_MATCH, n_test * N_PERM, replace = TRUE),
                      nrow = N_PERM, ncol = n_test)

# For each epigenome, compute permuted overlap counts
perm_counts <- matrix(0L, nrow = N_PERM, ncol = n_epi)
colnames(perm_counts) <- epi_names
for (e in seq_len(n_epi)) {
  po <- pool_overlaps[[e]]  # [n_test x K]
  counts <- integer(N_PERM)
  for (i in seq_len(n_test)) {
    counts <- counts + po[i, sampled_idx[, i]]
  }
  perm_counts[, e] <- counts
  if (e %% 10 == 0) cat(sprintf("    %d/%d epigenomes done\n", e, n_epi))
}

################################################
# ==== SECTION 9: Binomial Test ===============#
################################################

#Let's first get the results in a data table to compute the binomial test on it

results <- data.table(
  eid = epi_names,
  observed = as.integer(observed_overlap[epi_names]),
  n_test = n_test,
  expected = as.numeric(colMeans(perm_counts[, epi_names])),
  expected_prop = as.numeric(colMeans(perm_counts[, epi_names])) / n_test,
  observed_prop = as.numeric(observed_overlap[epi_names]) / n_test
)

# Binomial test: observed successes out of n_test, with expected_prob
n <- n_test  # scalar, no name collision

results[, binomial_p := mapply(function(obs, prob) {
  if (is.na(prob) || prob == 0 || prob == 1) return(NA)
  binom.test(as.integer(obs), n, prob, alternative = "greater")$p.value
}, observed, expected_prop)]

# Permutation p-value
results[, perm_p := sapply(epi_names, function(e) {
  (sum(perm_counts[, e] >= observed_overlap[e]) + 1) / (N_PERM + 1)
})]
results[, neg_log10_p := -log10(binomial_p)]
results[, enrichment := observed_prop / expected_prop]

# Add tissue metadata
results[, tissue := sapply(eid, function(e) {
  m <- eid_meta[get(names(eid_meta)[1]) == e]
  if (nrow(m) > 0 && !is.na(std_name_col)) return(m[[std_name_col]])
  return(e)
})]

results[, group := sapply(eid, function(e) {
  m <- eid_meta[get(names(eid_meta)[1]) == e]
  if (nrow(m) > 0 && !is.na(group_col)) return(m[[group_col]])
  return("Unknown")
})]

# Sort by significance
results <- results[order(binomial_p)]
fwrite(results, file.path(OUTPUT_DIR, "results_binomial.tsv"))
