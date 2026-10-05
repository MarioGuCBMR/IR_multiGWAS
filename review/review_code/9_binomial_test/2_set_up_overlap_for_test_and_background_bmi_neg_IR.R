##############
#INTRODUCTION#
##############

#This code runs an R script that sets up the data for to replicate the  
#Lotta et al. 2017 Fig 3d Replication: Matched-Locus Binomial Enrichment
#Most 


# Method (Lotta et al. Online Methods):
#   - Background: GWAS Catalog SNPs (p<5e-8, European), pruned at r2<0.1
#   - Test: 53 IR-associated lead SNPs + r2>0.8 proxies
#   - Matching: proxy count + locus span + TSS distance
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

TEST_SNPS_FILE <- paste(input_path, "/neg_63_ir_input.txt", sep = "")

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
OUTPUT_DIR <- "/maps/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/review/review_output/5_enrichment_analyses/2_binomial_lotta/output/neg_63_ir"

dir.create(OUTPUT_DIR)

######################## 
#STEP 1: Load Test Loci#
########################

test_snps <- fread(TEST_SNPS_FILE)
setnames(test_snps, tolower(names(test_snps)))
# Expected columns: chr, pos, rsid (optional)
if (!"chr" %in% names(test_snps) || !"pos" %in% names(test_snps)) {
  stop("TEST_SNPS_FILE must have 'chr' and 'pos' columns (hg19 coordinates)")
}
test_snps[, chr := as.integer(chr)]
test_snps[, pos := as.integer(pos)]
test_snps <- test_snps[chr %in% 1:22 & !is.na(pos)]
cat(sprintf("  %d test SNPs loaded (hg19)\n", nrow(test_snps)))

# Look up UK10K r2>0.8 proxies from GARFIELD tags/r08
cat("  Looking up UK10K r2>0.8 proxies from GARFIELD...\n")
test_loci <- list()
n_found <- 0
for (i in seq_len(nrow(test_snps))) {
  ch <- test_snps$chr[i]
  p <- test_snps$pos[i]
  tags_file <- file.path(GARFIELD_DIR, "tags", "r08", paste0("chr", ch))
  if (!file.exists(tags_file)) next
  tags <- fread(tags_file, header = FALSE, showProgress = FALSE)
  setnames(tags, c("position", "proxies"))
  row <- tags[position == p]
  if (nrow(row) == 0) {
    # SNP not in UK10K — use just the lead SNP as its own locus
    proxy_pos <- p
  } else {
    proxy_pos <- as.integer(unlist(strsplit(row$proxies, ",")))
    proxy_pos <- c(p, proxy_pos)  # include lead SNP
  }
  locus_span <- max(proxy_pos) - min(proxy_pos)
  # TSS distance from maftssd
  maftssd_file <- file.path(GARFIELD_DIR, "maftssd", paste0("chr", ch))
  tss_dist <- NA
  if (file.exists(maftssd_file)) {
    maftssd <- fread(maftssd_file, header = FALSE, showProgress = FALSE)
    setnames(maftssd, c("position", "maf", "tss_dist"))
    mrow <- maftssd[position == p]
    if (nrow(mrow) > 0) tss_dist <- mrow$tss_dist
  }
  test_loci[[length(test_loci) + 1]] <- data.table(
    locus_id = i,
    chr = ch,
    lead_pos = p,
    proxy_count = length(proxy_pos),
    locus_span = locus_span,
    tss_dist = tss_dist,
    proxy_pos = list(proxy_pos)
  )
  n_found <- n_found + 1
}
test_loci_dt <- rbindlist(test_loci)
cat(sprintf("  %d/%d test loci characterized\n", n_found, nrow(test_snps)))
fwrite(test_loci_dt[, .(locus_id, chr, lead_pos, proxy_count, locus_span, tss_dist)],
       file.path(OUTPUT_DIR, "test_loci_features.tsv"))

################################################
#SECTION 2: Build Background from GWAS Catalog #
################################################

gwas <- fread(GWAS_CATALOG, na.strings = c("", "NA", "-"), quote = "")
cat(sprintf("  %d total associations\n", nrow(gwas)))

# Join with ancestry to filter European studies
ancestry <- fread(GWAS_ANCESTRY)

# Find European studies (column name may vary)
european_studies <- unique(ancestry$"STUDY ACCESSION"[which(ancestry$"BROAD ANCESTRAL CATEGORY" == "European")])

# --- Filter GWAS Catalog ---
# European studies + p < 5e-8 + autosomes with valid positions
bg <- gwas[`STUDY ACCESSION` %in% european_studies]
bg <- bg[PVALUE_MLOG > -log10(P_THRESHOLD)]   # p < 5e-8
bg <- bg[`CHR_ID` %in% as.character(1:22)]    # autosomes only, skip X/Y/MT/multi
bg[, chr := as.integer(`CHR_ID`)]
bg[, pos := as.integer(`CHR_POS`)]
bg <- bg[!is.na(pos)]

# Keep most significant per SNP position (deduplicate)
bg <- bg[order(-PVALUE_MLOG)]
bg <- bg[!duplicated(paste(chr, pos))]
bg <- bg[, .(chr, pos, pval_mlog = PVALUE_MLOG, snp = SNPS)]
cat(sprintf("  %d unique European GWAS-significant SNPs (hg38)
", nrow(bg)))

## liftOver hg38 -> hg19
cat("  Performing liftOver hg38 -> hg19...\n")
chain <- import.chain(LIFTOVER_CHAIN)
gr_hg38 <- GRanges(seqnames = paste0("chr", bg$chr),
                   ranges = IRanges(start = bg$pos, end = bg$pos),
                   snp = bg$snp, pval_mlog = bg$pval_mlog)

gr_hg19 <- unlist(liftOver(gr_hg38, chain))
bg_hg19 <- as.data.table(gr_hg19)
bg_hg19[, chr := as.integer(sub("chr", "", as.character(seqnames)))]
bg_hg19[, pos := start]
bg_hg19 <- bg_hg19[chr %in% 1:22]
bg_hg19 <- bg_hg19[, .(chr, pos, pval_mlog, snp)]
cat(sprintf("  %d SNPs after liftOver to hg19\n", nrow(bg_hg19))) #conversion is correct
fwrite(bg_hg19, file.path(OUTPUT_DIR, "background_snps_european_p5e8_hg19.tsv"))

#######################################################
# ==== SECTION 3: Prune Background at r2 < 0.1 =======#
#######################################################

#To run this we actually need to install the binaries first.
#Let's do this:

#install.packages("genetics.binaRies", repos = c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org"))
plink_bin <- genetics.binaRies::get_plink_binary()

#Good, now we can actually run this

clump_input <- dplyr::tibble(
    rsid = bg_hg19$snp,
    pval = 10^(-bg_hg19$pval_mlog),
    id   = "background"
)

clump_input = clump_input[which(duplicated(clump_input$rsid) == FALSE),] #as expected, the deduplication worked above

clumped <- ieugwasr::ld_clump_local(
    dat       = clump_input,
    clump_kb  = 1000,
    clump_r2  = 0.1,
    clump_p   = 1,           # keep all as index SNPs (already filtered to p<5e-8)
    bfile     = "/maps/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/review/review_raw_data/binomial_data/EUR",
    plink_bin = genetics.binaRies::get_plink_binary()
)

#Now let's recover the data

bg_pruned <- bg_hg19[which(bg_hg19$snp%in%clumped$rsid),]

#Let's remove Lotta, we want another background!
#Maybe we will remove this section depending on whether we 
#can replicate the results
#Lotta does not do this, but it makes sense to remove redundancy

test_proxy_chrpos <- unlist(mapply(function(pp, ch) {
  paste0(ch, ":", pp)
}, test_loci_dt$proxy_pos, test_loci_dt$chr))

bg_pruned$chr_pos = paste0(bg_pruned$chr, ":", bg_pruned$pos)

bg_pruned = bg_pruned[which(!(bg_pruned$chr_pos%in%test_proxy_chrpos)),]

cat(sprintf("  %d background loci after pruning\n", nrow(bg_pruned)))
fwrite(bg_pruned, file.path(OUTPUT_DIR, "background_pruned.tsv"))

######################################################
# ==== SECTION 4: Characterize Background Loci ======#
######################################################

#We are gonna characterize the background loci according to proxy count, span, TSS
#To properly compare the overlap with annotations!

bg_loci <- list()
for (ch in 1:22) {
  ch_snps <- bg_pruned[chr == ch]
  if (nrow(ch_snps) == 0) next
  # Load tags/r08 and maftssd for this chromosome
  tags_file <- file.path(GARFIELD_DIR, "tags", "r08", paste0("chr", ch))
  maftssd_file <- file.path(GARFIELD_DIR, "maftssd", paste0("chr", ch))
  tags <- if (file.exists(tags_file)) fread(tags_file, header = FALSE, showProgress = FALSE) else NULL
  maftssd <- if (file.exists(maftssd_file)) fread(maftssd_file, header = FALSE, showProgress = FALSE) else NULL
  if (!is.null(tags)) setnames(tags, c("position", "proxies"))
  if (!is.null(maftssd)) setnames(maftssd, c("position", "maf", "tss_dist"))

  for (j in seq_len(nrow(ch_snps))) {
    p <- ch_snps$pos[j]
    # Get r2>0.8 proxies
    proxy_pos <- p
    if (!is.null(tags)) {
      row <- tags[position == p]
      if (nrow(row) > 0) {
        proxy_pos <- c(p, as.integer(unlist(strsplit(row$proxies, ","))))
      }
    }
    # TSS distance
    tss_dist <- NA
    if (!is.null(maftssd)) {
      mrow <- maftssd[position == p]
      if (nrow(mrow) > 0) tss_dist <- mrow$tss_dist
    }
    bg_loci[[length(bg_loci) + 1]] <- data.table(
      locus_id = ch_snps$snp[j],
      chr = ch,
      lead_pos = p,
      proxy_count = length(proxy_pos),
      locus_span = max(proxy_pos) - min(proxy_pos),
      tss_dist = tss_dist,
      proxy_pos = list(proxy_pos)
    )
  }
}
bg_loci_dt <- rbindlist(bg_loci)
cat(sprintf("  %d background loci characterized\n", nrow(bg_loci_dt)))
fwrite(bg_loci_dt[, .(locus_id, chr, lead_pos, proxy_count, locus_span, tss_dist)],
       file.path(OUTPUT_DIR, "background_loci_features.tsv"))

##################################################
# ==== SECTION 5: Match Test to Background ======#
##################################################

# Standardize features across background

features <- c("proxy_count", "locus_span", "tss_dist")
bg_feat <- bg_loci_dt[, ..features]
bg_scaled <- scale(bg_feat)
# Scale test loci using background mean/sd
test_feat <- test_loci_dt[, ..features]
test_scaled <- scale(test_feat, center = attr(bg_scaled, "scaled:center"),
                     scale = attr(bg_scaled, "scaled:scale"))

# For each test locus, find K nearest background loci (Euclidean distance)
matched_pools <- vector("list", nrow(test_loci_dt))
for (i in seq_len(nrow(test_loci_dt))) {
  dists <- sqrt(rowSums((bg_scaled - matrix(test_scaled[i, ], nrow = nrow(bg_scaled),
             ncol = length(features), byrow = TRUE))^2))
  k <- min(K_MATCH, length(dists))
  matched_pools[[i]] <- order(dists)[1:k]
}
cat(sprintf("  Matched %d test loci to pools of %d background loci each\n",
            length(matched_pools), K_MATCH))

saveRDS(matched_pools, file.path(OUTPUT_DIR, "matched_nearest_bg_loci_2_test.RDS"))

##################################################
# ==== SECTION 6: Load Roadmap Annotations ======#
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

################################################
# ==== SECTION 7: Precompute Overlaps =========#
################################################

cat("\n[7/10] Precomputing overlaps (background loci x epigenomes)...\n")
# For each background locus, check overlap with each epigenome's enhancers
n_bg <- nrow(bg_loci_dt)
n_epi <- length(enhancer_grs)
epi_names <- names(enhancer_grs)

# Build GRanges for all background loci (all proxy positions per locus)
bg_gr_list <- lapply(seq_len(n_bg), function(i) {
  pp <- bg_loci_dt$proxy_pos[[i]]
  ch <- bg_loci_dt$chr[i]
  GRanges(seqnames = paste0("chr", ch),
          ranges = IRanges(start = pp, end = pp))
}) #this works

# Overlap matrix: [n_bg x n_epi] binary
overlap_bg <- matrix(0L, nrow = n_bg, ncol = n_epi)
colnames(overlap_bg) <- epi_names
for (e in seq_len(n_epi)) {
  print(e)
  gr_enh <- enhancer_grs[[e]]
  if (length(gr_enh) == 0) next
  for (i in seq_len(n_bg)) {
    #print(i)
    hits <- findOverlaps(bg_gr_list[[i]], gr_enh)
    if (length(hits) > 0) overlap_bg[i, e] <- 1L
  }
}
cat(sprintf("  Background overlap matrix: %d x %d\n", nrow(overlap_bg), ncol(overlap_bg)))

saveRDS(bg_gr_list, file.path(OUTPUT_DIR, "bg_4_overlap.RDS"))
saveRDS(overlap_bg,  file.path(OUTPUT_DIR, "observed_bg_overlap.RDS"))

# Also compute observed overlap for test loci
cat("  Computing test locus overlaps...\n")
test_gr_list <- lapply(seq_len(nrow(test_loci_dt)), function(i) {
  pp <- test_loci_dt$proxy_pos[[i]]
  ch <- test_loci_dt$chr[i]
  GRanges(seqnames = paste0("chr", ch),
          ranges = IRanges(start = pp, end = pp))
})
observed_overlap <- integer(n_epi)
names(observed_overlap) <- epi_names
for (e in seq_len(n_epi)) {
  print(e)
  gr_enh <- enhancer_grs[[e]]
  if (length(gr_enh) == 0) next
  for (i in seq_len(nrow(test_loci_dt))) {
    hits <- findOverlaps(test_gr_list[[i]], gr_enh)
    if (length(hits) > 0) observed_overlap[e] <- observed_overlap[e] + 1L
  }
}

#Let's save the data just in case:

saveRDS(test_gr_list, file.path(OUTPUT_DIR, "test_4_overlap.RDS"))
saveRDS(observed_overlap,  file.path(OUTPUT_DIR, "observed_test_overlap.RDS"))
