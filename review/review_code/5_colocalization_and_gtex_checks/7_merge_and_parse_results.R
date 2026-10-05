##############
#INTRODUCTION#
##############

#This code leverages all the results from co-localization per lead, makes a supplementary table and, finally, compares the results!!

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

###################
#Loading functions#
###################

#' max() that returns NA when every value is NA (instead of -Inf with a warning)
.safe_max <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}

#' Convert empty strings ("") to NA_character_
.empty_to_na <- function(x) {
  x[x == ""] <- NA_character_
  x
}

compute_r2 <- function(lead, proxy, token, pop) {
  
  #STEP 0: return empty or 1 in special cases that do not require teting
  
  if (is.na(proxy)) return(NA_real_)
  if (lead == proxy) return(1.0)  # perfect LD with itself
  
  #STEP 1: run a trycatch the LDPair
  
  r2 <- tryCatch({
    ld <- LDlinkR::LDpair(lead, proxy, pop = pop, token = token)
    as.numeric(ld$r2)
  }, error = function(e) {
    warning("LDpair failed for ", lead, " vs ", proxy, ": ",
            conditionMessage(e), call. = FALSE)
    NA_real_
  })
  
  #STEP 3: return the info!
  
  return(r2)
}

add_coloc_to_irvariants_eqtls <- function(ir_df, coloc_df, tissue_name, ensdb = NULL, ir_snp_col = "Lead IR SNP", token = "04cad4ca4374", pop = "EUR", ld_delay = 0.5) {
  
  # STEP 0: prepare the column names according to the tissue
  col_pp       <- paste0(tissue_name, "_colo_pp")
  col_reg      <- paste0(tissue_name, "_regional_pp")
  col_gene     <- paste0(tissue_name, "_coloc_gene")
  col_sym      <- paste0(tissue_name, "_coloc_symbol")
  col_proxy    <- paste0(tissue_name, "_proxy")
  col_proxy_pp <- paste0(tissue_name, "_proxy_pp")
  col_proxy_r2 <- paste0(tissue_name, "_proxy_r2")
  col_high_conf <- paste0(tissue_name, "_high_conf_gene")
  
  # STEP 1: make sure that we have the right columns in the coloc res
  required_ir <- ir_snp_col
  required_colo <- c("lead_snp", "traits", "posterior_prob",
                     "regional_prob", "candidate_snp",
                     "posterior_explained_by_snp")
  
  if (!required_ir %in% names(ir_df))
    stop("ir_df must contain column '", required_ir, "'")
  missing_cols <- setdiff(required_colo, names(coloc_df))
  if (length(missing_cols) > 0)
    stop("coloc_df missing required columns: ",
         paste(missing_cols, collapse = ", "))
  
  # STEP 2: Deduplicate coloc results, though we have already done this before
  coloc_df <- unique(coloc_df)
  
  # STEP 3: Extract ENSG IDs from traits
  # traits looks like "iradjbmi, ENSG00000260063" or "None"
  coloc_df$ensg_list <- lapply(coloc_df$traits, function(x) {
    if (is.na(x)) return(character(0))
    unlist(regmatches(x, gregexpr("ENSG[0-9]{11}", x)))
  })
  coloc_df$has_gene <- lengths(coloc_df$ensg_list) > 0
  
  #################################################################
  # STEP 4: get co-localization results in a more structured shape#
  ##################################################################
  
  coloc_rows <- coloc_df[coloc_df$has_gene, ]
  
  if (nrow(coloc_rows) > 0) {
    
    #########################################################################
    # STEP 4.1: Compute LD r2 between lead_snp and candidate_snp via LDlinkR#
    #########################################################################
    
    if (!is.null(token)) {
      r2_cache <- new.env(hash = TRUE, parent = emptyenv())
      n_pairs <- nrow(coloc_rows)
      message("Computing LD r2 for ", n_pairs, " (lead, proxy) pair(s) via LDlinkR...")
      
      coloc_rows$proxy_r2 <- sapply(seq_len(n_pairs), function(i) {
        if (i %% 10 == 0) message("  ...pair ", i, " of ", n_pairs)
        compute_r2(
          lead   = coloc_rows$lead_snp[i],
          proxy  = coloc_rows$candidate_snp[i],
          token  = token,
          pop    = "EUR")
      })
      
    } else {
      warning("No LDlinkR token provided — proxy_r2 and high_conf_gene will be NA.",
              call. = FALSE)
      coloc_rows$proxy_r2 <- NA_real_
    }
    
    ########################################################################
    # STEP 4.2: Map ENSG -> gene symbol (batch, then per-row for high-conf)#
    ########################################################################
    
    all_ensgs <- unique(unlist(coloc_rows$ensg_list))
    
    #Let's transform them to symbols with EnsDb:
    
    if (is.null(ensdb)) {
      if (!requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE))
        stop("Install EnsDb.Hsapiens.v86: BiocManager::install('EnsDb.Hsapiens.v86')")
      ensdb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
    }
    
    #Let's make a map
    
    symbol_map <- AnnotationDbi::mapIds(
      ensdb,
      keys      = all_ensgs,
      column    = "SYMBOL",
      keytype   = "GENEID",
      multiVals = "first"
    )
    
    #Let's update coloc_rows so that we have the data fully mapped (not convering or pasting data per lead yet)
    coloc_rows$gene_symbol <- vapply(
      coloc_rows$ensg_list,
      function(ensgs) paste(unname(symbol_map[ensgs]), collapse = ";"),
      character(1)
    )
    
    ########################################
    # STEP 4.6: Flag high-confidence genes #
    ########################################
    
    # Criteria: posterior_prob > 0.5 AND regional_prob > 0.5 AND r2 > 0.8
    coloc_rows$high_conf <-
      (!is.na(coloc_rows$posterior_prob)   & coloc_rows$posterior_prob > 0.5) &
      (!is.na(coloc_rows$regional_prob)    & coloc_rows$regional_prob > 0.5) &
      (!is.na(coloc_rows$proxy_r2)         & coloc_rows$proxy_r2 > 0.8)
    
    ##############################################
    # STEP 5: Collapse to comma-separated columns#
    ##############################################
    
    coloc_summary <- coloc_rows |>
      dplyr::group_by(lead_snp) |>
      dplyr::summarise(
        !!col_pp        := paste(sprintf("%.4f", posterior_prob), collapse = ";"),
        !!col_reg       := paste(sprintf("%.4f", regional_prob), collapse = ";"),
        !!col_gene      := paste(unique(unlist(ensg_list)), collapse = ";"),
        !!col_proxy     := paste(unique(na.omit(candidate_snp)), collapse = ";"),
        !!col_proxy_pp  := paste(unique(na.omit(posterior_explained_by_snp)), collapse = ";"),
        !!col_proxy_r2  := paste(sprintf("%.4f", proxy_r2), collapse = ";"),
        !!col_sym       := paste(gene_symbol, collapse = ";"),
        !!col_high_conf := {
          hc <- gene_symbol[high_conf]
          if (length(hc) == 0 || all(is.na(hc))) NA_character_
          else paste(hc, collapse = ";")
        },
        .groups = "drop"
      )
    
    # empty strings -> NA
    coloc_summary[[col_pp]]        <- .empty_to_na(coloc_summary[[col_pp]])
    coloc_summary[[col_reg]]       <- .empty_to_na(coloc_summary[[col_reg]])
    coloc_summary[[col_gene]]      <- .empty_to_na(coloc_summary[[col_gene]])
    coloc_summary[[col_proxy]]     <- .empty_to_na(coloc_summary[[col_proxy]])
    coloc_summary[[col_proxy_pp]]  <- .empty_to_na(coloc_summary[[col_proxy_pp]])
    coloc_summary[[col_proxy_r2]]  <- .empty_to_na(coloc_summary[[col_proxy_r2]])
    
  } else {
    # No colocalization rows at all — create empty summary
    coloc_summary <- data.frame(lead_snp = character(0), stringsAsFactors = FALSE)
    coloc_summary[[col_pp]]         <- character(0)
    coloc_summary[[col_reg]]        <- character(0)
    coloc_summary[[col_gene]]       <- character(0)
    coloc_summary[[col_proxy]]      <- character(0)
    coloc_summary[[col_proxy_pp]]   <- character(0)
    coloc_summary[[col_proxy_r2]]   <- character(0)
    coloc_summary[[col_sym]]        <- character(0)
    coloc_summary[[col_high_conf]]  <- character(0)
  }
  
  #################################################
  # STEP 6: add the info to the original dataframe#
  #################################################
  
  join_by <- stats::setNames("lead_snp", ir_snp_col)
  
  if (nrow(coloc_summary) > 0) {
    ir_df <- ir_df |>
      dplyr::left_join(
        coloc_summary[, c("lead_snp", col_pp, col_reg, col_gene,
                          col_proxy, col_proxy_pp, col_proxy_r2,
                          col_sym, col_high_conf)],
        by = join_by
      )
  } else {
    ir_df[[col_pp]]         <- NA_character_
    ir_df[[col_reg]]        <- NA_character_
    ir_df[[col_gene]]       <- NA_character_
    ir_df[[col_proxy]]      <- NA_character_
    ir_df[[col_proxy_pp]]   <- NA_character_
    ir_df[[col_proxy_r2]]   <- NA_character_
    ir_df[[col_sym]]        <- NA_character_
    ir_df[[col_high_conf]]  <- NA_character_
  }
  
  return(ir_df)
}

add_coloc_to_irvariants_eqtls_cond <- function(ir_df, coloc_df, tissue_name, ensdb = NULL, ir_snp_col = "Lead IR SNP", token = "04cad4ca4374", pop = "EUR", ld_delay = 0.5) {
  
  # STEP 0: prepare the column names according to the tissue
  col_pp       <- paste0(tissue_name, "_colo_pp")
  col_reg      <- paste0(tissue_name, "_regional_pp")
  col_gene     <- paste0(tissue_name, "_coloc_gene")
  col_sym      <- paste0(tissue_name, "_coloc_symbol")
  col_proxy    <- paste0(tissue_name, "_proxy")
  col_proxy_pp <- paste0(tissue_name, "_proxy_pp")
  col_proxy_r2 <- paste0(tissue_name, "_proxy_r2")
  col_high_conf <- paste0(tissue_name, "_high_conf_gene")
  
  ####################################################################
  # STEP 1: make sure that we have the right columns in the coloc res#
  ####################################################################
  
  required_ir <- ir_snp_col
  required_colo <- c("lead_snp", "traits", "posterior_prob",
                     "regional_prob", "candidate_snp",
                     "posterior_explained_by_snp")
  
  if (!required_ir %in% names(ir_df))
    stop("ir_df must contain column '", required_ir, "'")
  missing_cols <- setdiff(required_colo, names(coloc_df))
  if (length(missing_cols) > 0)
    stop("coloc_df missing required columns: ",
         paste(missing_cols, collapse = ", "))
  
  #############################################################################
  # STEP 2: Deduplicate coloc results, though we have already done this before#
  #############################################################################
  
  coloc_df <- unique(coloc_df)
  
  #######################################
  # STEP 3: Extract ENSG IDs from traits#
  #######################################
  
  # traits looks like "iradjbmi, ENSG00000260063" or "None"
  coloc_df$ensg_list <- lapply(coloc_df$traits, function(x) {
    if (is.na(x)) return(character(0))
    unlist(regmatches(x, gregexpr("ENSG[0-9]{11}", x)))
  })
  coloc_df$has_gene <- lengths(coloc_df$ensg_list) > 0
  
  #################################################################
  # STEP 4: get co-localization results in a more structured shape#
  ##################################################################
  
  coloc_rows <- coloc_df[coloc_df$has_gene, ]
  
  if (nrow(coloc_rows) > 0) {
    
    #########################################################################
    # STEP 4.1: Compute LD r2 between lead_snp and candidate_snp via LDlinkR#
    #########################################################################
    
    if (!is.null(token)) {
      r2_cache <- new.env(hash = TRUE, parent = emptyenv())
      n_pairs <- nrow(coloc_rows)
      message("Computing LD r2 for ", n_pairs, " (lead, proxy) pair(s) via LDlinkR...")
      
      coloc_rows$proxy_r2 <- sapply(seq_len(n_pairs), function(i) {
        if (i %% 10 == 0) message("  ...pair ", i, " of ", n_pairs)
        compute_r2(
          lead   = coloc_rows$lead_snp[i],
          proxy  = coloc_rows$candidate_snp[i],
          token  = token,
          pop    = "EUR")
      })
      
    } else {
      warning("No LDlinkR token provided — proxy_r2 and high_conf_gene will be NA.",
              call. = FALSE)
      coloc_rows$proxy_r2 <- NA_real_
    }
    
    ########################################################################
    # STEP 4.2: Map ENSG -> gene symbol (batch, then per-row for high-conf)#
    ########################################################################
    
    all_ensgs <- unique(unlist(coloc_rows$ensg_list))
    
    #Let's transform them to symbols with EnsDb:
    
    if (is.null(ensdb)) {
      if (!requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE))
        stop("Install EnsDb.Hsapiens.v86: BiocManager::install('EnsDb.Hsapiens.v86')")
      ensdb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
    }
    
    #Let's make a map
    
    symbol_map <- AnnotationDbi::mapIds(
      ensdb,
      keys      = all_ensgs,
      column    = "SYMBOL",
      keytype   = "GENEID",
      multiVals = "first"
    )
    
    #Let's update coloc_rows so that we have the data fully mapped (not convering or pasting data per lead yet)
    coloc_rows$gene_symbol <- vapply(
      coloc_rows$ensg_list,
      function(ensgs) paste(unname(symbol_map[ensgs]), collapse = ";"),
      character(1)
    )
    
    ########################################
    # STEP 4.6: Flag high-confidence genes #
    ########################################
    
    # Criteria: posterior_prob > 0.5 AND regional_prob > 0.5 AND r2 > 0.8
    coloc_rows$high_conf <-
      (!is.na(coloc_rows$posterior_prob)   & coloc_rows$posterior_prob > 0.5) &
      (!is.na(coloc_rows$regional_prob)    & coloc_rows$regional_prob > 0.5) &
      (!is.na(coloc_rows$proxy_r2)         & coloc_rows$proxy_r2 > 0.8)
    
    ##############################################
    # STEP 5: Collapse to comma-separated columns#
    ##############################################
    
    #Let's add the genes + lead SNPs to avoid issues:
    
    coloc_rows$traits = unlist(str_remove_all(coloc_rows$traits, "iradjbmi, "))
    
    coloc_summary <- coloc_rows |>
      dplyr::group_by(lead_snp) |>
      dplyr::summarise(
        !!col_pp        := paste(sprintf("%.4f", posterior_prob), collapse = ";"),
        !!col_reg       := paste(sprintf("%.4f", regional_prob), collapse = ";"),
        !!col_gene      := paste(traits, collapse = ";"),
        !!col_proxy     := paste(unique(na.omit(candidate_snp)), collapse = ";"),
        !!col_proxy_pp  := paste(unique(na.omit(posterior_explained_by_snp)), collapse = ";"),
        !!col_proxy_r2  := paste(sprintf("%.4f", proxy_r2), collapse = ";"),
        !!col_sym       := paste(gene_symbol, collapse = ";"),
        !!col_high_conf := {
          hc <- gene_symbol[high_conf]
          if (length(hc) == 0 || all(is.na(hc))) NA_character_
          else paste(hc, collapse = ";")
        },
        .groups = "drop"
      )
    
    # empty strings -> NA
    coloc_summary[[col_pp]]        <- .empty_to_na(coloc_summary[[col_pp]])
    coloc_summary[[col_reg]]       <- .empty_to_na(coloc_summary[[col_reg]])
    coloc_summary[[col_gene]]      <- .empty_to_na(coloc_summary[[col_gene]])
    coloc_summary[[col_proxy]]     <- .empty_to_na(coloc_summary[[col_proxy]])
    coloc_summary[[col_proxy_pp]]  <- .empty_to_na(coloc_summary[[col_proxy_pp]])
    coloc_summary[[col_proxy_r2]]  <- .empty_to_na(coloc_summary[[col_proxy_r2]])
    
  } else {
    # No colocalization rows at all — create empty summary
    coloc_summary <- data.frame(lead_snp = character(0), stringsAsFactors = FALSE)
    coloc_summary[[col_pp]]         <- character(0)
    coloc_summary[[col_reg]]        <- character(0)
    coloc_summary[[col_gene]]       <- character(0)
    coloc_summary[[col_proxy]]      <- character(0)
    coloc_summary[[col_proxy_pp]]   <- character(0)
    coloc_summary[[col_proxy_r2]]   <- character(0)
    coloc_summary[[col_sym]]        <- character(0)
    coloc_summary[[col_high_conf]]  <- character(0)
  }
  
  #################################################
  # STEP 6: add the info to the original dataframe#
  #################################################
  
  join_by <- stats::setNames("lead_snp", ir_snp_col)
  
  if (nrow(coloc_summary) > 0) {
    ir_df <- ir_df |>
      dplyr::left_join(
        coloc_summary[, c("lead_snp", col_pp, col_reg, col_gene,
                          col_proxy, col_proxy_pp, col_proxy_r2,
                          col_sym, col_high_conf)],
        by = join_by
      )
  } else {
    ir_df[[col_pp]]         <- NA_character_
    ir_df[[col_reg]]        <- NA_character_
    ir_df[[col_gene]]       <- NA_character_
    ir_df[[col_proxy]]      <- NA_character_
    ir_df[[col_proxy_pp]]   <- NA_character_
    ir_df[[col_proxy_r2]]   <- NA_character_
    ir_df[[col_sym]]        <- NA_character_
    ir_df[[col_high_conf]]  <- NA_character_
  }
  
  return(ir_df)
}

add_coloc_to_irvariants_sqtls <- function(ir_df, coloc_df, tissue_name, ensdb = NULL, ir_snp_col = "Lead IR SNP", token = "04cad4ca4374", pop = "EUR", ld_delay = 0.5) {
  
  # STEP 0: prepare the column names according to the tissue
  col_pp       <- paste0(tissue_name, "_colo_pp")
  col_reg      <- paste0(tissue_name, "_regional_pp")
  col_gene     <- paste0(tissue_name, "_coloc_gene")
  col_sym      <- paste0(tissue_name, "_coloc_symbol")
  col_proxy    <- paste0(tissue_name, "_proxy")
  col_proxy_pp <- paste0(tissue_name, "_proxy_pp")
  col_proxy_r2 <- paste0(tissue_name, "_proxy_r2")
  col_high_conf <- paste0(tissue_name, "_high_conf_gene")
  
  ####################################################################
  # STEP 1: make sure that we have the right columns in the coloc res#
  ####################################################################
  
  required_ir <- ir_snp_col
  required_colo <- c("lead_snp", "traits", "posterior_prob",
                     "regional_prob", "candidate_snp",
                     "posterior_explained_by_snp")
  
  if (!required_ir %in% names(ir_df))
    stop("ir_df must contain column '", required_ir, "'")
  missing_cols <- setdiff(required_colo, names(coloc_df))
  if (length(missing_cols) > 0)
    stop("coloc_df missing required columns: ",
         paste(missing_cols, collapse = ", "))
  
  #############################################################################
  # STEP 2: Deduplicate coloc results, though we have already done this before#
  #############################################################################
  
  coloc_df <- unique(coloc_df)
  
  #######################################
  # STEP 3: Extract ENSG IDs from traits#
  #######################################
  
  # traits looks like "iradjbmi, ENSG00000260063" or "None"
  coloc_df$ensg_list <- lapply(coloc_df$traits, function(x) {
    if (is.na(x)) return(character(0))
    unlist(regmatches(x, gregexpr("ENSG[0-9]{11}", x)))
  })
  coloc_df$has_gene <- lengths(coloc_df$ensg_list) > 0
  
  #################################################################
  # STEP 4: get co-localization results in a more structured shape#
  ##################################################################
  
  coloc_rows <- coloc_df[coloc_df$has_gene, ]
  
  if (nrow(coloc_rows) > 0) {
    
    #########################################################################
    # STEP 4.1: Compute LD r2 between lead_snp and candidate_snp via LDlinkR#
    #########################################################################
    
    if (!is.null(token)) {
      r2_cache <- new.env(hash = TRUE, parent = emptyenv())
      n_pairs <- nrow(coloc_rows)
      message("Computing LD r2 for ", n_pairs, " (lead, proxy) pair(s) via LDlinkR...")
      
      coloc_rows$proxy_r2 <- sapply(seq_len(n_pairs), function(i) {
        if (i %% 10 == 0) message("  ...pair ", i, " of ", n_pairs)
        compute_r2(
          lead   = coloc_rows$lead_snp[i],
          proxy  = coloc_rows$candidate_snp[i],
          token  = token,
          pop    = "EUR")
      })
      
    } else {
      warning("No LDlinkR token provided — proxy_r2 and high_conf_gene will be NA.",
              call. = FALSE)
      coloc_rows$proxy_r2 <- NA_real_
    }
    
    ########################################################################
    # STEP 4.2: Map ENSG -> gene symbol (batch, then per-row for high-conf)#
    ########################################################################
    
    all_ensgs <- unique(unlist(coloc_rows$ensg_list))
    
    #Let's transform them to symbols with EnsDb:
    
    if (is.null(ensdb)) {
      if (!requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE))
        stop("Install EnsDb.Hsapiens.v86: BiocManager::install('EnsDb.Hsapiens.v86')")
      ensdb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
    }
    
    #Let's make a map
    
    symbol_map <- AnnotationDbi::mapIds(
      ensdb,
      keys      = all_ensgs,
      column    = "SYMBOL",
      keytype   = "GENEID",
      multiVals = "first"
    )
    
    #Let's update coloc_rows so that we have the data fully mapped (not convering or pasting data per lead yet)
    coloc_rows$gene_symbol <- vapply(
      coloc_rows$ensg_list,
      function(ensgs) paste(unname(symbol_map[ensgs]), collapse = ";"),
      character(1)
    )
    
    ########################################
    # STEP 4.6: Flag high-confidence genes #
    ########################################
    
    # Criteria: posterior_prob > 0.5 AND regional_prob > 0.5 AND r2 > 0.8
    coloc_rows$high_conf <-
      (!is.na(coloc_rows$posterior_prob)   & coloc_rows$posterior_prob > 0.5) &
      (!is.na(coloc_rows$regional_prob)    & coloc_rows$regional_prob > 0.5) &
      (!is.na(coloc_rows$proxy_r2)         & coloc_rows$proxy_r2 > 0.8)
    
    ##############################################
    # STEP 5: Collapse to comma-separated columns#
    ##############################################
    
    #Let's add the genes + lead SNPs to avoid issues:
    
    coloc_rows$traits = unlist(str_remove_all(coloc_rows$traits, "iradjbmi, "))
    
    coloc_summary <- coloc_rows |>
      dplyr::group_by(lead_snp) |>
      dplyr::summarise(
        !!col_pp        := paste(sprintf("%.4f", posterior_prob), collapse = ";"),
        !!col_reg       := paste(sprintf("%.4f", regional_prob), collapse = ";"),
        !!col_gene      := paste(traits, collapse = ";"),
        !!col_proxy     := paste(unique(na.omit(candidate_snp)), collapse = ";"),
        !!col_proxy_pp  := paste(unique(na.omit(posterior_explained_by_snp)), collapse = ";"),
        !!col_proxy_r2  := paste(sprintf("%.4f", proxy_r2), collapse = ";"),
        !!col_sym       := paste(gene_symbol, collapse = ";"),
        !!col_high_conf := {
          hc <- gene_symbol[high_conf]
          if (length(hc) == 0 || all(is.na(hc))) NA_character_
          else paste(hc, collapse = ";")
        },
        .groups = "drop"
      )
    
    # empty strings -> NA
    coloc_summary[[col_pp]]        <- .empty_to_na(coloc_summary[[col_pp]])
    coloc_summary[[col_reg]]       <- .empty_to_na(coloc_summary[[col_reg]])
    coloc_summary[[col_gene]]      <- .empty_to_na(coloc_summary[[col_gene]])
    coloc_summary[[col_proxy]]     <- .empty_to_na(coloc_summary[[col_proxy]])
    coloc_summary[[col_proxy_pp]]  <- .empty_to_na(coloc_summary[[col_proxy_pp]])
    coloc_summary[[col_proxy_r2]]  <- .empty_to_na(coloc_summary[[col_proxy_r2]])
    
  } else {
    # No colocalization rows at all — create empty summary
    coloc_summary <- data.frame(lead_snp = character(0), stringsAsFactors = FALSE)
    coloc_summary[[col_pp]]         <- character(0)
    coloc_summary[[col_reg]]        <- character(0)
    coloc_summary[[col_gene]]       <- character(0)
    coloc_summary[[col_proxy]]      <- character(0)
    coloc_summary[[col_proxy_pp]]   <- character(0)
    coloc_summary[[col_proxy_r2]]   <- character(0)
    coloc_summary[[col_sym]]        <- character(0)
    coloc_summary[[col_high_conf]]  <- character(0)
  }
  
  #################################################
  # STEP 6: add the info to the original dataframe#
  #################################################
  
  join_by <- stats::setNames("lead_snp", ir_snp_col)
  
  if (nrow(coloc_summary) > 0) {
    ir_df <- ir_df |>
      dplyr::left_join(
        coloc_summary[, c("lead_snp", col_pp, col_reg, col_gene,
                          col_proxy, col_proxy_pp, col_proxy_r2,
                          col_sym, col_high_conf)],
        by = join_by
      )
  } else {
    ir_df[[col_pp]]         <- NA_character_
    ir_df[[col_reg]]        <- NA_character_
    ir_df[[col_gene]]       <- NA_character_
    ir_df[[col_proxy]]      <- NA_character_
    ir_df[[col_proxy_pp]]   <- NA_character_
    ir_df[[col_proxy_r2]]   <- NA_character_
    ir_df[[col_sym]]        <- NA_character_
    ir_df[[col_high_conf]]  <- NA_character_
  }
  
  return(ir_df)
}


add_coloc_to_irvariants_pqtls <- function(ir_df, coloc_df, tissue_name, ensdb = NULL, ir_snp_col = "Lead IR SNP", token = "04cad4ca4374", pop = "EUR", ld_delay = 0.5) {
  
  # STEP 0: prepare the column names according to the tissue
  col_pp       <- paste0(tissue_name, "_colo_pp")
  col_reg      <- paste0(tissue_name, "_regional_pp")
  col_gene     <- paste0(tissue_name, "_coloc_gene")
  col_sym      <- paste0(tissue_name, "_coloc_symbol")
  col_proxy    <- paste0(tissue_name, "_proxy")
  col_proxy_pp <- paste0(tissue_name, "_proxy_pp")
  col_proxy_r2 <- paste0(tissue_name, "_proxy_r2")
  col_high_conf <- paste0(tissue_name, "_high_conf_gene")
  
  # STEP 1: make sure that we have the right columns in the coloc res
  required_ir <- ir_snp_col
  required_colo <- c("lead_snp", "traits", "posterior_prob",
                     "regional_prob", "candidate_snp",
                     "posterior_explained_by_snp")
  
  if (!required_ir %in% names(ir_df))
    stop("ir_df must contain column '", required_ir, "'")
  missing_cols <- setdiff(required_colo, names(coloc_df))
  if (length(missing_cols) > 0)
    stop("coloc_df missing required columns: ",
         paste(missing_cols, collapse = ", "))
  
  # STEP 2: Deduplicate coloc results, though we have already done this before
  coloc_df <- unique(coloc_df)
  
  # STEP 3: Extract ENSG IDs from traits
  # traits looks like "iradjbmi, ENSG00000260063" or "None"
  coloc_df$gene_symbol <- ifelse(
    is.na(coloc_df$traits) | coloc_df$traits == "None",
    NA_character_,
    str_trim(str_remove(coloc_df$traits, "^iradjbmi,\\s*"))
  )
  coloc_df$has_gene <- ifelse(is.na(coloc_df$gene_symbol), FALSE, TRUE)
  
  #################################################################
  # STEP 4: get co-localization results in a more structured shape#
  ##################################################################
  
  coloc_rows <- coloc_df[coloc_df$has_gene, ]
  
  if (nrow(coloc_rows) > 0) {
    
    #########################################################################
    # STEP 4.1: Compute LD r2 between lead_snp and candidate_snp via LDlinkR#
    #########################################################################
    
    if (!is.null(token)) {
      r2_cache <- new.env(hash = TRUE, parent = emptyenv())
      n_pairs <- nrow(coloc_rows)
      message("Computing LD r2 for ", n_pairs, " (lead, proxy) pair(s) via LDlinkR...")
      
      coloc_rows$proxy_r2 <- sapply(seq_len(n_pairs), function(i) {
        if (i %% 10 == 0) message("  ...pair ", i, " of ", n_pairs)
        compute_r2(
          lead   = coloc_rows$lead_snp[i],
          proxy  = coloc_rows$candidate_snp[i],
          token  = token,
          pop    = "EUR")
      })
      
    } else {
      warning("No LDlinkR token provided — proxy_r2 and high_conf_gene will be NA.",
              call. = FALSE)
      coloc_rows$proxy_r2 <- NA_real_
    }
    
    ########################################
    # STEP 4.6: Flag high-confidence genes #
    ########################################
    
    # Criteria: posterior_prob > 0.5 AND regional_prob > 0.5 AND r2 > 0.8
    coloc_rows$high_conf <-
      (!is.na(coloc_rows$posterior_prob)   & coloc_rows$posterior_prob > 0.5) &
      (!is.na(coloc_rows$regional_prob)    & coloc_rows$regional_prob > 0.5) &
      (!is.na(coloc_rows$proxy_r2)         & coloc_rows$proxy_r2 > 0.8)
    
    ##############################################
    # STEP 5: Collapse to comma-separated columns#
    ##############################################
    
    coloc_summary <- coloc_rows |>
      dplyr::group_by(lead_snp) |>
      dplyr::summarise(
        !!col_pp        := paste(sprintf("%.4f", posterior_prob), collapse = ";"),
        !!col_reg       := paste(sprintf("%.4f", regional_prob), collapse = ":"),
        # !!col_gene      := paste(unique(unlist(ensg_list)), collapse = ";"),
        !!col_proxy     := paste(unique(na.omit(candidate_snp)), collapse = ";"),
        !!col_proxy_pp  := paste(unique(na.omit(posterior_explained_by_snp)), collapse = ";"),
        !!col_proxy_r2  := paste(sprintf("%.4f", proxy_r2), collapse = ";"),
        !!col_sym       := paste(gene_symbol, collapse = ";"),
        !!col_high_conf := {
          hc <- gene_symbol[high_conf]
          if (length(hc) == 0 || all(is.na(hc))) NA_character_
          else paste(hc, collapse = ";")
        },
        .groups = "drop"
      )
    
    # empty strings -> NA
    coloc_summary[[col_pp]]        <- .empty_to_na(coloc_summary[[col_pp]])
    coloc_summary[[col_reg]]       <- .empty_to_na(coloc_summary[[col_reg]])
    # coloc_summary[[col_gene]]      <- .empty_to_na(coloc_summary[[col_gene]])
    coloc_summary[[col_proxy]]     <- .empty_to_na(coloc_summary[[col_proxy]])
    coloc_summary[[col_proxy_pp]]  <- .empty_to_na(coloc_summary[[col_proxy_pp]])
    coloc_summary[[col_proxy_r2]]  <- .empty_to_na(coloc_summary[[col_proxy_r2]])
    
  } else {
    # No colocalization rows at all — create empty summary
    coloc_summary <- data.frame(lead_snp = character(0), stringsAsFactors = FALSE)
    coloc_summary[[col_pp]]         <- character(0)
    coloc_summary[[col_reg]]        <- character(0)
    coloc_summary[[col_gene]]       <- character(0)
    coloc_summary[[col_proxy]]      <- character(0)
    coloc_summary[[col_proxy_pp]]   <- character(0)
    coloc_summary[[col_proxy_r2]]   <- character(0)
    coloc_summary[[col_sym]]        <- character(0)
    coloc_summary[[col_high_conf]]  <- character(0)
  }
  
  #################################################
  # STEP 6: add the info to the original dataframe#
  #################################################
  
  join_by <- stats::setNames("lead_snp", ir_snp_col)
  
  if (nrow(coloc_summary) > 0) {
    ir_df <- ir_df |>
      dplyr::left_join(
        coloc_summary[, c("lead_snp", col_pp, col_reg, 
                          col_proxy, col_proxy_pp, col_proxy_r2,
                          col_sym, col_high_conf)],
        by = join_by
      )
  } else {
    ir_df[[col_pp]]         <- NA_character_
    ir_df[[col_reg]]        <- NA_character_
    # ir_df[[col_gene]]       <- NA_character_
    ir_df[[col_proxy]]      <- NA_character_
    ir_df[[col_proxy_pp]]   <- NA_character_
    ir_df[[col_proxy_r2]]   <- NA_character_
    ir_df[[col_sym]]        <- NA_character_
    ir_df[[col_high_conf]]  <- NA_character_
  }
  
  return(ir_df)
}


##############
#Loading data#
##############

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/review/")

ir_variants = readxl::read_xlsx("manuscript/10082026/Supplementary Tables_11August2026.xlsx", sheet = 3)
ir_variants = as.data.frame(ir_variants)

#Let's first select the columns that we are interested in and later parse them a bit better to make the data easy to follow

ir_variants_ = ir_variants[,seq(1,12),]
colnames(ir_variants_) = ir_variants_[2,]
ir_variants_ = ir_variants_[3:284,]

#And let's just parse the info one second:  

ir_variants_$`Novel/Reported`[which(is.na(ir_variants_$`Novel/Reported`))] = "Novel"

#Let's load now the proxies!!

proxies <- fread("../output/4_ir_loci_discovery/4_proxies/proxies_build_37_and_38_4_282_variants.txt")
proxies <- proxies[which(proxies$query_snp_rsid%in%ir_variants_$`Lead IR SNP`),]

length(which(duplicated(proxies$query_snp_rsid) == FALSE)) #282/282

#Let's add chr_pos just in case:

proxies$chr_pos <- paste("chr", proxies$chr, ":", proxies$pos_hg19, sep= "") #let's add the chr_pos 

###########################################################################
#Alright, let's add the raw data sets and see how we can compile this info#
###########################################################################

asat_adipo = readRDS("review_output/4_colocalization_and_gtex_checks/raw_asat_adipoxpress_coloc_df.RDS")
asat_adipo_cond = readRDS("review_output/4_colocalization_and_gtex_checks/raw_asat_adipoxpress_cond_coloc_df.RDS")

asat_gtex = readRDS("review_output/4_colocalization_and_gtex_checks/raw_asat_gtex_coloc_df.RDS")
asat_gtex_sqtl = readRDS("review_output/4_colocalization_and_gtex_checks/raw_asat_gtex_sqtl_coloc_df.RDS")

vat_gtex = readRDS("review_output/4_colocalization_and_gtex_checks/raw_vat_gtex_coloc_df.RDS")
vat_gtex_sqtl = readRDS("review_output/4_colocalization_and_gtex_checks/raw_vat_gtex_sqtl_coloc_df.RDS")

muscle_gtex = readRDS("review_output/4_colocalization_and_gtex_checks/raw_muscle_gtex_coloc_df.RDS")
muscle_gtex_sqtl = readRDS("review_output/4_colocalization_and_gtex_checks/raw_muscle_gtex_sqtl_coloc_df.RDS")

liver_gtex = readRDS("review_output/4_colocalization_and_gtex_checks/raw_liver_gtex_coloc_df.RDS")
liver_gtex_sqtl = readRDS("review_output/4_colocalization_and_gtex_checks/raw_liver_gtex_sqtl_coloc_df.RDS")

pqtl = readRDS("review_output/4_colocalization_and_gtex_checks/raw_pqtl_coloc_df.RDS")

##################################################################################################################################################
#Let's clean the datasets for the time being - we are gonna need different versions of the same function depending on the structure and data type#
##################################################################################################################################################

#1. Let's first do this for ASAT cis-eQTLs for AdipoXpress meta-analyses

#Removing loci that did not co-localize
asat_adipo = asat_adipo[which(asat_adipo$traits != "None" & is.na(asat_adipo$traits) == FALSE),]
#Removing duplicates that may happen in busy locus
asat_adipo = asat_adipo[which(duplicated(asat_adipo) == FALSE),] #sometimes hits that are very close end up duplicating co-localization results. Hence why we need to do this
#Using our function to add the data to the IR dataframe:
ir_variants_ = add_coloc_to_irvariants_eqtls(ir_df=ir_variants_, coloc_df = asat_adipo, tissue_name="asat_adipoxpress_marginal_eqtl")

#2. Now the conditional co-localization results from Adipoxpress too:

#Removing loci that did not co-localize
asat_adipo_cond = asat_adipo_cond[which(asat_adipo_cond$traits != "None" & is.na(asat_adipo_cond$traits) == FALSE),]
#Removing duplicates that may happen in busy locus
asat_adipo_cond = asat_adipo_cond[which(duplicated(asat_adipo_cond) == FALSE),] #sometimes hits that are very close end up duplicating co-localization results. Hence why we need to do this
#Using our function to add the data to the IR dataframe:
ir_variants_ = add_coloc_to_irvariants_eqtls_cond(ir_df=ir_variants_, coloc_df = asat_adipo_cond, tissue_name="asat_adipoxpress_conditional_eqtl")

#3. ASAT from cis-eQTL from GTEx:

#Removing loci that did not co-localize
asat_gtex = asat_gtex[which(asat_gtex$traits != "None" & is.na(asat_gtex$traits) == FALSE),]
#Removing duplicates that may happen in busy locus
asat_gtex = asat_gtex[which(duplicated(asat_gtex) == FALSE),] #sometimes hits that are very close end up duplicating co-localization results. Hence why we need to do this
#Using our function to add the data to the IR dataframe:
ir_variants_ = add_coloc_to_irvariants_eqtls(ir_df=ir_variants_, coloc_df = asat_gtex, tissue_name="asat_gtex_marginal_eqtl")

#4. ASAT from cis-sQTL from GTEx:

#Removing loci that did not co-localize
asat_gtex_sqtl = asat_gtex_sqtl[which(asat_gtex_sqtl$traits != "None" & is.na(asat_gtex_sqtl$traits) == FALSE),]
#Removing duplicates that may happen in busy locus
asat_gtex_sqtl = asat_gtex_sqtl[which(duplicated(asat_gtex_sqtl) == FALSE),] #sometimes hits that are very close end up duplicating co-localization results. Hence why we need to do this
#Using our function to add the data to the IR dataframe:
ir_variants_ = add_coloc_to_irvariants_sqtls(ir_df=ir_variants_, coloc_df = asat_gtex_sqtl, tissue_name="asat_gtex_marginal_sqtl")

#5. Now vat from cis-eQTL from GTEx:

#Removing loci that did not co-localize
vat_gtex = vat_gtex[which(vat_gtex$traits != "None" & is.na(vat_gtex$traits) == FALSE),]
#Removing duplicates that may happen in busy locus
vat_gtex = vat_gtex[which(duplicated(vat_gtex) == FALSE),] #sometimes hits that are very close end up duplicating co-localization results. Hence why we need to do this
#Using our function to add the data to the IR dataframe:
ir_variants_ = add_coloc_to_irvariants_eqtls(ir_df=ir_variants_, coloc_df = vat_gtex, tissue_name="vat_gtex_marginal_eqtl")

#6. Now vat cis-sQTL from GTEx

#Removing loci that did not co-localize
vat_gtex_sqtl = vat_gtex_sqtl[which(vat_gtex_sqtl$traits != "None" & is.na(vat_gtex_sqtl$traits) == FALSE),]
#Removing duplicates that may happen in busy locus
vat_gtex_sqtl = vat_gtex_sqtl[which(duplicated(vat_gtex_sqtl) == FALSE),] #sometimes hits that are very close end up duplicating co-localization results. Hence why we need to do this
#Using our function to add the data to the IR dataframe:
ir_variants_ = add_coloc_to_irvariants_sqtls(ir_df=ir_variants_, coloc_df = vat_gtex_sqtl, tissue_name="vat_gtex_marginal_sqtl")

#7. Now Liver from cis-eQTL from GTEx:

#Removing loci that did not co-localize
liver_gtex = liver_gtex[which(liver_gtex$traits != "None" & is.na(liver_gtex$traits) == FALSE),]
#Removing duplicates that may happen in busy locus
liver_gtex = liver_gtex[which(duplicated(liver_gtex) == FALSE),] #sometimes hits that are very close end up duplicating co-localization results. Hence why we need to do this
#Using our function to add the data to the IR dataframe:
ir_variants_ = add_coloc_to_irvariants_eqtls(ir_df=ir_variants_, coloc_df = liver_gtex, tissue_name="liver_gtex_marginal_eqtl")

#8. Now vat cis-sQTL from GTEx

#Removing loci that did not co-localize
liver_gtex_sqtl = liver_gtex_sqtl[which(liver_gtex_sqtl$traits != "None" & is.na(liver_gtex_sqtl$traits) == FALSE),]
#Removing duplicates that may happen in busy locus
liver_gtex_sqtl = liver_gtex_sqtl[which(duplicated(liver_gtex_sqtl) == FALSE),] #sometimes hits that are very close end up duplicating co-localization results. Hence why we need to do this
#Using our function to add the data to the IR dataframe:
ir_variants_ = add_coloc_to_irvariants_sqtls(ir_df=ir_variants_, coloc_df = liver_gtex_sqtl, tissue_name="liver_gtex_marginal_sqtl")

#9. Now Muscle from cis-eQTL from GTEx:

#Removing loci that did not co-localize
muscle_gtex = muscle_gtex[which(muscle_gtex$traits != "None" & is.na(muscle_gtex$traits) == FALSE),]
#Removing duplicates that may happen in busy locus
muscle_gtex = muscle_gtex[which(duplicated(muscle_gtex) == FALSE),] #sometimes hits that are very close end up duplicating co-localization results. Hence why we need to do this
#Using our function to add the data to the IR dataframe:
ir_variants_ = add_coloc_to_irvariants_eqtls(ir_df=ir_variants_, coloc_df = muscle_gtex, tissue_name="muscle_gtex_marginal_eqtl")

#10. Now Muscle cis-sQTL from GTEx

#Removing loci that did not co-localize
muscle_gtex_sqtl = muscle_gtex_sqtl[which(muscle_gtex_sqtl$traits != "None" & is.na(muscle_gtex_sqtl$traits) == FALSE),]
#Removing duplicates that may happen in busy locus
muscle_gtex_sqtl = muscle_gtex_sqtl[which(duplicated(muscle_gtex_sqtl) == FALSE),] #sometimes hits that are very close end up duplicating co-localization results. Hence why we need to do this
#Using our function to add the data to the IR dataframe:
ir_variants_ = add_coloc_to_irvariants_sqtls(ir_df=ir_variants_, coloc_df = muscle_gtex_sqtl, tissue_name="muscle_gtex_marginal_sqtl")

#Now the pQTLs too:

#Removing loci that did not co-localize
pqtl = pqtl[which(pqtl$traits != "None" & is.na(pqtl$traits) == FALSE),]
#Removing duplicates that may happen in busy locus
pqtl = pqtl[which(duplicated(pqtl) == FALSE),] #sometimes hits that are very close end up duplicating co-localization results. Hence why we need to do this
#Using our function to add the data to the IR dataframe:
ir_variants_ = add_coloc_to_irvariants_pqtls(ir_df=ir_variants_, coloc_df = pqtl, tissue_name="ukbb_pqtl_marginal")

############################################################################
#Let's clean some of the columns post-hoc to be able to have a decent table#
############################################################################

cleaning_repeats = function(high_conf_vect, id_vect){
  
  final_vect= c()
  
  for(index in seq(1, length(high_conf_vect))){
    
    #STEP 0: make dummy examples:
    #high_conf_vect =  ir_variants_$asat_adipoxpress_marginal_eqtl_high_conf_gene[277]
    #id_vect =  ir_variants_$asat_adipoxpress_marginal_eqtl_coloc_gene[225]
    
    #STEP 1: let's skip if NA
    
    if(is.na(high_conf_vect[index])){
      
      final_vect <- c(final_vect, high_conf_vect[index])
      
    } else {
      
      #STEP 1: split both vectors to have the exact info in different indexable positions:
      high_conf_vect_ = unlist(str_split(high_conf_vect[index], "[;]"))
      id_vect_ = unlist(str_split(id_vect[index], "[;]"))
      
      #STEP 2: find the NA indexes and replace them:
      
      index_na= which(high_conf_vect_ == "NA")
      high_conf_vect_[index_na] = id_vect_[index_na]
      
      #STEP 3: remove redundancy, if any. Likely in conditional and sQTLs
      high_conf_vect_ = unique(high_conf_vect_)
      high_conf_vect_ = high_conf_vect_[order(high_conf_vect_, decreasing = FALSE)]
      
      #STEP 4: collapse and return
      high_conf_vect_clean=paste(high_conf_vect_, collapse = ";")
      
      final_vect <- c(final_vect, high_conf_vect_clean)
      
    }
    
  } #for loop
  
  return(final_vect)
  
} #function 

clean_data_asat_adipoxpress_marginal= cleaning_repeats(high_conf_vect=ir_variants_$asat_adipoxpress_marginal_eqtl_high_conf_gene, id_vect=ir_variants_$asat_adipoxpress_marginal_eqtl_coloc_gene)
clean_data_asat_adipoxpress_conditional= cleaning_repeats(high_conf_vect=ir_variants_$asat_adipoxpress_conditional_eqtl_high_conf_gene, id_vect=ir_variants_$asat_adipoxpress_conditional_eqtl_coloc_gene)

clean_data_asat_gtex_eqtl= cleaning_repeats(high_conf_vect=ir_variants_$asat_gtex_marginal_eqtl_high_conf_gene, id_vect=ir_variants_$asat_gtex_marginal_eqtl_coloc_gene)
clean_data_asat_gtex_sqtl= cleaning_repeats(high_conf_vect=ir_variants_$asat_gtex_marginal_sqtl_high_conf_gene, id_vect=ir_variants_$asat_gtex_marginal_sqtl_coloc_gene)

clean_data_vat_gtex_eqtl= cleaning_repeats(high_conf_vect=ir_variants_$vat_gtex_marginal_eqtl_high_conf_gene, id_vect=ir_variants_$vat_gtex_marginal_eqtl_coloc_gene)
clean_data_vat_gtex_sqtl= cleaning_repeats(high_conf_vect=ir_variants_$vat_gtex_marginal_sqtl_high_conf_gene, id_vect=ir_variants_$vat_gtex_marginal_sqtl_coloc_gene)

clean_data_liver_gtex_eqtl= cleaning_repeats(high_conf_vect=ir_variants_$liver_gtex_marginal_eqtl_high_conf_gene, id_vect=ir_variants_$liver_gtex_marginal_eqtl_coloc_gene)
clean_data_liver_gtex_sqtl= cleaning_repeats(high_conf_vect=ir_variants_$liver_gtex_marginal_sqtl_high_conf_gene, id_vect=ir_variants_$liver_gtex_marginal_sqtl_coloc_gene)

clean_data_muscle_gtex_eqtl= cleaning_repeats(high_conf_vect=ir_variants_$muscle_gtex_marginal_eqtl_high_conf_gene, id_vect=ir_variants_$muscle_gtex_marginal_eqtl_coloc_gene)
clean_data_muscle_gtex_sqtl= cleaning_repeats(high_conf_vect=ir_variants_$muscle_gtex_marginal_sqtl_high_conf_gene, id_vect=ir_variants_$muscle_gtex_marginal_sqtl_coloc_gene)

#############################################################################
#Alright, let's update the data and get it ready to become a freaking table!#
#############################################################################

ir_variants_$asat_adipoxpress_marginal_eqtl_high_conf_gene=clean_data_asat_adipoxpress_marginal
ir_variants_$asat_adipoxpress_conditional_eqtl_high_conf_gene=clean_data_asat_adipoxpress_conditional

ir_variants_$asat_gtex_marginal_eqtl_high_conf_gene=clean_data_asat_gtex_eqtl
ir_variants_$asat_gtex_marginal_sqtl_high_conf_gene=clean_data_asat_gtex_sqtl

ir_variants_$vat_gtex_marginal_eqtl_high_conf_gene=clean_data_vat_gtex_eqtl
ir_variants_$vat_gtex_marginal_sqtl_high_conf_gene=clean_data_vat_gtex_sqtl

ir_variants_$liver_gtex_marginal_eqtl_high_conf_gene=clean_data_liver_gtex_eqtl
ir_variants_$liver_gtex_marginal_sqtl_high_conf_gene=clean_data_liver_gtex_sqtl

ir_variants_$muscle_gtex_marginal_eqtl_high_conf_gene=clean_data_muscle_gtex_eqtl
ir_variants_$muscle_gtex_marginal_sqtl_high_conf_gene=clean_data_muscle_gtex_sqtl

fwrite(ir_variants_, "review_output/4_colocalization_and_gtex_checks/parsed_results_all_colocs.txt")

##############################################################################################
#Alright!!! We have the co-localization results, let's go and try to get the comparisons done#
##############################################################################################

ir_variants_ = fread("review_output/4_colocalization_and_gtex_checks/parsed_results_all_colocs.txt")

#Get the genes involved:

coloc_sat = ir_variants_[which(is.na(ir_variants_$asat_adipoxpress_marginal_eqtl_high_conf_gene) == FALSE |
                               is.na(ir_variants_$asat_adipoxpress_conditional_eqtl_high_conf_gene) == FALSE |
                               is.na(ir_variants_$asat_gtex_marginal_eqtl_high_conf_gene) == FALSE |
                               is.na(ir_variants_$asat_gtex_marginal_sqtl_high_conf_gene) == FALSE),]

unique_genes = c(unlist(str_split(ir_variants_$asat_adipoxpress_marginal_eqtl_high_conf_gene, ";")),
                 unlist(str_split(ir_variants_$asat_adipoxpress_conditional_eqtl_high_conf_gene, ";")),
                 unlist(str_split(ir_variants_$asat_gtex_marginal_eqtl_high_conf_gene, ";")),
                 unlist(str_split(ir_variants_$asat_gtex_marginal_sqtl_high_conf_gene, ";")))

unique_genes=unique_genes[which(is.na(unique_genes) == FALSE)]
unique_genes=unique(unique_genes)

###############
#Now the pQTLs#
###############

pqtl = ir_variants_[which(is.na(ir_variants_$ukbb_pqtl_marginal_high_conf_gene) == FALSE),]
