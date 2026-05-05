################################################################################
##  Differential Abundance Pipeline
##  Author: Camila Alvarez-Silva
##  Methods: LinDA | ANCOMBC2 | LM | edgeR | metagenomeSeq
##  Covariates: Age + BMI (extend via `covariates` argument)
##  Input : named list of phyloseq objects
################################################################################

# ── 0. Dependencies ────────────────────────────────────────────────────────────
required_pkgs <- c(
  "phyloseq", "MicrobiomeStat", "ANCOMBC", "edgeR",
  "metagenomeSeq", "lme4", "lmerTest", "broom.mixed",
  "tidyverse", "ggplot2", "ggrepel", "patchwork")

invisible(lapply(required_pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE))
    stop(paste("Package not installed:", p,
               "\nInstall via BiocManager::install() or install.packages()"))
  library(p, character.only = TRUE)
}))


# ── 1. Helper utilities ────────────────────────────────────────────────────────


#' Safe tax_table accessor — returns a data.frame always
safe_tax <- function(ps) {
  if (!is.null(tax_table(ps))) {
    as.data.frame(tax_table(ps))
  } else {
    data.frame(OTU = taxa_names(ps), row.names = taxa_names(ps))
  }
}

#' Standardise output to a common schema
#'   feature | log2FC | pval | padj | diff_robust | method | dataset
standardise <- function(df, method, dataset,
                        feature_col, log2fc_col, pval_col, padj_col, diff_robust) {
  df <- as.data.frame(df)
  out <- data.frame(
    feature  = df[[feature_col]],
    log2FC   = df[[log2fc_col]],
    pval     = df[[pval_col]],
    padj     = df[[padj_col]],
    diff_robust = df[[diff_robust]],
    method   = method,
    dataset  = dataset,
    stringsAsFactors = FALSE
  )
  out[!is.na(out$feature), ]
}


# ── 2. Individual DA wrappers ──────────────────────────────────────────────────

## 2a. LinDA ───────────────────────────────────────────────────────────────────
run_linda <- function(ps, formula, dataset_name,
                      p_adj_method = "fdr", alpha = 0.05) {
  
  message("[LinDA] ", dataset_name)
  
  # transform taxonomy phyloseq object from counts to relative abundances
  if(dataset_name== "Taxonomy"){ps <- transform_sample_counts(ps, function(x) x / sum(x))}
 
   tryCatch({
     
     meta <- as(sample_data(ps), "data.frame")
     rownames(meta)<-as.character(meta$DNA_ID) 
     otu <- as(otu_table(ps), "matrix")  # features x samples 
     
     otu <- otu[,rownames(meta)]  #data frame or matrix representing observed OTU table. Row: taxa; column: samples
     
     res<-LinDA::linda(otu.tab = (otu),
                                meta=meta,
                                alpha = alpha,
                                formula = formula  ,
                                p.adj.method  = p_adj_method,
                                type = "proportion")
     
    # Extract the first variable of interest (first term after ~)
    term <- names(res$output)[1]
    df   <- res$output[[term]]
    df$feature <- rownames(df)
    df$diff_robust <- NA
    
    standardise(df, "LinDA", dataset_name,
                "feature", "log2FoldChange", "pvalue", "padj","diff_robust")
  }, error = function(e) {
    warning("[LinDA] Failed on ", dataset_name, ": ", e$message)
    NULL
  })
}

## 2b. ANCOMBC2 ────────────────────────────────────────────────────────────────
run_ancombc2 <- function(ps, fixed_formula, dataset_name,
                         p_adj_method = "fdr", alpha = 0.05,
                         rand_formula = NULL) {
  
  message("[ANCOMBC2] ", dataset_name)
  tryCatch({
    res <- ancombc2(
      data         = ps,
      fix_formula  = fixed_formula,
      rand_formula = rand_formula,
      p_adj_method = p_adj_method,
      alpha        = alpha,
      prv_cut      = 0.10,
      n_cl         = 1
    )
    df   <- res$res
    # Identify the first lfc_ column
    lfc_col  <- grep("^lfc_",  names(df), value = TRUE)[2]
    pval_col <- grep("^p_",    names(df), value = TRUE)[2]
    padj_col <- grep("^q_",    names(df), value = TRUE)[2]
    diff_robust<- grep("^diff_robust",    names(df), value = TRUE)[2]
    data.frame(
      feature = df$taxon,
      log2FC  = df[[lfc_col]],
      pval    = df[[pval_col]],
      padj    = df[[padj_col]],
      diff_robust = df[[diff_robust]],
      method  = "ANCOMBC2",
      dataset = dataset_name,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    warning("[ANCOMBC2] Failed on ", dataset_name, ": ", e$message)
    NULL
  })
}

## 2c. Linear  Model (lmerTest) ──────────────────────────────────────────
run_lmm <- function(ps, variable_col, covariates = covariates,
                    random_effect = NULL, dataset_name,
                    p_adj_method = "fdr", alpha = 0.05) {
  message("[LM] ", dataset_name)
  tryCatch({
    # CLR-transform counts
    otu  <- as.data.frame(otu_table(ps))
    if (!taxa_are_rows(ps)) otu <- t(otu)
    otu  <- otu + min(otu[otu != 0])/2                         # pseudo-count
    clr  <- t(apply(otu, 2, function(x) log(x) - mean(log(x))))  # samples x taxa
    
    sdata <- as.data.frame(sample_data(ps))
    fixed <- paste(c(variable_col, covariates), collapse = " + ")
    formula_str <- if (!is.null(random_effect)) {
      paste0("taxa_val ~ ", fixed, " + (1|", random_effect, ")")
    } else {
      paste0("taxa_val ~ ", fixed)
    }
    
    results <- lapply(colnames(clr), function(tax) {
      
      df_tmp <- cbind(sdata, taxa_val = clr[, tax])
      fit <- tryCatch(
        if (!is.null(random_effect)) {
          lmerTest::lmer(as.formula(formula_str), data = df_tmp, REML = FALSE)
        } else {
          lm(as.formula(formula_str), data = df_tmp)
        },
        error = function(e) NULL
      )
      
      if (is.null(fit)) return(NULL)
      coef_tbl <- broom.mixed::tidy(fit, effects = "fixed")
      row <- coef_tbl[coef_tbl$term %like% variable_col, ]
      if (nrow(row) == 0) return(NULL)
      data.frame(feature = tax,
                 log2FC  = row$estimate, #Estimate not log2FC
                 pval    = row$p.value,
                 stringsAsFactors = FALSE)
    })
    results <- dplyr::bind_rows(Filter(Negate(is.null), results))
    results$padj <- p.adjust(results$pval, method = p_adj_method)
    results$diff_robust <- NA
    results$method  <- "LMM"
    results$dataset <- dataset_name
    results
  }, error = function(e) {
    warning("[LMM] Failed on ", dataset_name, ": ", e$message)
    NULL
  })
}

## 2d. edgeR ───────────────────────────────────────────────────────────────────
run_edger <- function(ps, variable_col, covariates = covariates,
                      dataset_name, p_adj_method = "fdr", alpha = 0.05) {
  message("[edgeR] ", dataset_name)
  tryCatch({
    
    otu <- as(otu_table(ps), "matrix")
    if (!taxa_are_rows(ps)) otu <- t(otu)
    
    sdata <- as.data.frame(unclass(sample_data(ps)))
    rownames(sdata) <- as.character(sdata$DNA_ID)
    
    common <- intersect(colnames(otu), rownames(sdata))
    
    otu <- otu[, common, drop = FALSE]
    sdata <- sdata[common, , drop = FALSE]
    
    otu <- otu[, order(colnames(otu))]
    sdata <- sdata[order(rownames(sdata)), ]
    
    stopifnot(all(colnames(otu) == rownames(sdata)))

    design_formula <- as.formula(
      paste("~", paste(c(variable_col, covariates), collapse = " + "))
    )
    design <- model.matrix(design_formula, data = sdata)

    
    # Build DGEList
    dge <- DGEList(counts = otu)
    dge <- edgeR::calcNormFactors(dge)    
    dge <- edgeR::estimateDisp(dge, design)
    fit <- glmQLFit(dge, design)
    
    # Test the variable of interest (2nd column after intercept)
    coef_idx <- grep(variable_col, colnames(design))[1]
    qlf <- glmQLFTest(fit, coef = coef_idx)
    tbl <- topTags(qlf, n = Inf, adjust.method = p_adj_method)$table
    tbl$feature <- rownames(tbl)
    tbl$diff_robust <- NA
    standardise(tbl, "edgeR", dataset_name,
                "feature", "logFC", "PValue", "FDR","diff_robust")
  }, error = function(e) {
    warning("[edgeR] Failed on ", dataset_name, ": ", e$message)
    NULL
  })
}

## 2f. metagenomeSeq fitFeatureModel ──────────────────────────────────────────
run_mgs <- function(ps, variable_col, covariates = covariates,
                    dataset_name, p_adj_method = "fdr", alpha = 0.05) {
  message("[metagenomeSeq fitFeatureModel] ", dataset_name)
  tryCatch({
    # Convert phyloseq → MRexperiment
    mre <- phyloseq_to_metagenomeSeq(ps)
    mre <- cumNorm(mre, p = cumNormStat(mre))
    
    sdata  <- as.data.frame(unclass(sample_data(ps)))
    pd     <- AnnotatedDataFrame(sdata)
    pData(mre) <- sdata
    
  
    # Not covariates 
    design_formula <- as.formula(
      paste("~", variable_col)
    )
    mod <- model.matrix(design_formula, data = sdata)
    
    fit <- fitFeatureModel(mre, mod)
    tbl <- MRcoefs(fit, number = nrow(assayData(mre)$counts),
                   adjustMethod = p_adj_method)
    tbl$feature <- rownames(tbl)
    tbl$diff_robust <- NA
    # MRcoefs returns: logFC, se, pvalues, adjPvalues
    standardise(tbl, "metagenomeSeq", dataset_name,
                "feature", "logFC", "pvalues", "adjPvalues","diff_robust")
  }, error = function(e) {
    warning("[metagenomeSeq] Failed on ", dataset_name, ": ", e$message)
    NULL
  })
}


# ── 3. Master runner ───────────────────────────────────────────────────────────

#' Run all DA methods on a list of phyloseq objects
#'
#' @param ps_list       Named list of phyloseq objects
#' @param variable_col  Primary variable of interest (string, must be in sample_data)
#' @param covariates    Character vector of covariate column names
#' @param random_effect Optional random effect column for LMM
#' @param alpha         Significance threshold
#' @param p_adj_method  Adjustment method (default "BH")
#' @param outdir        Directory to save result tables
#' @return List with elements: all_results, significant, consensus
run_da_pipeline <- function(ps_list=ps_list,
                            variable_col=variable_col,
                            covariates    = covariates,
                            random_effect = NULL,
                            alpha         = 0.05,
                            p_adj_method  = "fdr",
                            outdir        = "DA_results",
                            filterCondition=filterCondition) {
  
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  # Build formula strings
  fixed_vars     <- paste(c(variable_col, covariates), collapse = " + ")
  formula_str    <- paste0("~ ", fixed_vars)            # for LinDA / ANCOMBC
  fixed_formula  <- fixed_vars                          # for ANCOMBC2
  
  all_results <- list()
  
  for (ds_name in names(ps_list)) {
    
    ps <- ps_list[[ds_name]]

    
    #  Filtering: remove samples that are NA  in BMI
    if(filterCondition == FALSE){
    ps<- subset_samples(ps, !is.na(sample_data(ps)$bmi_calculated))
    
    otu <- as(otu_table(ps), "matrix")
    meta <- as(sample_data(ps), "data.frame")
    ps <- prune_taxa(taxa_sums(ps) > 0, ps)
    otu <- as(otu_table(ps), "matrix")
    meta <- as(sample_data(ps), "data.frame")
    
    #  Filtering: remove sparse/low-abundance features
    
    prev_thresh <- 0.1    # keep features in ≥10% of samples
    keep <- which(rowMeans(otu > 0) >= prev_thresh )
    otu <- otu[keep, , drop = FALSE]
    ps<-prune_taxa(taxa_names(ps) %in% names(keep), ps)
    
    }else{
      #Filter phyloseq object for DiseaseStatus to run relapse analysis
      ps<-prune_samples(sample_data(ps)$DiseaseStatus  =="BreastCancer",ps)
      ps <- subset_samples(ps, 
                           !is.na(relapse) &
                             !is.na(bmi_calculated) &
                             !is.na(ki67_ihc) &
                             !is.na(bmi_calculated) &
                             !is.na(ln_tumorpositive.binary) &
                             !is.na(tumor_size)
      )
      
      ps <- prune_taxa(taxa_sums(ps) > 0, ps)
      otu <- as(otu_table(ps), "matrix")
      meta <- as(sample_data(ps), "data.frame")
      
      #  Filtering: remove sparse/low-abundance features
      
      prev_thresh <- 0.1    # keep features in ≥10% of samples
      keep <- which(rowMeans(otu > 0) >= prev_thresh )
      otu <- otu[keep, , drop = FALSE]
      ps<-prune_taxa(taxa_names(ps) %in% names(keep), ps)
      
    }
    

    
    #------------------------------------------------------------
    message("\n══ Dataset: ", ds_name, " ══")
    
    res_linda   <- run_linda(ps, formula_str, ds_name, p_adj_method ,  alpha = 0.1)
    res_lmm     <- run_lmm(ps, variable_col, covariates, random_effect,
                           ds_name, p_adj_method, alpha)
    
    if(ds_name == "Taxonomy"){
    res_anc2    <- run_ancombc2(ps, fixed_formula, ds_name, p_adj_method, alpha,
                                rand_formula = if (!is.null(random_effect))
                                  paste0("~1|", random_effect) else NULL)
    res_edger   <- run_edger(ps, variable_col, covariates, ds_name, p_adj_method, alpha)
    res_mgs     <- run_mgs(ps, variable_col, covariates, ds_name, p_adj_method, alpha)
    
    combined <- dplyr::bind_rows(
      res_linda, res_anc2, res_lmm, res_mgs )
    }else{combined <- dplyr::bind_rows( res_linda, res_lmm)}
   
    
    all_results[[ds_name]] <- combined
    
    # Save per-dataset raw results
    write.csv(combined,
              file = file.path(outdir, paste0(ds_name, "_all_methods.csv")),
              row.names = FALSE)
  }
  
  # ── 4. Filter significant features ──────────────────────────────────────────
  all_df <- dplyr::bind_rows(all_results)
  
  significant <- all_df %>%
    dplyr::filter(!is.na(padj) & padj < alpha ) %>%
    dplyr::arrange(dataset, method, padj) 
  
  significant<-dplyr::filter(significant, is.na(diff_robust) | diff_robust!=FALSE  )
 
  
  write.csv(significant,
            file = file.path(outdir,  "significant_all_methods.csv"),
            row.names = FALSE)
  
  # ── 5. Consensus features ────────────────────────────────────────────────────
  n_methods <- length(unique(all_df$method))
  
  consensus <- significant %>%
    dplyr::group_by(dataset, feature) %>%
    dplyr::summarise(
      n_methods_sig  = dplyr::n_distinct(method),
      methods_agreed = paste(sort(unique(method)), collapse = ", "),
      mean_log2FC    = mean(log2FC, na.rm = TRUE),
      sd_log2FC      = sd(log2FC,   na.rm = TRUE),
      min_padj       = min(padj,    na.rm = TRUE),
      .groups        = "drop"
    ) %>%
    dplyr::arrange(dataset, dplyr::desc(n_methods_sig), min_padj)
  
  write.csv(consensus,
            file = file.path(outdir, "consensus_features.csv"),
            row.names = FALSE)
  
  message("\n✔ Results saved to: ", outdir)
  message("  • Significant hits  : ", nrow(significant), " rows")

list(
    all_results = all_df,
    significant = significant,
    consensus   = consensus
  )
}

