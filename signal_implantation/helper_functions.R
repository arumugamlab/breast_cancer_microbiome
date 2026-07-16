# ------------------------------------------------------------
# Phyloseq helper
# ------------------------------------------------------------

#' Filter taxa in a phyloseq object by prevalence
#'
#' Keeps only taxa present in at least a specified fraction of samples.
#' Prevalence is defined as the fraction of samples with abundance strictly
#' greater than zero.
#'
#' @param ps A \code{phyloseq} object.
#' @param prevalence_threshold Numeric scalar between 0 and 1 giving the minimum
#'   prevalence required for a taxon to be retained.
#'
#' @return A filtered \code{phyloseq} object containing only taxa with
#'   prevalence greater than or equal to \code{prevalence_threshold}.
#'
#' @examples
#' # ps_filt <- filter_taxa_by_prevalence(ps, prevalence_threshold = 0.1)
filter_taxa_by_prevalence <- function(ps, prevalence_threshold) {
  # Check that the input is a phyloseq object.
  if (!inherits(ps, "phyloseq")) {
    stop("ps must be a phyloseq object.", call. = FALSE)
  }
  
  # Check that the threshold is a single numeric value in [0, 1].
  if (!is.numeric(prevalence_threshold) || length(prevalence_threshold) != 1 ||
      prevalence_threshold < 0 || prevalence_threshold > 1) {
    stop("prevalence_threshold must be a numeric value between 0 and 1.", call. = FALSE)
  }
  
  # Extract the OTU table as a matrix and ensure taxa are rows.
  otu <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) {
    otu <- t(otu)
  }
  
  # Compute prevalence as the fraction of samples with non-zero abundance.
  prevalence <- rowSums(otu > 0, na.rm = TRUE) / ncol(otu)
  
  # Identify taxa that pass the prevalence threshold.
  keep_taxa <- names(prevalence)[prevalence >= prevalence_threshold]
  
  # Subset the phyloseq object to retained taxa only.
  ps_filtered <- prune_taxa(keep_taxa, ps)
  
  ps_filtered
}


#' Add a small in-panel legend to a facet-grid upper-bound plot
#'
#' Returns ggplot layers that can be added with `+` to an existing plot.
#' The legend is drawn only inside the facet identified by `legend_test`
#' and `legend_stratum`.
#'
#' @param legend_test Character scalar. Method/facet row where the legend should appear.
#' @param legend_stratum Character scalar. Abundance stratum/facet column where the legend should appear.
#' @param x Numeric vector of length 2 giving the start and end x positions of legend lines.
#' @param y_top Numeric scalar. y position for the first legend line.
#' @param y_gap Numeric scalar. Vertical gap between the first and second legend lines.
#' @param text_gap Numeric scalar. Horizontal gap between line end and label text.
#' @param threshold_label Character scalar. Label for the blue threshold line.
#' @param observed_label Character scalar. Label for the black observed-envelope line.
#' @param threshold_color Character scalar. Color for the threshold line.
#' @param observed_color Character scalar. Color for the observed-envelope line.
#' @param threshold_linetype Character scalar/integer. Linetype for the threshold line.
#' @param observed_linetype Character scalar/integer. Linetype for the observed-envelope line.
#' @param line_width Numeric scalar. Line width for legend segments.
#' @param text_size Numeric scalar. Text size for legend labels.
#'
#' @return A list of ggplot layers that can be added with `+`.
#'
add_inpanel_threshold_legend <- function(
    legend_test,
    legend_stratum,
    test_levels = NULL,
    x_pos = c(0.12, 0.36),
    y_top = 0.78,
    y_gap = 0.10,
    text_gap = 0.04,
    threshold_label = "80% detection barrier",
    observed_label = "Observed effect q95",
    threshold_color = "#2C6BEA",
    observed_color = "black",
    threshold_linetype = "31",
    observed_linetype = "solid",
    line_width = 0.9,
    text_size = 2.5
) {
  stopifnot(length(x_pos) == 2)
  
  seg_df <- tibble::tibble(
    test = factor(c(legend_test, legend_test), levels = test_levels),
    col_label = "Abundance stratum",
    abund_stratum = factor(c(legend_stratum, legend_stratum),
                           levels = c("Low", "Medium", "High")),
    x = c(x_pos[1], x_pos[1]),
    xend = c(x_pos[2], x_pos[2]),
    y = c(y_top, y_top - y_gap),
    yend = c(y_top, y_top - y_gap),
    line_group = c("threshold", "observed")
  )
  
  txt_df <- tibble::tibble(
    test = factor(c(legend_test, legend_test), levels = test_levels),
    col_label = "Abundance stratum",
    abund_stratum = factor(c(legend_stratum, legend_stratum),
                           levels = c("Low", "Medium", "High")),
    x = c(x_pos[2] + text_gap, x_pos[2] + text_gap),
    y = c(y_top, y_top - y_gap),
    label = c(threshold_label, observed_label)
  )
  
  list(
    ggplot2::geom_segment(
      data = seg_df %>% dplyr::filter(line_group == "threshold"),
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE,
      colour = threshold_color,
      linetype = threshold_linetype,
      linewidth = line_width
    ),
    ggplot2::geom_segment(
      data = seg_df %>% dplyr::filter(line_group == "observed"),
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE,
      colour = observed_color,
      linetype = observed_linetype,
      linewidth = line_width
    ),
    ggplot2::geom_label(
      data = txt_df,
      ggplot2::aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 0.5,
      size = text_size,
      colour = "black",
      fill = "white",
      linewidth = 0.2,
      label.r = grid::unit(0.08, "lines"),
      label.padding = grid::unit(0.12, "lines")
    )
  )
}

#' Add a small in-panel legend to a facet-grid p-value plot
#'
#' Returns ggplot layers that can be added with `+` to an existing plot.
#' The legend is drawn only inside the facet identified by `legend_test`
#' and `legend_stratum`.
#'
#' @param legend_test Character scalar. Method/facet row where the legend should appear.
#' @param legend_stratum Character scalar. Abundance stratum/facet column where the legend should appear.
#' @param x Numeric vector of length 2 giving the start and end x positions of legend lines.
#' @param y_top Numeric scalar. y position for the first legend line.
#' @param y_gap Numeric scalar. Vertical gap between the first and second legend lines.
#' @param text_gap Numeric scalar. Horizontal gap between line end and label text.
#' @param threshold_label Character scalar. Label for the blue threshold line.
#' @param observed_label Character scalar. Label for the black observed-envelope line.
#' @param threshold_color Character scalar. Color for the threshold line.
#' @param observed_color Character scalar. Color for the observed-envelope line.
#' @param threshold_linetype Character scalar/integer. Linetype for the threshold line.
#' @param observed_linetype Character scalar/integer. Linetype for the observed-envelope line.
#' @param line_width Numeric scalar. Line width for legend segments.
#' @param text_size Numeric scalar. Text size for legend labels.
#'
#' @return A list of ggplot layers that can be added with `+`.
#'
add_inpanel_pval_legend <- function(
    legend_test,
    legend_stratum,
    test_levels = NULL,
    x_pos = c(0.12, 0.36),
    y_top = 0.78,
    y_gap = 0.10,
    text_gap = 0.04,
    threshold_label = "80% detection barrier",
    observed_label = "Observed effect q95",
    threshold_color = "#2C6BEA",
    observed_color = "black",
    threshold_linetype = "31",
    observed_linetype = "solid",
    line_width = 0.9,
    text_size = 2.5
) {
  stopifnot(length(x_pos) == 2)
  
  seg_df <- tibble::tibble(
    test = factor(c(legend_test, legend_test), levels = test_levels),
    abund_stratum = factor(c(legend_stratum, legend_stratum),
                           levels = c("Low", "Medium", "High")),
    x = c(x_pos[1], x_pos[1]),
    xend = c(x_pos[2], x_pos[2]),
    y = c(y_top, y_top - y_gap),
    yend = c(y_top, y_top - y_gap),
    line_group = c("threshold", "observed")
  )
  
  txt_df <- tibble::tibble(
    test = factor(c(legend_test, legend_test), levels = test_levels),
    abund_stratum = factor(c(legend_stratum, legend_stratum),
                           levels = c("Low", "Medium", "High")),
    x = c(x_pos[2] + text_gap, x_pos[2] + text_gap),
    y = c(y_top, y_top - y_gap),
    label = c(threshold_label, observed_label)
  )
  
  list(
    ggplot2::geom_segment(
      data = seg_df %>% dplyr::filter(line_group == "threshold"),
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE,
      colour = threshold_color,
      linetype = threshold_linetype,
      linewidth = line_width
    ),
    ggplot2::geom_segment(
      data = seg_df %>% dplyr::filter(line_group == "observed"),
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE,
      colour = observed_color,
      linetype = observed_linetype,
      linewidth = line_width
    ),
    ggplot2::geom_label(
      data = txt_df,
      ggplot2::aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 0.5,
      size = text_size,
      colour = "black",
      fill = "white",
      linewidth = 0.2,
      label.r = grid::unit(0.08, "lines"),
      label.padding = grid::unit(0.12, "lines")
    )
  )
}


# ------------------------------------------------------------
# Detectability model comparison using AIC and 5-fold CV AUC
# ------------------------------------------------------------

#' Prepare a marker-level detectability table for model comparison
#'
#' Filters a calibration table to implanted markers with complete data for
#' detectability modeling. This is intended for comparing simple logistic
#' regression models of detection status using effect size, prevalence, and
#' baseline abundance.
#'
#' @param master_calib A data frame or tibble containing marker-level
#'   calibration output across methods and replicates.
#' @param method Optional character vector of method names to retain. If NULL,
#'   all methods are kept.
#' @param target_col Character scalar naming the implanted-marker indicator
#'   column.
#' @param good_col Character scalar naming the detection-status column.
#' @param effect_col Character scalar naming the absolute effect-size column.
#' @param prevalence_col Character scalar naming the prevalence column.
#' @param abundance_col Character scalar naming the baseline abundance column.
#'
#' @return A tibble filtered to implanted markers with complete values in the
#'   required columns, with standardized modeling columns added.
prepare_detectability_model_data <- function(
    master_calib,
    method = NULL,
    target_col = "is_marker",
    good_col = "detected",
    effect_col = "abs_effect_size",
    prevalence_col = "prevalence_full",
    abundance_col = "log10_abundance_nonzero_median_baseline"
) {
  req <- c("test", target_col, good_col, effect_col, prevalence_col, abundance_col)
  miss <- setdiff(req, colnames(master_calib))
  if (length(miss) > 0) {
    stop("master_calib is missing required columns: ",
         paste(shQuote(miss), collapse = ", "))
  }
  
  out <- master_calib %>%
    dplyr::filter(.data[[target_col]]) %>%
    dplyr::filter(
      !is.na(.data[[good_col]]),
      !is.na(.data[[effect_col]]),
      !is.na(.data[[prevalence_col]]),
      !is.na(.data[[abundance_col]])
    ) %>%
    dplyr::mutate(
      good_bin = as.integer(.data[[good_col]]),
      effect = as.numeric(.data[[effect_col]]),
      prevalence = as.numeric(.data[[prevalence_col]]),
      abundance = as.numeric(.data[[abundance_col]])
    )
  
  if (!is.null(method)) {
    out <- out %>% dplyr::filter(test %in% method)
  }
  
  out
}


#' Return formulas for nested detectability models
#'
#' Defines the candidate logistic regression formulas used for AIC comparison
#' and cross-validated AUC estimation.
#'
#' @return A named list of formulas.
get_detectability_model_formulas <- function(model_type = "linear") {
  if (model_type == "nonlinear") {
    list(
      #M0 = good_bin ~ te(effect, prevalence, k = c(3, 4)) + s(abundance, k = 4, bs = "ts")
      M0 = good_bin ~ s(effect, k = 4, bs = "ts"),
      
      M1 = good_bin ~
        s(effect, k = 4, bs = "ts") +
        s(prevalence, k = 4, bs = "ts"),
      
      M2 = good_bin ~
        s(effect, k = 4, bs = "ts") +
        s(abundance, k = 4, bs = "ts"),
      
      M3 = good_bin ~
        s(effect, k = 4, bs = "ts") +
        s(prevalence, k = 4, bs = "ts") +
        s(abundance, k = 4, bs = "ts"),
      
      M4 = good_bin ~
        s(effect, k = 4, bs = "ts") +
        s(prevalence, k = 4, bs = "ts") +
        s(abundance, k = 4, bs = "ts") +
        ti(effect, prevalence, k = c(4, 4), bs = c("ts", "ts")),
      
      M5 = good_bin ~
        s(effect, k = 4, bs = "ts") +
        s(prevalence, k = 4, bs = "ts") +
        s(abundance, k = 4, bs = "ts") +
        ti(effect, abundance, k = c(4, 4), bs = c("ts", "ts")),
      
      M6 = good_bin ~
        s(effect, k = 4, bs = "ts") +
        s(prevalence, k = 4, bs = "ts") +
        s(abundance, k = 4, bs = "ts") +
        ti(effect, prevalence, k = c(4, 4), bs = c("ts", "ts")) +
        ti(effect, abundance, k = c(4, 4), bs = c("ts", "ts")),
      
      M7 = good_bin ~
        s(effect, k = 4, bs = "ts") +
        s(prevalence, k = 4, bs = "ts") +
        s(abundance, k = 4, bs = "ts") +
        ti(effect, prevalence, k = c(4, 4), bs = c("ts", "ts")) +
        ti(effect, abundance, k = c(4, 4), bs = c("ts", "ts")) +
        ti(prevalence, abundance, k = c(4, 4), bs = c("ts", "ts"))
    )
  } else {
    list(
      M0 = good_bin ~ effect,
      M1 = good_bin ~ effect + prevalence,
      M2 = good_bin ~ effect + abundance,
      M3 = good_bin ~ effect + prevalence + abundance,
      M4 = good_bin ~ effect + prevalence + abundance + effect:prevalence,
      M5 = good_bin ~ effect + prevalence + abundance + effect:abundance,
      M6 = good_bin ~ effect + prevalence + abundance + effect:prevalence + effect:abundance,
      M7 = good_bin ~ effect + prevalence + abundance + effect:prevalence + effect:abundance + prevalence:abundance
    )
  }
}


#' Fit nested detectability models for one method
#'
#' Fits a small set of nested logistic regression models to assess whether
#' prevalence and baseline abundance improve explanation of marker detectability
#' beyond absolute effect size.
#'
#' @param df A tibble containing one method only, as returned by
#'   \code{prepare_detectability_model_data()}.
#' @param test_name Optional method name used only in messages.
#' @param formulas Optional named list of formulas. If NULL, defaults to
#'   \code{get_detectability_model_formulas()}.
#'
#' @return A named list of fitted \code{gam} objects.
fit_detectability_models_one_method <- function(df, test_name = NULL, formulas = NULL, model_type = "linear") {
  library(mgcv)
  if (!all(c("good_bin", "effect", "prevalence", "abundance") %in% colnames(df))) {
    stop("Input data frame does not contain the required modeling columns.")
  }
  
  if (is.null(formulas)) {
    formulas <- get_detectability_model_formulas(model_type=model_type)
  }
  
  # Require both outcome classes to be present.
  if (length(unique(df$good_bin)) < 2) {
    stop("Method ", test_name %||% "", " does not contain both detected and missed markers.")
  }
  
  purrr::map(
    formulas,
    ~ mgcv::gam(
      formula = .x,
      data = df,
      family = stats::binomial(),
      method = "REML"
    )
  )
}


#' Summarize AIC across detectability models for one method
#'
#' Computes AIC for the nested detectability models and reports delta AIC
#' relative to the effect-size-only baseline model.
#'
#' @param model_list Named list of fitted \code{gam} objects returned by
#'   \code{fit_detectability_models_one_method()}.
#' @param method_name Character scalar naming the method.
#'
#' @return A tibble with one row per model and columns for AIC and delta AIC.
summarize_detectability_model_aic_one_method <- function(model_list, method_name) {
  aic_tbl <- tibble::tibble(
    test = method_name,
    model = names(model_list),
    aic = purrr::map_dbl(model_list, stats::AIC)
  )
  
  aic0 <- aic_tbl$aic[aic_tbl$model == "M0"]
  
  aic_tbl %>%
    dplyr::mutate(delta_aic_vs_M0 = aic - aic0) %>%
    dplyr::arrange(aic)
}


#' Create grouped fold assignments for detectability cross-validation
#'
#' Assigns cross-validation folds at the simulation-group level so that all rows
#' from the same group are kept in the same fold. This avoids leakage arising
#' from splitting highly related rows from the same implanted simulation across
#' training and test folds.
#'
#' If a binary outcome column is supplied, groups are ordered by their mean
#' outcome before assignment to promote rough balance of detection prevalence
#' across folds. This is only a heuristic, but it usually gives more stable
#' folds than purely random group assignment.
#'
#' @param data A data frame containing at least a grouping column and, optionally,
#'   an outcome column.
#' @param group_col Character scalar naming the column that identifies the
#'   simulation group. All rows with the same group value will be assigned to the
#'   same fold.
#' @param y_col Optional character scalar naming a binary outcome column coded as
#'   0/1. If provided, groups are ordered by mean outcome before fold assignment.
#' @param k Integer number of folds.
#' @param seed Optional random seed for reproducibility.
#'
#' @return An integer vector of fold assignments of length \code{nrow(data)},
#'   taking values in \code{1:k}.
make_grouped_folds <- function(
    data,
    group_col = "feature",
    y_col = "good_bin",
    k = 5,
    seed = 1
) {
  if (!group_col %in% colnames(data)) {
    stop("Grouping column not found: ", shQuote(group_col))
  }
  
  if (!is.null(y_col) && !y_col %in% colnames(data)) {
    stop("Outcome column not found: ", shQuote(y_col))
  }
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Build one row per group.
  group_df <- data.frame(
    group = unique(data[[group_col]]),
    stringsAsFactors = FALSE
  )
  
  # Optionally summarize the binary outcome by group to help balance folds.
  if (!is.null(y_col)) {
    y_summary <- stats::aggregate(
      data[[y_col]],
      by = list(group = data[[group_col]]),
      FUN = mean
    )
    colnames(y_summary)[2] <- "y_mean"
    
    n_summary <- stats::aggregate(
      data[[y_col]],
      by = list(group = data[[group_col]]),
      FUN = length
    )
    colnames(n_summary)[2] <- "n_rows"
    
    group_df <- merge(group_df, y_summary, by = "group", all.x = TRUE)
    group_df <- merge(group_df, n_summary, by = "group", all.x = TRUE)
  } else {
    n_summary <- stats::aggregate(
      rep(1, nrow(data)),
      by = list(group = data[[group_col]]),
      FUN = length
    )
    colnames(n_summary)[2] <- "n_rows"
    group_df <- merge(group_df, n_summary, by = "group", all.x = TRUE)
    group_df$y_mean <- NA_real_
  }
  
  # Randomize order first to avoid deterministic artifacts when ties occur.
  group_df <- group_df[sample(seq_len(nrow(group_df))), , drop = FALSE]
  
  # If outcome is available, sort groups by mean outcome (descending) so the
  # round-robin assignment spreads easy/hard groups more evenly across folds.
  if (!all(is.na(group_df$y_mean))) {
    ord <- order(group_df$y_mean, decreasing = TRUE, na.last = TRUE)
    group_df <- group_df[ord, , drop = FALSE]
  }
  
  # Assign groups to folds in round-robin fashion.
  group_df$fold <- rep(seq_len(k), length.out = nrow(group_df))
  
  # Map group-level folds back to row-level folds.
  fold_map <- setNames(group_df$fold, group_df$group)
  folds <- unname(fold_map[data[[group_col]]])
  
  as.integer(folds)
}

#' Compute binary AUC from observed labels and predicted probabilities
#'
#' Computes ROC AUC using the rank-based Mann-Whitney formulation.
#'
#' @param y_true Binary numeric vector coded as 0/1.
#' @param y_prob Numeric vector of predicted probabilities.
#'
#' @return Numeric scalar AUC, or \code{NA_real_} if undefined.
compute_binary_auroc <- function(y_true, y_prob) {
  keep <- stats::complete.cases(y_true, y_prob)
  y_true <- as.integer(y_true[keep])
  y_prob <- as.numeric(y_prob[keep])
  
  n1 <- sum(y_true == 1)
  n0 <- sum(y_true == 0)
  
  # AUC is undefined if one class is absent.
  if (n1 == 0 || n0 == 0) {
    return(NA_real_)
  }
  
  r <- rank(y_prob, ties.method = "average")
  auroc <- (sum(r[y_true == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  
  as.numeric(auroc)
}

#' Compute binary AUPRC from observed labels and predicted probabilities
#'
#' Computes area under the precision-recall curve using a threshold sweep over
#' predicted probabilities.
#'
#' @param y_true Binary numeric vector coded as 0/1.
#' @param y_prob Numeric vector of predicted probabilities.
#'
#' @return Numeric scalar AUPRC, or \code{NA_real_} if undefined.
compute_binary_auprc <- function(y_true, y_prob) {
  keep <- stats::complete.cases(y_true, y_prob)
  y_true <- as.integer(y_true[keep])
  y_prob <- as.numeric(y_prob[keep])
  
  n_pos <- sum(y_true == 1)
  if (n_pos == 0) {
    return(NA_real_)
  }
  
  ord <- order(y_prob, decreasing = TRUE)
  y_true <- y_true[ord]
  y_prob <- y_prob[ord]
  
  tp <- cumsum(y_true == 1)
  fp <- cumsum(y_true == 0)
  
  precision <- tp / (tp + fp)
  recall <- tp / n_pos
  
  # Add anchor point at recall = 0 with precision = 1
  recall <- c(0, recall)
  precision <- c(1, precision)
  
  # Trapezoidal area under PR curve
  sum((recall[-1] - recall[-length(recall)]) *
        (precision[-1] + precision[-length(precision)]) / 2)
}

#' Compute k-fold cross-validated AUC for one method and one model set
#'
#' Fits each candidate logistic model on training folds and evaluates AUC on the
#' held-out fold. Fold creation is stratified by detection status.
#'
#' @param df A tibble containing one method only, as returned by
#'   \code{prepare_detectability_model_data()}.
#' @param method_name Character scalar naming the method.
#' @param formulas Optional named list of formulas. If NULL, defaults to
#'   \code{get_detectability_model_formulas()}.
#' @param k Number of folds for cross-validation.
#' @param seed Optional random seed for fold generation.
#'
#' @return A list with components:
#'   \itemize{
#'     \item \code{fold_auc}: tibble with one row per model per fold
#'     \item \code{summary_auc}: tibble with mean/sd AUC by model
#'   }
cv_detectability_auc_one_method <- function(
    df,
    method_name,
    model_type = "linear",
    formulas = NULL,
    k = 5,
    seed = 1
) {
  if (is.null(formulas)) {
    formulas <- get_detectability_model_formulas(model_type=model_type)
  }
  
  if (length(unique(df$good_bin)) < 2) {
    stop("Method ", method_name, " does not contain both detected and missed markers.")
  }
  
  group_col = "feature"
  folds <- make_grouped_folds(
    data = df,
    group_col = group_col,
    y_col = "good_bin",
    k = k,
    seed = seed
  )
  
  # QC: each grouping unit must be assigned to exactly one fold
  group_fold_map <- stats::aggregate(
    folds,
    by = list(group = df[[group_col]]),
    FUN = function(x) length(unique(x))
  )
  
  bad_groups <- group_fold_map$group[group_fold_map$x != 1]
  if (length(bad_groups) > 0) {
    stop(
      "Grouping QC failed: some ", group_col,
      " values were assigned to multiple folds. Examples: ",
      paste(utils::head(bad_groups, 10), collapse = ", ")
    )
  }
  
  fold_auc <- purrr::map_dfr(seq_len(k), function(f) {
    train_df <- df[folds != f, , drop = FALSE]
    test_df  <- df[folds == f, , drop = FALSE]
    
    purrr::imap_dfr(
      formulas,
      function(frm, model_name) {
        fit <- mgcv::gam(
          formula = frm,
          data = train_df,
          family = stats::binomial(),
          method = "REML"
        )
        
        # Predict probabilities on the held-out fold.
        pred <- stats::predict(fit, newdata = test_df, type = "response")
        auroc <- compute_binary_auroc(test_df$good_bin, pred)
        auprc <- compute_binary_auprc(test_df$good_bin, pred)
        
        tibble::tibble(
          test = method_name,
          model = model_name,
          fold = f,
          auroc = auroc,
          auprc = auprc,
          n_train = nrow(train_df),
          n_test = nrow(test_df),
          positives_test = sum(test_df$good_bin == 1),
          negatives_test = sum(test_df$good_bin == 0)
        )
      })
  })
  
  summary_auc <- fold_auc %>%
    dplyr::group_by(test, model) %>%
    dplyr::summarise(
      mean_auroc = mean(auroc, na.rm = TRUE),
      sd_auroc = stats::sd(auroc, na.rm = TRUE),
      mean_auprc = mean(auprc, na.rm = TRUE),
      sd_auprc = stats::sd(auprc, na.rm = TRUE),
      n_folds = sum(!is.na(auroc)),
      .groups = "drop"
    ) %>%
  dplyr::arrange(dplyr::desc(mean_auroc), model)
  
  list(
    fold_auc = fold_auc,
    summary_auc = summary_auc
  )
}


#' Compare detectability models across methods using AIC and 5-fold CV AUC
#'
#' For each method, fits nested logistic regression models of marker
#' detectability, compares them by AIC, and estimates 5-fold cross-validated
#' AUC. This helps assess both which predictors matter and whether the resulting
#' models predict detectability well enough to drive a smooth detectability
#' boundary.
#'
#' @param master_calib A data frame or tibble containing marker-level
#'   calibration output across methods and replicates.
#' @param methods Optional character vector of methods to analyze. If NULL,
#'   all methods present in \code{master_calib} are used.
#' @param target_col Character scalar naming the implanted-marker indicator
#'   column.
#' @param good_col Character scalar naming the detection-status column.
#' @param effect_col Character scalar naming the absolute effect-size column.
#' @param prevalence_col Character scalar naming the prevalence column.
#' @param abundance_col Character scalar naming the baseline abundance column.
#' @param min_rows Minimum number of rows required to fit models for a method.
#' @param k_folds Number of folds for cross-validation.
#' @param seed Optional random seed for fold generation.
#' @param verbose Logical; if TRUE, emit messages when methods are skipped.
#'
#' @return A list with components:
#'   \itemize{
#'     \item \code{data}: prepared marker-level data used for modeling
#'     \item \code{models}: named list of fitted-model lists by method
#'     \item \code{aic_long}: long-format tibble with AIC summaries
#'     \item \code{aic_wide}: wide-format tibble for quick comparison
#'     \item \code{best_model}: best-supported model per method by AIC
#'     \item \code{cv_fold_auc}: fold-level AUC results
#'     \item \code{cv_summary_auc}: summarized AUC by method and model
#'   }
compare_detectability_models_aic_cv <- function(
    master_calib,
    methods = NULL,
    model_type = "linear",
    target_col = "is_marker",
    good_col = "detected",
    effect_col = "abs_effect_size",
    prevalence_col = "prevalence_full",
    abundance_col = "log10_abundance_nonzero_median_baseline",
    min_rows = 50,
    k_folds = 5,
    seed = 1,
    verbose = TRUE
) {
  dat <- prepare_detectability_model_data(
    master_calib = master_calib,
    method = methods,
    target_col = target_col,
    good_col = good_col,
    effect_col = effect_col,
    prevalence_col = prevalence_col,
    abundance_col = abundance_col
  )
  
  formulas <- get_detectability_model_formulas(model_type=model_type)
  method_names <- sort(unique(dat$test))
  
  models <- list()
  aic_long_list <- list()
  cv_fold_auc_list <- list()
  cv_summary_auc_list <- list()
  
  for (m in method_names) {
    dfm <- dat %>% dplyr::filter(test == m)
    
    # Skip methods with too few rows or only one outcome class.
    if (nrow(dfm) < min_rows) {
      if (isTRUE(verbose)) {
        message("Skipping ", m, ": fewer than ", min_rows, " rows.")
      }
      next
    }
    if (length(unique(dfm$good_bin)) < 2) {
      if (isTRUE(verbose)) {
        message("Skipping ", m, ": detected outcome has only one class.")
      }
      next
    }
    
    # Fit full-data models for AIC comparison.
    fit_list <- fit_detectability_models_one_method(
      df = dfm,
      test_name = m,
      formulas = formulas
    )
    
    models[[m]] <- fit_list
    aic_long_list[[m]] <- summarize_detectability_model_aic_one_method(fit_list, m)
    
    # Estimate k-fold CV AUC for the same model set.
    cv_res <- cv_detectability_auc_one_method(
      df = dfm,
      method_name = m,
      formulas = formulas,
      k = k_folds,
      seed = seed
    )
    
    cv_fold_auc_list[[m]] <- cv_res$fold_auc
    cv_summary_auc_list[[m]] <- cv_res$summary_auc
  }
  
  aic_long <- dplyr::bind_rows(aic_long_list)
  
  aic_wide <- aic_long %>%
    dplyr::select(test, model, aic, delta_aic_vs_M0) %>%
    tidyr::pivot_wider(
      names_from = model,
      values_from = c(aic, delta_aic_vs_M0),
      names_sep = "_"
    )
  
  best_model_aic <- aic_long %>%
    dplyr::group_by(test) %>%
    dplyr::slice_min(order_by = aic, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::rename(best_model = model, best_aic = aic)
  
  cv_fold_auc <- dplyr::bind_rows(cv_fold_auc_list)
  
  cv_summary_auc <- dplyr::bind_rows(cv_summary_auc_list) %>%
    dplyr::group_by(test) %>%
    dplyr::mutate(
      delta_auroc_vs_M0 = mean_auroc - mean_auroc[model == "M0"],
      delta_auprc_vs_M0 = mean_auprc - mean_auprc[model == "M0"]
    ) %>%
    dplyr::ungroup()
  
  best_model_auroc <- cv_summary_auc %>%
    dplyr::group_by(test) %>%
    dplyr::mutate(mean_auroc = round(mean_auroc, 2)) %>%
    dplyr::slice_min(order_by = tibble(-mean_auroc, model), n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::rename(best_model = model)

  best_model_auprc <- cv_summary_auc %>%
    dplyr::group_by(test) %>%
    dplyr::mutate(mean_auprc = round(mean_auprc, 2)) %>%
    dplyr::slice_min(order_by = tibble(-mean_auprc, model), n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::rename(best_model = model)
  
  list(
    data = dat,
    models = models,
    aic_long = aic_long,
    aic_wide = aic_wide,
    best_model_aic = best_model_aic,
    cv_fold_auc = cv_fold_auc,
    cv_summary_auc = cv_summary_auc,
    best_model_auroc = best_model_auroc,
    best_model_auprc = best_model_auprc
  )
}


plot_pseudomiami_detectability_abundance <- function(
    res_bound,
    methods = NULL,
    abundance_levels = c("Low", "Medium", "High"),
    pval_ymax = 4,                 # default = -log10(1e-4)
    top_fraction = 0.5,
    gap_fraction = 0.02,           # continuous; converted to touching margins
    effect_expand = 0.05,
    point_size_failed = 1.3,
    point_alpha_failed = 0.65,
    point_size_positive = 2.4,
    boundary_linewidth = 1.0,
    divider_linewidth = 0.7,
    failed_label = "Species failing the association test",
    positive_label = "Positive taxon",
    boundary_label = "80% detectability boundary (stratum)",
    failed_colour = "#9E9E9E",
    positive_colour = "#D55E00",
    boundary_colour = "#9C1C6B",
    method_strip_width = 0.045,
    method_strip_fill = "grey95",
    method_strip_linewidth = 0.8,
    show_super_header = TRUE,
    super_header_label = "Abundance stratum",
    figure_width = 8,
    figure_height = 10,
    match_header_to_method_strip = TRUE,
    annotate_detectability = TRUE,
    annotation_method = NULL,
    annotation_stratum = "High",
    annotation_x_upper = 0.90,
    annotation_x_lower = 0.90,
    annotation_text_upper = "signals less likely\nto be detected",
    annotation_text_lower = "signals more likely\nto be detected",
    annotation_size = 2.7,
    annotation_colour = "grey20",
    annotation_angle = 0,
    outer_margin_pt = 5.5
) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required.", call. = FALSE)
  }
  if (!requireNamespace("grid", quietly = TRUE)) {
    stop("Package 'grid' is required.", call. = FALSE)
  }
  if (!requireNamespace("gtable", quietly = TRUE)) {
    stop("Package 'gtable' is required.", call. = FALSE)
  }
  
  stopifnot(is.list(res_bound))
  stopifnot(all(c("real_df", "thresh_df") %in% names(res_bound)))
  stopifnot(is.numeric(top_fraction), length(top_fraction) == 1,
            top_fraction > 0, top_fraction < 1)
  stopifnot(is.numeric(gap_fraction), length(gap_fraction) == 1,
            gap_fraction >= 0, gap_fraction < 0.5)
  
  real_df <- res_bound$real_df
  thresh_df <- res_bound$thresh_df
  
  if (is.null(methods)) {
    methods <- unique(real_df$test)
  }
  
  real_df <- real_df %>%
    dplyr::filter(test %in% methods)
  
  thresh_df <- thresh_df %>%
    dplyr::filter(test %in% methods)
  
  if (!"prevalence" %in% names(real_df)) {
    if ("prevalence_full" %in% names(real_df)) {
      real_df$prevalence <- real_df$prevalence_full
    } else {
      stop("real_df must contain 'prevalence' or 'prevalence_full'.", call. = FALSE)
    }
  }
  
  if (!"prev_mid" %in% names(thresh_df)) {
    if ("prevalence" %in% names(thresh_df)) {
      thresh_df$prev_mid <- thresh_df$prevalence
    } else {
      stop("thresh_df must contain 'prev_mid' or 'prevalence'.", call. = FALSE)
    }
  }
  
  if (!"p_raw" %in% names(real_df)) {
    if ("pval" %in% names(real_df)) {
      real_df$p_raw <- real_df$pval
    } else {
      stop("real_df must contain 'p_raw' or 'pval'.", call. = FALSE)
    }
  }
  
  if (!"detected" %in% names(real_df)) {
    real_df$detected <- FALSE
  }
  
  real_df <- real_df %>%
    dplyr::mutate(
      neglog10_p = -log10(pmax(p_raw, .Machine$double.xmin)),
      neglog10_p = pmin(neglog10_p, pval_ymax),
      test = factor(test, levels = methods),
      abund_stratum = factor(abund_stratum, levels = abundance_levels),
      legend_key = ifelse(detected, positive_label, failed_label)
    )
  
  thresh_df <- thresh_df %>%
    dplyr::mutate(
      test = factor(test, levels = methods),
      abund_stratum = factor(abund_stratum, levels = abundance_levels),
      legend_key = boundary_label
    )
  
  has_positive <- any(real_df$detected, na.rm = TRUE)
  
  panel_real_max <- real_df %>%
    dplyr::group_by(test, abund_stratum) %>%
    dplyr::summarise(real_max = max(abs_effect_size, na.rm = TRUE), .groups = "drop")
  
  panel_thresh_max <- thresh_df %>%
    dplyr::group_by(test, abund_stratum) %>%
    dplyr::summarise(thresh_max = max(threshold, na.rm = TRUE), .groups = "drop")
  
  panel_max <- dplyr::full_join(
    panel_real_max, panel_thresh_max,
    by = c("test", "abund_stratum")
  ) %>%
    dplyr::mutate(
      real_max = dplyr::coalesce(real_max, 0),
      thresh_max = dplyr::coalesce(thresh_max, 0),
      panel_effect_max = pmax(real_max, thresh_max),
      panel_effect_max = dplyr::if_else(
        is.finite(panel_effect_max) & panel_effect_max > 0,
        panel_effect_max * (1 + effect_expand),
        1
      )
    ) %>%
    dplyr::select(test, abund_stratum, panel_effect_max)
  
  real_df <- real_df %>%
    dplyr::left_join(panel_max, by = c("test", "abund_stratum"))
  
  thresh_df <- thresh_df %>%
    dplyr::left_join(panel_max, by = c("test", "abund_stratum"))
  
  real_fail <- real_df %>% dplyr::filter(!detected)
  real_pos  <- real_df %>% dplyr::filter(detected)
  
  if (is.null(annotation_method)) {
    annotation_method <- tail(methods, 1)
  }
  
  colour_values <- c(
    setNames(failed_colour, failed_label),
    setNames(boundary_colour, boundary_label)
  )
  if (has_positive) {
    colour_values <- c(colour_values, setNames(positive_colour, positive_label))
  }
  
  colour_breaks <- c(boundary_label, failed_label)
  if (has_positive) {
    colour_breaks <- c(colour_breaks, positive_label)
  }
  
  # Continuous gap control via touching margins only.
  # 0 => no gap. Small positive => tiny separation.
  gap_pt <- 60 * gap_fraction
  
  common_colour_scale <- ggplot2::scale_colour_manual(
    values = colour_values,
    breaks = colour_breaks,
    drop = FALSE,
    name = NULL
  )
  
  common_theme <- ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.5),
      strip.background = ggplot2::element_blank(),
      legend.position = "none"
    )
  
  make_top_panel <- function(method_i, stratum_i, show_y = FALSE) {
    df_fail <- real_fail %>%
      dplyr::filter(test == method_i, abund_stratum == stratum_i)
    
    df_pos <- real_pos %>%
      dplyr::filter(test == method_i, abund_stratum == stratum_i)
    
    p <- ggplot2::ggplot() +
      ggplot2::geom_point(
        data = df_fail,
        ggplot2::aes(x = prevalence, y = neglog10_p, colour = legend_key),
        size = point_size_failed,
        alpha = point_alpha_failed
      )
    
    if (nrow(df_pos) > 0) {
      p <- p +
        ggplot2::geom_point(
          data = df_pos,
          ggplot2::aes(x = prevalence, y = neglog10_p, colour = legend_key),
          size = point_size_positive
        )
    }
    
    p +
      common_colour_scale +
      ggplot2::scale_x_continuous(
        limits = c(0.05, 1.0),
        expand = ggplot2::expansion(mult = c(0.02, 0.04))
      ) +
      ggplot2::scale_y_continuous(
        limits = c(0, pval_ymax),
        breaks = seq(0, pval_ymax, by = 1),
        expand = ggplot2::expansion(mult = c(0, 0.02))
      ) +
      ggplot2::labs(
        x = NULL,
        y = if (show_y) expression(-log[10](p)) else NULL
      ) +
      common_theme +
      ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        plot.margin = grid::unit(c(outer_margin_pt, outer_margin_pt, gap_pt / 2, outer_margin_pt), "pt")
      ) +
      {
        if (!show_y) {
          ggplot2::theme(
            axis.title.y = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank()
          )
        } else {
          ggplot2::theme(axis.text.y = ggplot2::element_text(vjust = 0))
        }
      }
  }
  
  make_bottom_panel <- function(method_i, stratum_i, show_y = FALSE, show_x = FALSE, annotate = FALSE) {
    df_fail <- real_fail %>%
      dplyr::filter(test == method_i, abund_stratum == stratum_i)
    
    df_pos <- real_pos %>%
      dplyr::filter(test == method_i, abund_stratum == stratum_i)
    
    df_thr <- thresh_df %>%
      dplyr::filter(test == method_i, abund_stratum == stratum_i) %>%
      dplyr::filter(
        is.finite(prev_mid),
        is.finite(threshold),
        is.finite(panel_effect_max)
      ) %>%
      dplyr::arrange(prev_mid)
    
    ymax_i <- unique(df_thr$panel_effect_max)
    if (length(ymax_i) == 0 || !is.finite(ymax_i[1])) ymax_i <- 1 else ymax_i <- ymax_i[1]
    
    df_thr <- df_thr %>%
      dplyr::filter(
        prev_mid >= 0.05,
        prev_mid <= 1.0,
        threshold >= 0,
        threshold <= ymax_i
      )
    
    p <- ggplot2::ggplot() +
      ggplot2::geom_hline(yintercept = 0, linewidth = divider_linewidth, colour = "grey35") +
      ggplot2::geom_point(
        data = df_fail,
        ggplot2::aes(x = prevalence, y = abs_effect_size, colour = legend_key),
        size = point_size_failed,
        alpha = point_alpha_failed
      ) +
      ggplot2::geom_line(
        data = df_thr,
        ggplot2::aes(x = prev_mid, y = threshold, colour = legend_key),
        linewidth = boundary_linewidth,
        lineend = "round"
      )
    
    if (nrow(df_pos) > 0) {
      p <- p +
        ggplot2::geom_point(
          data = df_pos,
          ggplot2::aes(x = prevalence, y = abs_effect_size, colour = legend_key),
          size = point_size_positive
        )
    }
    
    if (annotate && nrow(df_thr) >= 2) {
      y_lower <- stats::approx(
        x = df_thr$prev_mid,
        y = df_thr$threshold,
        xout = annotation_x_lower,
        rule = 2
      )$y
      y_upper <- stats::approx(
        x = df_thr$prev_mid,
        y = df_thr$threshold,
        xout = annotation_x_upper,
        rule = 2
      )$y
      
      p <- p +
        ggplot2::annotate(
          "text",
          x = annotation_x_upper,
          y = y_upper * 0.70,
          label = annotation_text_upper,
          hjust = 1,
          vjust = 0.5,
          size = annotation_size,
          colour = annotation_colour,
          lineheight = 0.95,
          angle = annotation_angle
        ) +
        ggplot2::annotate(
          "text",
          x = annotation_x_lower,
          y = min(y_lower * 1.4, ymax_i * 0.92),
          label = annotation_text_lower,
          hjust = 1,
          vjust = 1,
          size = annotation_size,
          colour = annotation_colour,
          lineheight = 0.95,
          angle = annotation_angle
        )
    }
    
    p +
      common_colour_scale +
      ggplot2::scale_x_continuous(
        limits = c(0.05, 1.0),
        expand = ggplot2::expansion(mult = c(0.02, 0.04))
      ) +
      ggplot2::scale_y_reverse(
        limits = c(ymax_i, 0),
        expand = ggplot2::expansion(mult = c(0.02, 0))
      ) +
      ggplot2::labs(
        x = if (show_x) "Species prevalence" else NULL,
        y = if (show_y) "abs(eff. size)" else NULL
      ) +
      common_theme +
      ggplot2::theme(
        plot.margin = grid::unit(c(gap_pt / 2, outer_margin_pt, outer_margin_pt, outer_margin_pt), "pt")
      ) +
      {
        if (!show_y) {
          ggplot2::theme(
            axis.title.y = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank()
          )
        } else {
          ggplot2::theme(
            axis.text.y = ggplot2::element_text(vjust = 1)
          )
        }
      } +
      {
        if (!show_x) ggplot2::theme(
          axis.title.x = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank()
        )
      }
  }
  
  make_method_strip <- function(label) {
    patchwork::wrap_elements(
      full = grid::grobTree(
        grid::rectGrob(
          width = grid::unit(1, "npc"),
          height = grid::unit(1, "npc"),
          gp = grid::gpar(
            fill = method_strip_fill,
            col = "black",
            lwd = method_strip_linewidth
          )
        ),
        grid::textGrob(
          label = label,
          rot = 270,
          gp = grid::gpar(fontsize = 11)
        )
      )
    )
  }
  
  make_header_cell <- function(label) {
    patchwork::wrap_elements(
      full = grid::grobTree(
        grid::rectGrob(
          gp = grid::gpar(fill = "grey95", col = "black", lwd = 0.8)
        ),
        grid::textGrob(
          label = label,
          gp = grid::gpar(fontsize = 11)
        )
      )
    )
  }
  
  make_blank_header <- function() {
    patchwork::wrap_elements(
      full = grid::nullGrob()
    )
  }
  
  create_xtick_plot <- function(
    abundance_levels,
    left_gutter_width = -0.05,
    method_strip_width = 0.045,
    x_breaks = c(0.25, 0.50, 0.75, 1.00),
    x_limits = c(0, 1),
    tick_y0 = 0.6,
    tick_y1 = 0.95,
    axis_y = 0.95,
    text_y = 0.1,
    line_size = 0.3,
    text_size = 9,
    label_fun = function(x) formatC(x, format = "f", digits = 2)
  ) {
    stopifnot(length(x_limits) == 2, diff(x_limits) > 0)
    
    # Convert breaks to within-panel NPC coordinates
    tick_pos <- (x_breaks - x_limits[1]) / diff(x_limits)
    keep <- is.finite(tick_pos) & tick_pos >= 0 & tick_pos <= 1
    tick_pos <- tick_pos[keep]
    x_breaks <- x_breaks[keep]
    
    make_tick_cell <- function(
    padding_left = 0.03,
    padding_right = 0.03,
    y_top = 1,
    tick_length = 0.2
    ) {
      tick_df <- tibble::tibble(
        x = x_breaks,
        y = 0
      )
      
      ggplot2::ggplot(tick_df, ggplot2::aes(x = x, y = y)) +
        # horizontal axis line only from x = 0 to x = 1
        ggplot2::geom_segment(
          inherit.aes = FALSE,
          aes(x = 0, xend = 1, y = y_top, yend = y_top),
          linewidth = line_size,
          colour = "black"
        ) +
        # ticks
        ggplot2::geom_segment(
          data = tick_df,
          ggplot2::aes(
            x = x, xend = x,
            y = y_top-tick_length, yend = y_top
          ),
          inherit.aes = FALSE,
          linewidth = line_size,
          colour = "black"
        ) +
        # labels
        ggplot2::geom_text(
          data = tick_df,
          ggplot2::aes(x = x, y = 0, label = x),
          inherit.aes = FALSE,
          size = text_size / ggplot2::.pt,
          vjust = 0
        ) +
        ggplot2::scale_x_continuous(
          limits = c(0 - padding_left, 1 + padding_right),
          expand = c(0, 0)
        ) +
        ggplot2::scale_y_continuous(
          limits = c(0, 1),
          expand = c(0, 0)
        ) +
        ggplot2::theme_void() +
        ggplot2::theme(
          plot.margin = ggplot2::margin(0, 0, 0, 0)
        )
    }
    
    tick_cells <- lapply(seq_along(abundance_levels), function(i) make_tick_cell(padding_left = 0.03, padding_right = 0.03))
    
    xtick_row <- patchwork::wrap_plots(
      list(
        make_blank_header(),   # left gutter
        tick_cells[[1]],
        tick_cells[[2]],
        tick_cells[[3]],
        make_blank_header()    # right strip spacer
      ),
      nrow = 1
    ) + patchwork::plot_layout(
      widths = c(left_gutter_width, 1, 1, 1, method_strip_width)
    )
    
    xtick_row
  }
  
  # Dedicated legend plot to avoid collected-guide inconsistencies
  legend_df <- tibble::tibble(
    x = c(1, 2),
    y = c(1, 1),
    legend_key = c(boundary_label, failed_label)
  )
  if (has_positive) {
    legend_df <- dplyr::bind_rows(
      legend_df,
      tibble::tibble(x = 3, y = 1, legend_key = positive_label)
    )
  }
  
  legend_line_df <- tibble::tibble(
    x = c(0, 1),
    y = c(1, 1),
    legend_key = boundary_label
  )
  
  legend_point_df <- tibble::tibble(
    x = seq_len(if (has_positive) 2 else 1),
    y = 1,
    legend_key = c(failed_label, if (has_positive) positive_label)
  )
  
  legend_plot <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = legend_line_df,
      ggplot2::aes(x = x, y = y, colour = legend_key),
      linewidth = boundary_linewidth
    ) +
    ggplot2::geom_point(
      data = legend_point_df,
      ggplot2::aes(x = x, y = y, colour = legend_key),
      size = 3.2
    ) +
    ggplot2::geom_point(
      data = legend_df %>% dplyr::filter(legend_key != boundary_label),
      ggplot2::aes(x = x, y = y, colour = legend_key),
      size = c(rep(3.0, nrow(legend_df) - 1))
    ) +
    ggplot2::scale_colour_manual(
      values = colour_values,
      breaks = colour_breaks,
      drop = FALSE,
      name = NULL
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        nrow = 2,
        byrow = TRUE,
        override.aes = list(
          linetype = c(1, rep(0, length(colour_breaks) - 1)),
          shape = c(NA, rep(16, length(colour_breaks) - 1)),
          size = c(boundary_linewidth, rep(3.2, length(colour_breaks) - 1)),
          alpha = rep(1, length(colour_breaks))
        )
      )
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "top",
      legend.box = "vertical",
      legend.key.width = grid::unit(0.9, "lines"),
      legend.key.height = grid::unit(0.7, "lines"),
      legend.spacing.y = grid::unit(0, "pt"),
      legend.spacing.x = grid::unit(2, "pt"),
      legend.margin = ggplot2::margin(t = 6, r = 4, b = 6, l = 4),
      legend.box.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0)
    )
  
  # Header
  
  left_gutter_width <- -0.05

  if (isTRUE(show_super_header)) {
    super_header_main <- patchwork::wrap_elements(
      full = grid::grobTree(
        grid::rectGrob(
          gp = grid::gpar(fill = "grey95", col = "black", lwd = 0.8)
        ),
        grid::textGrob(
          label = super_header_label,
          gp = grid::gpar(fontsize = 11)
        )
      )
    )
    
    super_header <- patchwork::wrap_plots(
      list(
        make_blank_header(),
        super_header_main,
        make_blank_header()
      ),
      nrow = 1
    ) + patchwork::plot_layout(
      widths = c(left_gutter_width, 3, method_strip_width)
    )
  } else {
    super_header <- NULL
  }
  
  header_cells <- lapply(abundance_levels, make_header_cell)
  header_row <- patchwork::wrap_plots(
    list(
      make_blank_header(),   # left gutter
      header_cells[[1]],
      header_cells[[2]],
      header_cells[[3]],
      make_blank_header()    # right strip spacer
    ),
    nrow = 1
  ) + patchwork::plot_layout(
    widths = c(left_gutter_width, 1, 1, 1, method_strip_width)
  )
  
  row_plots <- vector("list", length(methods))
  
  for (i in seq_along(methods)) {
    method_i <- methods[i]
    tile_list <- vector("list", length(abundance_levels))
    
    for (j in seq_along(abundance_levels)) {
      stratum_j <- abundance_levels[j]
      
      p_top <- make_top_panel(
        method_i = method_i,
        stratum_i = stratum_j,
        show_y = (j == 1)
      )
      
      p_bottom <- make_bottom_panel(
        method_i = method_i,
        stratum_i = stratum_j,
        show_y = (j == 1),
        show_x = FALSE,
        annotate = isTRUE(annotate_detectability) &&
          identical(as.character(method_i), as.character(annotation_method)) &&
          identical(as.character(stratum_j), as.character(annotation_stratum))
      )
      
      tile_list[[j]] <- (
        p_top / p_bottom
      ) + patchwork::plot_layout(
        heights = c(top_fraction, 1 - top_fraction)
      )
    }
    
    row_patch <- tile_list[[1]]
    for (j in 2:length(tile_list)) {
      row_patch <- row_patch | tile_list[[j]]
    }
    row_plots[[i]] <- (
      row_patch | make_method_strip(as.character(method_i))
    ) + patchwork::plot_layout(
      widths = c(rep(1, length(tile_list)), method_strip_width)
    )
  }
  
  body <- patchwork::wrap_plots(row_plots, ncol = 1)
  
  # abund_stratum headers
  if (isTRUE(match_header_to_method_strip)) {
    n_cols <- length(abundance_levels)
    r <- figure_width / figure_height
    a <- (method_strip_width / (n_cols + method_strip_width)) * r
    
    if (isTRUE(show_super_header)) {
      sub_header_height <- a / (1 - 2 * a)
      super_header_height <- sub_header_height
    } else {
      sub_header_height <- a / (1 - a)
      super_header_height <- NULL
    }
  }
  
  xtick_plot <- create_xtick_plot(
    abundance_levels = abundance_levels,
    left_gutter_width = left_gutter_width,
    method_strip_width = method_strip_width,
    x_breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.00),
    x_limits = c(0, 1),
    text_size = 9
  )

  assembled <- if (isTRUE(show_super_header)) {
    (
      super_header / header_row / body / xtick_plot
    ) + patchwork::plot_layout(heights = c(super_header_height, sub_header_height, 0.98, 0.02))
  } else {
    (
      header_row / body / xtick_plot
    ) + patchwork::plot_layout(heights = c(sub_header_height, 0.95, 0.05))
  }
  
  xlab_plot <- cowplot::ggdraw() +
    cowplot::draw_label("Species prevalence", size = 12)
  
  # legend box
  n_cols <- length(abundance_levels)
  legend_width <- 1 / n_cols
  legend_grob <- gtable::gtable_filter(ggplot2::ggplotGrob(legend_plot), "guide-box")
  legend_grob <- grid::grobTree(
    grid::rectGrob(
      gp = grid::gpar(fill = "white", col = "black", lwd = 0.8)
    ),
    legend_grob
  )
  legend_plot <- cowplot::ggdraw() +
    cowplot::draw_grob(
      legend_grob,
      x = 1 - legend_width - 0.02,
      y = 0.28,
      width = legend_width,
      height = 0.58
      )
  
  final_plot <- cowplot::plot_grid(
    assembled,
    xlab_plot,
    legend_plot,
    ncol = 1,
    rel_heights = c(1, 0.02, 0.08)
  )
  
  return(final_plot)
  
}

#' Plot a compact recurrence volcano panel for one method
#'
#' @param real_df Real-data association table, e.g. res_bound_abund_recurrence$real_df.
#' @param method Method to plot, e.g. "LinDA" or "lm".
#' @param candidate_feature Feature to highlight.
#' @param method_col Method column.
#' @param feature_col Feature identifier column.
#' @param effect_col Signed effect-size column.
#' @param pval_col P-value column.
#' @param detected_col Logical detected/significant column.
#' @param p_floor Lower p-value floor before -log10 transformation.
#' @param show_p_threshold Whether to draw a nominal p-value threshold line.
#' @param p_threshold Nominal p-value threshold for horizontal line.
#' @param title Optional plot title. Defaults to method name.
#'
#' @return ggplot object.
#' @export
plot_recurrence_volcano <- function(
    real_df,
    method,
    candidate_feature = "meta_mOTU_v3_12389",
    method_col = "test",
    feature_col = "feature",
    effect_col = "effect_size",
    pval_col = "pval",
    detected_col = "detected",
    p_floor = 1e-300,
    show_p_threshold = FALSE,
    p_threshold = 0.05,
    title = NULL
) {
  required_cols <- c(
    method_col, feature_col, effect_col,
    pval_col, detected_col
  )
  
  missing_cols <- setdiff(required_cols, colnames(real_df))
  if (length(missing_cols) > 0) {
    stop(
      "real_df is missing columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  
  df <- real_df %>%
    dplyr::filter(.data[[method_col]] == method) %>%
    dplyr::mutate(
      .effect = .data[[effect_col]],
      .pval = pmax(.data[[pval_col]], p_floor),
      .neglog10p = -log10(.pval),
      .detected = as.logical(.data[[detected_col]]),
      .is_candidate = .data[[feature_col]] == candidate_feature
    ) %>%
    dplyr::filter(
      is.finite(.effect),
      is.finite(.neglog10p)
    )
  
  if (nrow(df) == 0) {
    stop("No rows available for method: ", method, call. = FALSE)
  }
  
  candidate_df <- df %>%
    dplyr::filter(.is_candidate)
  
  if (nrow(candidate_df) == 0) {
    warning("candidate_feature was not found for method: ", method)
  }
  
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .effect, y = .neglog10p)
  ) +
    ggplot2::geom_point(
      data = df %>% dplyr::filter(!.is_candidate),
      colour = "grey65",
      alpha = 0.75,
      size = 1.1
    ) +
    ggplot2::geom_point(
      data = df %>% dplyr::filter(.detected & !.is_candidate),
      colour = "grey30",
      alpha = 0.9,
      size = 1.4
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      colour = "grey45",
      linewidth = 0.35,
      linetype = "dashed"
    )
  
  if (show_p_threshold) {
    p <- p +
      ggplot2::geom_hline(
        yintercept = -log10(p_threshold),
        colour = "grey45",
        linewidth = 0.35,
        linetype = "dotted"
      )
  }
  
  p <- p +
    ggplot2::geom_point(
      data = candidate_df,
      ggplot2::aes(x = .effect, y = .neglog10p),
      inherit.aes = FALSE,
      colour = "#D55E00",
      size = 2.8
    ) +
    ggplot2::labs(
      title = title,
      x = "Effect size",
      y = expression(-log[10](p))
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      aspect.ratio = 1,
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 13,
        hjust = 0.5
      ),
      axis.title = ggplot2::element_text(size = 10),
      axis.text = ggplot2::element_text(size = 9.5),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey88"),
      plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 5.5)
    )
  
  p
}

