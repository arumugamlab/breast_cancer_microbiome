# A. Differential Abundance Pipeline

**Purpose:** Multi-method differential abundance (DA) analysis of microbiome count/TPM data stored as `phyloseq` objects, with automatic consensus calling across methods.

---

## Overview

This pipeline runs up to five complementary DA methods on one or more `phyloseq` datasets, harmonises their outputs into a common schema, filters significant hits, and identifies **consensus features** — features flagged by more than one method. This reduces false-discovery risk inherent to any single method.

```
ps_list (named list of phyloseq objects)
        │
        ▼
  run_da_pipeline()
        │
        ├─► LinDA          (all dataset types)
        ├─► LM / LMM       (all dataset types)
        ├─► ANCOMBC2       (Taxonomy only)
        ├─► edgeR          (Taxonomy only)
        └─► metagenomeSeq  (Taxonomy only)
                │
                ▼
        Harmonised results table
                │
                ├─► <dataset>_all_methods.csv   (per dataset, all features)
                ├─► significant_all_methods.csv (padj < alpha)
                └─► consensus_features.csv      (cross-method agreement)
```

---

## Dependencies

Install all packages before running. Bioconductor packages require `BiocManager`.

```r
# CRAN
install.packages(c("tidyverse", "ggplot2", "ggrepel", "patchwork",
                   "lme4", "lmerTest", "broom.mixed"))

# Bioconductor
BiocManager::install(c("phyloseq", "MicrobiomeStat", "ANCOMBC",
                       "edgeR", "metagenomeSeq", "LinDA"))
```

| Package | Version tested | Role |
|---|---|---|
| phyloseq | ≥ 1.44 | Data container |
| LinDA | ≥ 1.0 | LinDA DA test |
| ANCOMBC | ≥ 2.2 | ANCOMBC2 DA test |
| edgeR | ≥ 3.40 | Negative-binomial GLM |
| metagenomeSeq | ≥ 1.40 | Zero-inflated log-normal |
| lme4 / lmerTest | ≥ 1.1 | CLR-based linear (mixed) models |
| broom.mixed | ≥ 0.2 | Tidy model output |
| tidyverse | ≥ 2.0 | Data wrangling |

---

## Input

### `ps_list` — Named list of phyloseq objects

Each element is a `phyloseq` object with:

| Component | Required | Notes |
|---|---|---|
| `otu_table` | ✅ | Raw counts (Taxonomy) / TPM (Functional profiles); features × samples orientation |
| `sample_data` | ✅ | Must contain `variable_col`,and all `covariates` |
| `tax_table` | Optional | Used for annotation only |

**Example:**
```r
ps_list <- list(
  Taxonomy      = ps_mOTUs,
  KEGG_Pathway  = FA.KEGG_Pathway,
  KEGG_KO       = FA.KEGG_KO 
)
```

> **Dataset type matters:** datasets  named `"Taxonomy"` receive the full five-method suite (LinDA + ANCOMBC2 + LMM + edgeR + metagenomeSeq). All other names (Functional profiles, TPM normalized) receive LinDA + LMM only.

---

### `run_da_pipeline()` — Master function arguments

| Argument | Type | Default | Description |
|---|---|---|---|
| `ps_list` | named list | — | List of `phyloseq` objects (see above) |
| `variable_col` | character | — | Primary variable of interest; must be a column in `sample_data` (e.g. `"DiseaseStatus"`, `"relapse"`) |
| `covariates` | character vector | — | Adjustment covariates (e.g. `c("Age", "bmi_calculated")`). Included in all model formulae |
| `random_effect` | character or NULL | `NULL` | Column name for random intercept in LMM (e.g. `"PatientID"` for repeated-measures). Pass `NULL` for standard LM |
| `alpha` | numeric | `0.05` | Adjusted p-value significance threshold |
| `p_adj_method` | character | `"fdr"` | Multiple-testing correction method passed to all methods (any value accepted by `p.adjust`) |
| `outdir` | character | `"DA_results"` | Directory where all output CSVs are written. Created automatically if absent |
| `filterCondition` | logical | — | **`FALSE`**: runs on all samples, removing those with missing BMI. **`TRUE`**: subsets to `BreastCancer` samples only and additionally requires non-missing `relapse`, `ki67_ihc`, `ln_tumorpositive.binary`, and `tumor_size` |

**Minimal example:**
```r
results <- run_da_pipeline(
  ps_list        = ps_list,
  variable_col   = "DiseaseStatus",
  covariates     = c("Age", "bmi_calculated"),
  random_effect  = NULL,
  alpha          = 0.05,
  p_adj_method   = "fdr",
  outdir         = "DA_results",
  filterCondition = FALSE
)
```

---

### Required columns in `sample_data`

| Column | Used when | Description |
|---|---|---|
| `DNA_ID` | Always | Sample identifier; used to align OTU table and metadata |
| `bmi_calculated` | Always | Samples with `NA` are removed before analysis |
| `DiseaseStatus` | `filterCondition = TRUE` | Samples are subset to `"BreastCancer"` |
| `relapse` | `filterCondition = TRUE` | Must be non-missing |
| `ki67_ihc` | `filterCondition = TRUE` | Must be non-missing |
| `ln_tumorpositive.binary` | `filterCondition = TRUE` | Must be non-missing |
| `tumor_size` | `filterCondition = TRUE` | Must be non-missing |

---

## Methods

### Feature filtering (applied before all methods)

- Samples with `NA` in `bmi_calculated` (and additional columns when `filterCondition = TRUE`) are removed.
- Taxa with zero total counts are pruned.
- Taxa present in **< 10 % of samples** are removed (prevalence threshold).

### Transformations per method

| Method | Transformation | Notes |
|---|---|---|
| LinDA | Relative abundance (`Taxonomy` dataset only) | Proportion normalisation before LinDA |
| ANCOMBC2 | Internal (bias-correction) | Handles compositional bias natively |
| LM / LMM | CLR (centred log-ratio) after pseudo-count | Standard compositional approach |
| edgeR | TMM normalisation | Negative-binomial dispersion estimated per dataset |
| metagenomeSeq | Cumulative-sum scaling (CSS) | Quantile-based normalisation |

---

## Output

### Return value

`run_da_pipeline()` returns a named list with three elements:

| Element | Class | Description |
|---|---|---|
| `all_results` | `data.frame` | All features × all methods × all datasets, with unified schema |
| `significant` | `data.frame` | Subset where `padj < alpha` and `diff_robust != FALSE` |
| `consensus` | `data.frame` | One row per (dataset, feature), summarising cross-method agreement |

### Unified result schema (`all_results` / `significant`)

| Column | Type | Description |
|---|---|---|
| `feature` | character | Taxon / pathway identifier |
| `log2FC` | numeric | Log2 fold-change (or CLR estimate for LMM) |
| `pval` | numeric | Nominal p-value |
| `padj` | numeric | Adjusted p-value (method specified by `p_adj_method`) |
| `diff_robust` | logical / NA | ANCOMBC2 robustness flag; `NA` for methods that do not provide it |
| `method` | character | One of `LinDA`, `ANCOMBC2`, `LMM`, `edgeR`, `metagenomeSeq` |
| `dataset` | character | Dataset name as provided in `ps_list` |

### Consensus table schema

| Column | Type | Description |
|---|---|---|
| `dataset` | character | Dataset name |
| `feature` | character | Taxon / pathway identifier |
| `n_methods_sig` | integer | Number of methods that flagged this feature as significant |
| `methods_agreed` | character | Comma-separated list of agreeing methods |
| `mean_log2FC` | numeric | Mean log2FC across agreeing methods |
| `sd_log2FC` | numeric | Standard deviation of log2FC across methods |
| `min_padj` | numeric | Lowest adjusted p-value across methods |

### CSV files written to `outdir`

| File | Description |
|---|---|
| `<dataset>_all_methods.csv` | Full results for each dataset (one file per element in `ps_list`) |
| `significant_all_methods.csv` | All significant features across datasets and methods |
| `consensus_features.csv` | Cross-method consensus summary |


---

## Note

- `log2FC` from LMM reflects the **CLR-scale regression coefficient** for the variable of interest, not a true log2 fold-change. Interpret accordingly when comparing across methods.
