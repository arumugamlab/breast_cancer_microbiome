################################################################################
##  Sutterellaceae Relapse & Survival Analysis
##  Author : Camila Alvarez-Silva
##  Project: Gut microbiome and breast cancer recurrence
##
##  Sections:
##    0. Dependencies
##    1. Load & prepare survival metadata
##    2. Load & filter phyloseq object (breast cancer, complete cases)
##    3. OTU / CLR preparation
##    4. Logistic regression — odds ratios (OR)
##       4a. Binary presence/absence model
##       4b. Continuous CLR-abundance model
##    5. Survival analysis
##       5a. Build Surv object and annotate analysis data frame
##       5b. Cox proportional hazards — binary model
##       5c. Cox proportional hazards — CLR model
##       5d. Kaplan–Meier plot
################################################################################


# ── 0. Dependencies ────────────────────────────────────────────────────────────

library(readxl)
library(dplyr)
library(phyloseq)
library(compositions)   # clr()
library(survival)
library(survminer)      # ggsurvplot()


# ── 1. Load & prepare survival metadata ───────────────────────────────────────

# Load the clinical data  to match participant IDs to DNA_IDs 
clinical_data.filt<-readRDS("/projects/dk_breast_cancer-AUDIT/data/processed/metadata/Clean/Clinical_data.2026_updated.rds")
row.names(clinical_data.filt)<-clinical_data.filt$DNA_ID

# Load the survival data sheet and match participant IDs to DNA_IDs via lookup
SurvivalAnalysis <- read_excel("/maps/projects/arumugam/people/xhf865/Camila/BreastCancer/data/SurvivalAnalysis.xlsx") %>% as.data.frame()

# Map participant_id → DNA_ID using the clinical data reference table
SurvivalAnalysis$DNA_ID <- lookup_e(
  paste0("BC", SurvivalAnalysis$participant_id),   # formatted participant ID
  as.character(clinical_data.filt$participant_id),
  as.character(clinical_data.filt$DNA_ID)
)

# Keep only rows with a matched DNA_ID
SurvivalAnalysis <- filter(SurvivalAnalysis, !is.na(DNA_ID))
row.names(SurvivalAnalysis) <- as.character(SurvivalAnalysis$DNA_ID)

# Keep only breast cancer patients with a known relapse status
SurvivalAnalysis <- filter(SurvivalAnalysis, !is.na(relapse))

# Parse key date columns
SurvivalAnalysis$relapse_date  <- as.Date(SurvivalAnalysis$relapse_date,  format = "%Y-%m-%d")
SurvivalAnalysis$baseline_date <- as.Date(SurvivalAnalysis$enrolment,     format = "%Y-%m-%d")
SurvivalAnalysis$last_lookup   <- as.Date(SurvivalAnalysis$date_lookup,   format = "%Y-%m-%d")


# ── 2. Load & filter phyloseq (breast cancer, complete cases) ─────────────────

ps.motus <- readRDS("/projects/arumugam/people/xhf865/Camila/BreastCancer/phyloseq/ps.motus_raw.3.0.3.rds")

# Subset to breast cancer samples only
ps.motus <- prune_samples(
  sample_data(ps.motus)$DiseaseStatus == "BreastCancer",
  ps.motus
)

# Remove samples with missing values in any covariate used in modelling
ps.motus <- subset_samples(
  ps.motus,
  !is.na(relapse)                    &
    !is.na(bmi_calculated)           &
    !is.na(ki67_ihc)                 &
    !is.na(ln_tumorpositive.binary)  
)

# Drop taxa with zero total counts after sample filtering
ps.motus <- prune_taxa(taxa_sums(ps.motus) > 0, ps.motus)

# Extract OTU matrix (features × samples) and metadata
otu  <- as(otu_table(ps.motus), "matrix")
meta <- as(sample_data(ps.motus), "data.frame")

# Prevalence filter: retain taxa present in ≥ 10 % of samples
prev_thresh <- 0.10
keep        <- which(rowMeans(otu > 0) >= prev_thresh)
otu         <- otu[keep, , drop = FALSE]
ps.motus    <- prune_taxa(taxa_names(ps.motus) %in% names(keep), ps.motus)


# ── 3. OTU / CLR data preparation ─────────────────────────────────────────────

# Transpose OTU matrix to samples × features for modelling
otu_mat <- as.data.frame(t(otu_table(ps.motus)))

# Align metadata row order with OTU matrix
meta.df           <- as.data.frame(unclass(sample_data(ps.motus)))
row.names(meta.df) <- meta.df$DNA_ID
meta.df           <- meta.df[rownames(otu_mat), ]
stopifnot(all(rownames(meta.df) == rownames(otu_mat)))

# CLR transformation (centred log-ratio) with half-minimum pseudocount for zeros
min_val    <- min(otu_mat[otu_mat > 0])
pseudocount <- min_val / 2

otu_pos <- otu_mat
otu_pos[otu_pos == 0] <- pseudocount

otu_clr <- as.data.frame(apply(otu_pos, 2, compositions::clr))

# Align metadata to CLR matrix
meta.df <- meta.df[rownames(otu_clr), ]
stopifnot(all(rownames(meta.df) == rownames(otu_clr)))


# ── 4. Logistic regression — odds ratios ──────────────────────────────────────
#
#  Target taxon: Sutterellaceae (mOTU ID: meta_mOTU_v3_12389)
#  Outcome     : relapse.binary  (1 = relapsed, 0 = no relapse)
#  Covariates  : ki67_ihc, bmi_calculated, ln_tumorpositive.binary

## 4a. Binary model — Sutterellaceae presence / absence ────────────────────────

df_bin <- meta.df %>%
  mutate(Sutterellaceae_bin = otu_mat[, "meta_mOTU_v3_12389"] > 0)

table(df_bin$Sutterellaceae_bin)   # check prevalence in cohort

model_bin <- glm(
  relapse.binary ~ Sutterellaceae_bin + ki67_ihc + bmi_calculated + ln_tumorpositive.binary,
  data   = df_bin,
  family = binomial
)

# Exponentiated coefficients = odds ratios with 95 % CI
or_bin <- exp(cbind(OR = coef(model_bin), confint(model_bin)))
print(or_bin)


## 4b. Continuous model — Sutterellaceae CLR abundance ─────────────────────────

df_clr <- meta.df %>%
  mutate(Sutterellaceae = otu_clr$meta_mOTU_v3_12389)

model_clr <- glm(
  relapse.binary ~ Sutterellaceae + ki67_ihc + bmi_calculated + ln_tumorpositive.binary,
  data   = df_clr,
  family = binomial
)

or_clr <- exp(cbind(OR = coef(model_clr), confint(model_clr)))
print(or_clr)


# ── 5. Survival analysis ───────────────────────────────────────────────────────

## 5a. Build survival data frame ────────────────────────────────────────────────

# Event indicator: 1 = relapsed (from the relapse column), 0 = censored
SurvivalAnalysis$relapse_event <- SurvivalAnalysis$relapse

# Time-to-event (days):
#   • Relapsed  → baseline to relapse date
#   • Censored  → baseline to last follow-up date
SurvivalAnalysis$time_to_relapse <- ifelse(
  SurvivalAnalysis$relapse_event == 1,
  as.numeric(SurvivalAnalysis$relapse_date  - SurvivalAnalysis$baseline_date),
  as.numeric(SurvivalAnalysis$last_lookup   - SurvivalAnalysis$baseline_date)
)

# Align row order across all three data frames before merging columns
SurvivalAnalysis <- SurvivalAnalysis[rownames(otu_clr), ]
SurvivalAnalysis <- SurvivalAnalysis[rownames(meta.df),  ]
stopifnot(
  all(rownames(SurvivalAnalysis) == rownames(otu_clr)),
  all(rownames(SurvivalAnalysis) == rownames(meta.df))
)

# Attach microbial and clinical variables needed for Cox models
SurvivalAnalysis <- SurvivalAnalysis %>%
  mutate(
    Sutterellaceae          = otu_clr$meta_mOTU_v3_12389,  # CLR-transformed abundance
    Sutterellaceae_bin      = otu_mat[, "meta_mOTU_v3_12389"] > 0,  # presence/absence
    ki67_ihc                = meta.df$ki67_ihc,
    bmi_calculated          = meta.df$bmi_calculated,
    ln_tumorpositive.binary = meta.df$ln_tumorpositive.binary,
    tumor_size              = meta.df$tumor_size
  )

# Survival object
surv_obj <- Surv(
  time  = SurvivalAnalysis$time_to_relapse,
  event = SurvivalAnalysis$relapse_event
)


## 5b. Cox model 1 — Sutterellaceae presence / absence ─────────────────────────

cox_bin <- coxph(
  surv_obj ~ Sutterellaceae_bin + ki67_ihc + bmi_calculated + ln_tumorpositive.binary,
  data = SurvivalAnalysis
)
summary(cox_bin)


## 5c. Cox model 2 — Sutterellaceae CLR abundance ──────────────────────────────

cox_clr <- coxph(
  surv_obj ~ Sutterellaceae + ki67_ihc + bmi_calculated + ln_tumorpositive.binary,
  data = SurvivalAnalysis
)
summary(cox_clr)


## 5d. Kaplan–Meier plot ────────────────────────────────────────────────────────

fit_km <- survfit(surv_obj ~ Sutterellaceae_bin, data = SurvivalAnalysis)

km_plot <- ggsurvplot(
  fit_km,
  data       = SurvivalAnalysis,
  pval       = TRUE,               # log-rank p-value on plot
  risk.table = FALSE,
  palette    = c("#264653", "#9d0208"),   # absent = dark teal, present = deep red
  legend.labs = c("Absent", "Present"),
  title       = "Recurrence-free survival by Sutterellaceae presence"
)

print(km_plot)
