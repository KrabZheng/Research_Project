#!/usr/bin/env Rscript
# ddcq_pipeline_fixed.R -----------------------------------------------------
# Calculate ΔΔCq and plot log2 fold‑change for qPCR data
#
# Updates (2025‑07‑17, 4th rev):
#   • Re‑implements error bars via pre‑computed summary statistics (mean ± SEM)
#     instead of ggplot2::stat_summary(), eliminating the “attempt to apply
#     non‑function” issue some users encounter when helper objects are masked.
#   • Uses *base::mean* explicitly to avoid any shadowing.
#   • Keeps SVG output and ≥3‑replicate outlier trimming.
# --------------------------------------------------------------------------

library(tidyverse)

# ----------------------------- User settings -----------------------------
root_dir  <- "/Users/SofTicE/Desktop/Research Project/qPCR/Mi/"
well_file <- file.path(root_dir, "WellResultOmit.csv")     # WellResultOmit.csv
info_file <- file.path(root_dir, "qPCR_sample_info.tsv")   # qPCR_sample_info.tsv
out_svg1  <- file.path(root_dir, "ddCq_log2FC_by_gene.svg")          # all replicates
out_svg2  <- file.path(root_dir, "ddCq_log2FC_by_gene_no_outlier.svg") # trimmed

# ------------------------------- Load data -------------------------------
well <- read_csv(well_file, show_col_types = FALSE)
info <- read_tsv(info_file, show_col_types = FALSE)

# ------------------------- Clean Cq & metadata ---------------------------
well_clean <- well %>%
  mutate(Cq = na_if(Cq, "Undetermined"),   # convert “Undetermined” → NA
         Cq = as.numeric(Cq)) %>%            # ensure numeric
  filter(!is.na(Cq))                         # drop omitted wells

# Split metadata
sample_meta <- info %>%
  filter(Info_type == "Sample") %>%
  select(Sample = Info, Condition)

target_meta <- info %>%
  filter(Info_type == "Target") %>%
  select(Target = Info_name,
         Info_key,
         Info_key_1)

# Reference (housekeeping) genes
reference_genes <- target_meta %>%
  filter(Info_key  == "Reference" |
         Info_key_1 == "Reference") %>%
  distinct(Target) %>%
  pull()

# Helper for SEM -----------------------------------------------------------
sem <- function(x) stats::sd(x) / sqrt(length(x))

# ------------- Function to build ΔΔCq results & summary ------------------
make_results <- function(df_replicate) {
  # tech‑rep means
  sample_target_means <- df_replicate %>%
    group_by(Sample, Target) %>%
    summarise(Cq = base::mean(Cq), .groups = "drop")

  # sample‑wise reference mean
  ref_means <- sample_target_means %>%
    filter(Target %in% reference_genes) %>%
    group_by(Sample) %>%
    summarise(ref_mean = base::mean(Cq), .groups = "drop")

  # ΔCq
  dcq_tbl <- sample_target_means %>%
    inner_join(ref_means, by = "Sample") %>%
    mutate(dCq = Cq - ref_mean) %>%
    inner_join(sample_meta, by = "Sample")

  # Day0 baseline per gene
  baseline_tbl <- dcq_tbl %>%
    filter(Condition == "Day0") %>%
    group_by(Target) %>%
    summarise(baseline_dCq = base::mean(dCq), .groups = "drop")

  # ΔΔCq & log2FC
  results <- dcq_tbl %>%
    left_join(baseline_tbl, by = "Target") %>%
    mutate(ddCq   = dCq - baseline_dCq,
           log2FC = -ddCq)

  # Summary stats for error bars (mean ± SEM)
  summary_tbl <- results %>%
    group_by(Target, Condition) %>%
    summarise(mean_log2FC = base::mean(log2FC),
              sem_log2FC  = sem(log2FC),
              .groups = "drop")

  list(replicates = results, summary = summary_tbl)
}

# --------------------------- No‑trim pipeline ----------------------------
res_all <- make_results(well_clean)

plot_ddcq <- ggplot() +
  geom_point(data = res_all$replicates,
             aes(x = Condition, y = log2FC, group = Sample),
             alpha = 0.7,
             position = position_dodge(width = 0.25)) +
  geom_line(data = res_all$summary,
            aes(x = Condition, y = mean_log2FC, group = 1),
            linewidth = 1.1) +
  geom_errorbar(data = res_all$summary,
                aes(x = Condition,
                    ymin = mean_log2FC - sem_log2FC,
                    ymax = mean_log2FC + sem_log2FC),
                width = 0.25) +
  facet_wrap(~ Target, scales = "free_y") +
  labs(title = "Log₂ fold‑change (ΔΔCq) relative to Day 0",
       y = "log₂(FC)",
       x = NULL) +
  theme_bw() +
  theme(legend.position = "none")

ggsave(out_svg1, plot_ddcq, width = 9, height = 6)

# ------------------------- Outlier‑trim pipeline -------------------------
# Attach Condition to replicate‑level data
well_cond <- well_clean %>%
  left_join(sample_meta, by = "Sample")

# Helper: drop one furthest Cq if n ≥ 3 replicates
trim_furthest <- function(df) {
  if (nrow(df) < 3) return(df)
  med <- stats::median(df$Cq)
  df %>%
    mutate(dist = abs(Cq - med)) %>%
    slice(-which.max(dist)) %>%
    select(-dist)
}

well_trim <- well_cond %>%
  group_by(Sample, Target, Condition) %>%
  group_modify(~ trim_furthest(.x)) %>%
  ungroup()

res_trim <- make_results(well_trim)

plot_ddcq_trim <- ggplot() +
  geom_point(data = res_trim$replicates,
             aes(x = Condition, y = log2FC, group = Sample),
             alpha = 0.7,
             position = position_dodge(width = 0.25)) +
  geom_line(data = res_trim$summary,
            aes(x = Condition, y = mean_log2FC, group = 1),
            linewidth = 1.1) +
  geom_errorbar(data = res_trim$summary,
                aes(x = Condition,
                    ymin = mean_log2FC - sem_log2FC,
                    ymax = mean_log2FC + sem_log2FC),
                width = 0.25) +
  facet_wrap(~ Target, scales = "free_y") +
  labs(title = "Log₂ fold‑change (ΔΔCq) – outlier trimmed (≥3 reps)",
       y = "log₂(FC)",
       x = NULL) +
  theme_bw() +
  theme(legend.position = "none")

ggsave(out_svg2, plot_ddcq_trim, width = 9, height = 6)

# --------------------------------------------------------------------------
