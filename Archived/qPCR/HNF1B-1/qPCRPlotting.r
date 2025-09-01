# 0) Set working directory and global variables
setwd("/Users/SofTicE/Desktop/Research Project/qPCR/HNF1B-1") # path to your files, one-off

# 1) Load Libraries
library(readr)
library(dplyr)
library(ggplot2)

# 2) Read
well <- read_csv("WellResult.csv") # your qPCR output
info <- read_tsv("qPCR_sample_info.tsv") # your plate map

# 3) Remove Omit
well <- well %>% filter(!Omit)

# 4) Mapping Tables
sample_info <- info %>%
  filter(Info_type == "Sample") %>%
  select(Plate_Sample = Info,
         SampleName = Info_name, Patient, Condition)

target_info <- info %>%
  filter(Info_type == "Target") %>%
  select(Plate_Target = Info, Gene = Info_name, Role = Info_key)

# 5) Join mappings onto the well results
dat <- well %>%
  left_join(sample_info, by = c("Sample" = "Plate_Sample")) %>%
  left_join(target_info, by = c("Target" = "Plate_Target")) %>%
  select(Patient, Condition, Gene, Role, Cq)

# 5.5) Remove furthest outlier per Patient×Condition×Gene×Role when n ≥ 3
dat <- dat %>%
  group_by(Patient, Condition, Gene, Role) %>%
  # compute deviation from median, and count replicates
  mutate(
    dev     = abs(Cq - median(Cq, na.rm = TRUE)),
    n_rep   = sum(!is.na(Cq))
  ) %>%
  # if at least 3 reps, drop the one with max dev; otherwise keep all
  filter(!(n_rep >= 3 & dev == max(dev))) %>%
  ungroup() %>%
  select(-dev, -n_rep)

# 6) Mean Cq, SD, and counts per Patient×Condition×Gene and Role
mean_cq <- dat %>%
  group_by(Patient, Condition, Gene, Role) %>%
  summarise(
    mean_Cq = mean(Cq, na.rm = TRUE),
    sd_Cq   = sd(Cq,   na.rm = TRUE),
    n       = sum(!is.na(Cq)),
    .groups = "drop"
  )

# 7) Split reference and target gene Cq means
ref_cq <- mean_cq %>%
  filter(Role == "Reference") %>%
  select(Patient, Condition, ref_Cq = mean_Cq, sd_ref = sd_Cq, n_ref = n)

target_cq <- mean_cq %>%
  filter(Role != "Reference") %>%
  select(Patient, Condition, Gene,
         targ_Cq = mean_Cq,
         sd_targ = sd_Cq,
         n_targ = n)

# 8) Calculate ΔCq and its SE
deltacq <- target_cq %>%
  left_join(ref_cq, by = c("Patient", "Condition")) %>%
  mutate(
    delta_Cq    = targ_Cq - ref_Cq,
    sd_delta    = sqrt(sd_targ^2 + sd_ref^2),
    se_delta    = sd_delta / sqrt(pmin(n_targ, n_ref, na.rm = TRUE))
  )

# 9) Extract the control ΔCq and its SE for each patient & gene
control_delta <- deltacq %>%
  filter(Condition == "Control") %>%
  select(Patient, Gene, control_delta = delta_Cq, se_control = se_delta)

# 10) Compute ΔΔCq, propagate error, and log2 fold‐change
results <- deltacq %>%
  left_join(control_delta, by = c("Patient", "Gene")) %>%
  mutate(
    delta_delta_Cq = delta_Cq - control_delta,
    se_ddCt        = sqrt(se_delta^2 + se_control^2),
    log2FC         = -delta_delta_Cq
  )

# helper: compute grey-level (Y) and check spacing
check_palette_greyscale <- function(cols, threshold = 30) {
  # Rec. 601 luma: Y = 0.299 R + 0.587 G + 0.114 B
  lums <- sapply(cols, function(c) {
    rgb <- col2rgb(c)
    0.299 * rgb[1, ] + 0.587 * rgb[2, ] + 0.114 * rgb[3, ]
  })
  diffs <- combn(lums, 2, FUN = function(x) abs(diff(x)))
  if (any(diffs < threshold)) {
    warning(
      "Palette colours differ by < ", threshold,
      " in greyscale—some may be indistinguishable in B/W."
    )
  }
}

# 11) Define plotting function
plot_log2fc <- function(df, brewer_palette = "YlGnBu", grey_thresh = 30) {
  df2 <- df %>%
    filter(!is.na(log2FC) & Condition != "Control") # nolint

  # bail if no data or all-zero
  if (nrow(df2) == 0 || all(df2$log2FC == 0)) return(NULL)

  # determine how many distinct fills we need
  n_fills <- length(unique(df2$Condition))
  # grab that many colours from the chosen palette
  cols <- RColorBrewer::brewer.pal(n = n_fills, name = brewer_palette)

  # check their greyscale separation
  check_palette_greyscale(cols, threshold = grey_thresh) # nolint

  # now use scale_fill_manual instead of scale_fill_brewer
  ggplot(df2, aes(x = Condition, y = log2FC, fill = Condition)) + # nolint
    geom_hline(yintercept = 0, linetype = "solid") +
    geom_col(
      position = position_dodge(width = 0.8),
      width    = 0.7,
      color    = "black",
      size     = 0.3,
      show.legend = FALSE
    ) +
    geom_errorbar(
      aes(ymin = log2FC - se_ddCt, ymax = log2FC + se_ddCt), # nolint
      position = position_dodge(width = 0.8),
      width    = 0.2,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = cols) +
    facet_grid(Gene ~ Patient, scales = "fixed", space = "fixed") +
    coord_fixed(ratio = 1) +
    theme_bw(base_size = 12) +
    labs(
      x = "Condition",
      y = "Log₂ Fold-Change vs. Control"
    ) +
    theme(
      aspect.ratio       = 1,
      strip.background   = element_blank(),
      panel.grid.major   = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.x        = element_text(angle = 45, hjust = 1)
    )
}

# 12) Generate and save plot as SVG
p <- plot_log2fc(results, brewer_palette = "YlGnBu", grey_thresh = 30)
ggsave(
  filename = "qpcr_log2FC.svg",
  plot     = p,
  device   = "svg",
  width    = 2048,
  height   = 2048,
  units    = "px"
)