# --- Load libraries ---
library(tidyverse)
library(RColorBrewer)
# --- Paths ---
data_dir   <- "HNF1B"
well_file <- file.path(data_dir, "WellResult.csv")
info_file <- file.path(data_dir, "qPCR_sample_info.tsv")
# --- Read and filter ---
well <- read_csv(well_file) |>
  filter(!Omit)
info <- read_tsv(info_file)
# --- Mapping tables ---
sample_info <- info |>
  filter(Info_type == "Sample") |>
  select(Plate_Sample = Info, Patient, Condition)
target_info <- info |>
  filter(Info_type == "Target") |>
  select(Plate_Target = Info, Gene = Info_name, Role = Info_key)
# --- Summarise Cq and drop one outlier if >=3 reps ---
cq_summary <- well |>
  left_join(sample_info, by = c("Sample" = "Plate_Sample")) |>
  left_join(target_info, by = c("Target" = "Plate_Target")) |>
  transmute(Patient, Condition, Gene, Role, Cq) |>
  group_by(Patient, Condition, Gene, Role) |>
  filter(
    !(n() >= 3 &
        abs(Cq - median(Cq, na.rm = TRUE)) ==
          max(abs(Cq - median(Cq, na.rm = TRUE)), na.rm = TRUE))
  ) |>
  summarise(
    mean_Cq = mean(Cq, na.rm = TRUE),
    sd_Cq   = sd(Cq,   na.rm = TRUE),
    n       = sum(!is.na(Cq)),
    .groups = "drop"
  )
# --- Split into reference and target, compute ΔCq & ΔΔCq ---
ref_cq <- cq_summary |>
  filter(Role == "Reference") |>
  select(Patient, Condition, ref_Cq = mean_Cq, sd_ref = sd_Cq, n_ref = n)
target_cq <- cq_summary |>
  filter(Role != "Reference") |>
  select(Patient, Condition, Gene,
         targ_Cq = mean_Cq, sd_targ = sd_Cq, n_targ = n)
deltacq <- target_cq |>
  left_join(ref_cq, by = c("Patient", "Condition")) |>
  mutate(
    delta_Cq = targ_Cq - ref_Cq,
    se_delta = sqrt(sd_targ^2 + sd_ref^2) /
      sqrt(pmin(n_targ, n_ref, na.rm = TRUE))
  ) |>
  group_by(Patient, Gene) |>
  mutate(
    control_delta = delta_Cq[Condition == "Control"],
    se_control    = se_delta[Condition == "Control"]
  ) |>
  ungroup() |>
  filter(Condition != "Control") |>
  mutate(
    delta_delta_Cq = delta_Cq - control_delta,
    se_ddCt        = sqrt(se_delta^2 + se_control^2),
    log2FC         = -delta_delta_Cq
  )
# --- Greyscale check helper ---
check_palette_greyscale <- function(cols, threshold = 30) {
  lums <- col2rgb(cols) |>
    apply(2, function(rgb) sum(rgb * c(0.299, 0.587, 0.114)))
  if (any(combn(lums, 2, function(x) abs(diff(x))) < threshold)) {
    warning("Palette colours differ by <", threshold, " in greyscale")
  }
}
# --- Plotting function ---
plot_log2fc <- function(df, palette = "YlGnBu", threshold = 30) {
  df <- df |> filter(!is.na(log2FC))
  if (nrow(df) == 0 || all(df$log2FC == 0)) return(NULL)
  cols <- brewer.pal(n_distinct(df$Condition), palette)
  check_palette_greyscale(cols, threshold)
  ggplot(df, aes(Condition, log2FC, fill = Condition)) +
    geom_hline(yintercept = 0) +
    geom_col(position = position_dodge(0.8), width = 0.7, color = "black",
             size = 0.3, show.legend = FALSE) +
    geom_errorbar(aes(ymin = log2FC - se_ddCt, ymax = log2FC + se_ddCt),
                  position = position_dodge(0.8), width = 0.2) +
    scale_fill_manual(values = cols) +
    facet_grid(Gene ~ Patient, scales = "fixed", space = "fixed") +
    coord_fixed(ratio = 1) +
    theme_bw(base_size = 12) +
    labs(x = "Condition", y = "Log₂ Fold-Change vs. Control") +
    theme(
      aspect.ratio     = 1,
      strip.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x      = element_text(angle = 45, hjust = 1)
    )
}
# --- Generate and save plot ---
p <- plot_log2fc(deltacq)
if (!is.null(p)) {
  ggsave("qpcr_log2FC.svg", plot = p, device = "svg",
         width = 2048, height = 2048, units = "px", dpi = 300)
}