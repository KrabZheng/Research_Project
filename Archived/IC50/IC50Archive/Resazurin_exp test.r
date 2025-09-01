# IC50 analysis from sequential CSV with smooth dose–response plotting (v6)
# ============================================================================
# - No disk output → IC50 (`all_ic50`) and finely predicted curves (`all_preds`) in R.
# - Robust fitting: skips sparse curves (<3 non-zero doses).
# - CSE normalization: blank = row A, max = row H per BaseID (so B1 borrows B2, B3 borrows B4), ensuring distinct min/max per compound.
# - Viability raw = 100 − normalized signal, clamped 0–100%.
# - Smooth curve plotting:
#     • plot_across(): same compound across patients & conditions
#     • plot_within(): same compound within each patient across conditions
# ---------------------------------------------------------------------------

# Paths ─────────────────────────────────────────────────────────────────────
big_csv_path <- "IC50/Resazurin_Exp_30/20250704_Resazurin exp 30_CSE_H2O2_KBrO3_Menadione_RSL3_RSL3CSE_Erastin_ErastinCSE.csv"
layout_path  <- "IC50/Resazurin_Exp_30/Resazurin_exp_30_corrected.txt"

# Libraries ─────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(readr); library(tidyr)
  library(purrr);  library(stringr); library(tibble)
  library(drc);    library(ggplot2); library(dplyr)
})

# Helpers ───────────────────────────────────────────────────────────────────
safe_fit <- function(df) {
  doses <- unique(df$Concentration[df$Concentration > 0 & !is.na(df$Concentration)])
  if (length(doses) < 3) return(NULL)
  tryCatch(drm(Viability ~ Concentration, data = df, fct = LL.4()), error = function(e) NULL)
}

safe_params <- function(fit) {
  if (is.null(fit)) return(tibble(slope=NA, lower=NA, upper=NA, ic50=NA))
  cf <- coef(fit); tibble(slope=cf[1], lower=cf[2], upper=cf[3], ic50=cf[4])
}

# Load layout metadata ──────────────────────────────────────────────────────
layout_long <- read_tsv(layout_path, show_col_types=FALSE) %>%
  filter(str_detect(str_to_lower(Type), "treatment")) %>%
  pivot_longer(cols = `1`:`8`, names_to = "RowIdx", values_to = "Concentration") %>%
  mutate(Row = LETTERS[as.integer(RowIdx)]) %>%
  select(CompoundID = ID, CompoundName = Name, Row, Concentration, Unit)

# CSE mapping ───────────────────────────────────────────────────────────────
base_map <- c(B1 = "B2", B3 = "B4")

# Extract viability data ────────────────────────────────────────────────────
get_viability_data <- function(rows, layout, base_map) {
  map_dfr(rows, function(i) {
    block <- all_lines[i:(i+8)]
    plate <- suppressMessages(read_csv(I(paste(block, collapse="\n")), show_col_types=FALSE))
    parts <- str_split(names(plate)[1], "_", simplify=TRUE)
    patient <- parts[1]; group <- parts[2]; PlateID <- names(plate)[1]
    sig <- plate[,3:14]; colnames(sig) <- sprintf("%02d", 1:12)
    df <- sig %>%
      mutate(Row = LETTERS[1:8]) %>%
      pivot_longer(-Row, names_to = "Col", values_to = "Signal", names_transform=list(Col=as.integer)) %>%
      mutate(
        PlateID    = PlateID,
        Patient    = patient,
        Group      = group,
        ColGroup   = case_when(Col<=3~1, Col<=6~2, Col<=9~3, TRUE~4),
        CompoundID = paste0(Group, ColGroup),
        BaseID     = if_else(CompoundID %in% names(base_map), base_map[CompoundID], CompoundID)
      ) %>%
      left_join(layout, by=c("CompoundID","Row"))
    norm <- df %>%
      group_by(PlateID, BaseID) %>%
      summarise(
        Blank = mean(Signal[Row=="A"], na.rm=TRUE),
        Max   = mean(Signal[Row=="H"], na.rm=TRUE),
        .groups = "drop"
      )
    df2 <- df %>% left_join(norm, by=c("PlateID","BaseID")) %>%
      mutate(
        Signal_bc     = Signal - Blank,
        Signal_norm   = Signal_bc / (Max - Blank) * 100,
        Viability_raw = 100 - Signal_norm,
        Viability     = pmin(pmax(Viability_raw, 0), 100)
      ) %>%
      filter(!is.na(Concentration) & Concentration > 0)
    df2
  })
}

# IC50 computation ─────────────────────────────────────────────────────────
compute_ic50 <- function(viab) {
  viab %>%
    group_by(Patient,Group,CompoundID,CompoundName,Unit) %>%
    group_modify(~ safe_params(safe_fit(.x))) %>%
    ungroup() %>%
    rename(slope=slope, lower=lower, upper=upper, IC50=ic50)
}

# Predict smooth curves ─────────────────────────────────────────────────────
predict_curves <- function(viab, n=100) {
  viab %>%
    group_by(Patient,Group,CompoundID,CompoundName,Unit) %>%
    nest() %>%
    mutate(
      fit  = map(data, ~ safe_fit(.x)),
      pred = map2(fit,data,~{
        if (is.null(.x)) return(NULL)
        concs <- unique(.y$Concentration); concs <- concs[concs>0]
        if (length(concs) < 2) return(NULL)
        grid <- suppressWarnings(exp(seq(log(min(concs)), log(max(concs)), length.out=n)))
        tibble(Concentration=grid, Viability=predict(.x, newdata=data.frame(Concentration=grid)))
      })
    ) %>% select(-data,-fit) %>% unnest(pred)
}

# Plotting ─────────────────────────────────────────────────────────────────
plot_across <- function(preds, cmpd) {
  df <- filter(preds, CompoundID == cmpd)
  ggplot(df, aes(Concentration, Viability, color = paste(Patient, Group))) +
    geom_line() +
    scale_x_log10() +
    theme_minimal() +
    labs(
      title = sprintf("%s across patients & conditions", unique(df$CompoundName)),
      x     = sprintf("Conc [%s]", unique(df$Unit)),
      y     = "Viability (%)"
    )
}

plot_within <- function(preds, cmpd) {
  df <- filter(preds, CompoundID == cmpd)
  ggplot(df, aes(Concentration, Viability, color = Group)) +
    geom_line() +
    scale_x_log10() +
    facet_wrap(~ Patient) +
    theme_minimal() +
    labs(
      title = sprintf("%s by patient conditions", unique(df$CompoundName)),
      x     = sprintf("Conc [%s]", unique(df$Unit)),
      y     = "Viability (%)"
    )
}

# Run pipeline ────────────────────────────────────────────────────────────
all_lines    <- read_lines(big_csv_path)
rows         <- which(str_detect(all_lines, "_fluo")); rows <- rows[rows+8<=length(all_lines)]
all_viability<- get_viability_data(rows, layout_long, base_map)
all_ic50     <- compute_ic50(all_viability)
all_preds    <- predict_curves(all_viability)

compound_to_view <- "A3"
# Find corresponding ID if user prefers name lookup
if (compound_to_view %in% layout_long$CompoundName) {
  cmpd_id <- layout_long$CompoundID[layout_long$CompoundName == compound_to_view][1]
} else {
  cmpd_id <- compound_to_view
}
# Generate and display
print(plot_across(all_preds, cmpd_id))
print(plot_within(all_preds, cmpd_id))