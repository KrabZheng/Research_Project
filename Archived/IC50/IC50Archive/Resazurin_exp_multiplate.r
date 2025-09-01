# IC50 analysis from one sequential CSV – **interactive version (v2)**
# ============================================================================
# * No disk output* → IC50 results (`all_ic50`) live in your R session.
# * Robust fitting: skips sparse curves (<3 non-zero doses).
# * Custom CSE normalization: B1/B3 borrow blank/max from B2/B4.
# * Viability is defined as _100 − normalized signal_.
# * IC50 table includes slope, lower, upper parameters and concentration unit.
# ---------------------------------------------------------------------------

# Paths ── edit as needed -----------------------------------------------------
big_csv_path <- "/Users/SofTicE/Desktop/Research Project/IC50/20250704_Resazurin exp 30_CSE_H2O2_KBrO3_Menadione_RSL3_RSL3CSE_Erastin_ErastinCSE.csv"
layout_path  <- "/Users/SofTicE/Desktop/Research Project/IC50/Resazurin_exp_30_corrected.txt"

# ---- suppress messages and load libraries ----------------------------------
suppressPackageStartupMessages({
  library(readr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(tibble)
  library(drc)
  library(dplyr)
})

# ---------- helpers ---------------------------------------------------------
safe_fit <- function(df) {
  doses <- unique(df$Concentration[df$Concentration > 0 & !is.na(df$Concentration)])
  if (length(doses) < 3) return(NULL)
  tryCatch(drm(Viability ~ Concentration, data = df, fct = LL.4()), error = function(e) NULL)
}

safe_params <- function(fit_obj) {
  if (is.null(fit_obj)) return(tibble(slope = NA, lower = NA, upper = NA, ic50 = NA))
  cf <- coef(fit_obj)
  tibble(slope = cf[1], lower = cf[2], upper = cf[3], ic50 = cf[4])
}

# ---------- layout metadata -------------------------------------------------
layout_long <- read_tsv(layout_path, show_col_types = FALSE) %>%
  filter(str_detect(str_to_lower(Type), "treatment")) %>%
  pivot_longer(cols = `1`:`8`, names_to = "RowIdx", values_to = "Concentration") %>%
  mutate(Row = LETTERS[as.integer(RowIdx)]) %>%
  select(CompoundID = ID, CompoundName = Name, Row, Concentration, Unit)

# ---------- parse raw export ------------------------------------------------
all_lines   <- read_lines(big_csv_path)
header_rows <- which(str_detect(all_lines, "_fluo"))

# CSE mapping (define once)
base_map <- c(B1 = "B2", B3 = "B4")

process_block <- function(start) {
  block   <- all_lines[start:(start + 8)]
  plate   <- read_csv(I(paste(block, collapse = "\n")), show_col_types = FALSE)
  meta    <- str_split(names(plate)[1], "_", simplify = TRUE)
  patient <- meta[1]; group <- meta[2]

  # extract exactly the 12 fluorescence columns by position
  sig_mat <- plate[, 3:14]
  colnames(sig_mat) <- sprintf("%02d", 1:12)
  df <- sig_mat %>%
    mutate(Row = LETTERS[1:8]) %>%
    pivot_longer(-Row, names_to = "Col", values_to = "Signal", names_transform = list(Col = as.integer)) %>%
    mutate(
      PlateID    = names(plate)[1],
      Patient    = patient,
      Group      = group,
      ColGroup   = case_when(Col <= 3 ~ 1, Col <= 6 ~ 2, Col <= 9 ~ 3, TRUE ~ 4),
      CompoundID = paste0(Group, ColGroup),
      BaseID     = if_else(CompoundID %in% names(base_map), base_map[CompoundID], CompoundID)
    ) %>%
    left_join(layout_long, by = c("CompoundID", "Row"))

  # CSE normalization
  norm_df <- df %>%
    group_by(PlateID, BaseID) %>%
    summarise(
      Blank = mean(Signal[Row == "A"], na.rm = TRUE),
      Max   = mean(Signal[Row == "H"], na.rm = TRUE),
      .groups = "drop"
    )

  df2 <- df %>%
    left_join(norm_df, by = c("PlateID", "BaseID")) %>%
    mutate(
      Signal_bc   = Signal - Blank,
      Signal_norm = Signal_bc / (Max - Blank) * 100,
      Viability   = 100 - Signal_norm
    )

  # compute IC50
  df2 %>%
    filter(!is.na(Concentration) & Concentration > 0) %>%
    group_by(Patient, Group, CompoundID, CompoundName, Unit) %>%
    group_modify(~ safe_params(safe_fit(.x))) %>%
    ungroup() %>%
    rename(
      slope = slope,
      lower = lower,
      upper = upper,
      IC50  = ic50
    )
}

# ---------- run across all plates -------------------------------------------
all_ic50 <- map_dfr(header_rows, process_block)

# Output ---------------------------------------------------------------------
print(all_ic50)
