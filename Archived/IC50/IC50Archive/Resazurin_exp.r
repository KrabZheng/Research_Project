# Drug Response Curve Analysis Script (Final)
# ================================================
# Imports 96-well plate fluorescence data, merges metadata, normalizes signals,
# fits dose–response curves, computes IC50, and plots individual smooth curves.

# 1. Load libraries explicitly
library(readr)    # read_csv, read_tsv
library(tidyr)    # pivot_longer, unnest
library(stringr)  # str_split
library(ggplot2)  # plotting
library(drc)      # dose–response modeling
library(purrr)    # map functions
library(tibble)   # tibble
library(dplyr)    # data manipulation

# 2. Define file paths (adjust accordingly)
raw_data_file   <- "/Users/SofTicE/Desktop/Research Project/IC50/P676_A_fluo.csv"
info_sheet_file <- "/Users/SofTicE/Desktop/Research Project/IC50/Resazurin_exp_30.tsv"

# 3. Read raw plate data and extract plate metadata
raw_raw     <- read_csv(raw_data_file, col_names = FALSE)
plate_name  <- raw_raw[[1,1]]
group_letter <- str_split(plate_name, "_")[[1]][2]
# Column headers (wells 1–12)
header_vals <- raw_raw[1, 3:14] %>% as.character()

# 4. Reshape raw signals to long format and assign CompoundID
df_signals <- raw_raw %>%
  slice(2:9) %>%                                      # rows A–H
  rename_with(~ header_vals, .cols = 3:14) %>%
  mutate(Row = LETTERS[row_number()]) %>%             # A–H
  pivot_longer(
    cols      = all_of(header_vals),
    names_to  = "Col",
    values_to = "Signal",
    names_transform = list(Col = as.integer)
  ) %>%
  mutate(
    Signal     = as.numeric(Signal),
    Well       = paste0(Row, Col),
    ColGroup   = case_when(
      Col %in% 1:3   ~ 1,
      Col %in% 4:6   ~ 2,
      Col %in% 7:9   ~ 3,
      Col %in% 10:12 ~ 4
    ),
    CompoundID = paste0(group_letter, ColGroup)
  ) %>%
  select(Well, Row, Col, Signal, CompoundID)

# 5. Read metadata and reshape for plate compounds
info_df   <- read_tsv(info_sheet_file)
valid_ids <- paste0(group_letter, 1:4)
info_long <- info_df %>%
  filter(ID %in% valid_ids) %>%
  pivot_longer(
    cols     = `1`:`8`,
    names_to = "RowIndex",
    values_to= "Concentration"
  ) %>%
  mutate(
    Row          = LETTERS[as.integer(RowIndex)],
    CompoundID   = ID,
    CompoundName = Name,
    Unit         = Unit
  ) %>%
  select(CompoundID, CompoundName, Row, Concentration, Unit)

# 6. Merge signals with metadata and normalize per compound
df <- df_signals %>%
  left_join(info_long, by = c("CompoundID","Row"))
norm_stats <- df %>%
  group_by(CompoundID) %>%
  summarize(
    Blank     = mean(Signal[Concentration == 0], na.rm = TRUE),
    MaxSignal = mean(Signal[Row == "H"], na.rm = TRUE)
  )
df <- df %>%
  left_join(norm_stats, by = "CompoundID") %>%
  mutate(
    Signal_bc   = Signal - Blank,
    Signal_norm = Signal_bc / (MaxSignal - Blank) * 100
  )

# 7. Fit dose–response models excluding zero concentration
models <- df %>%
  filter(!is.na(Concentration) & Concentration > 0) %>%
  group_by(CompoundID, CompoundName) %>%
  nest() %>%
  mutate(
    fit = map(data, ~ suppressWarnings(
      drm(
        Signal_norm ~ Concentration,
        data = .,
        fct = LL.4(names = c("Slope","Lower","Upper","IC50"))
      )
    ))
  )

# 8. Extract IC50 values with robust error handling
ic50_df <- models %>%
  mutate(
    IC50 = map_dbl(fit, ~ {
      ed <- tryCatch(ED(.x, 50, interval = "delta"), error = function(e) NA)
      if (is.matrix(ed)) as.numeric(ed[1, "Estimate"]) else if (is.numeric(ed)) ed[1] else NA_real_
    })
  ) %>%
  select(CompoundID, CompoundName, IC50)

# 9. Generate smooth prediction curves manually using fitted parameters
df_pred <- models %>%
  filter(map_lgl(data, ~ {
    cs <- .x$Concentration; cs <- cs[!is.na(cs) & cs > 0]; length(unique(cs)) > 1
  })) %>%
  mutate(
    pred = map2(data, fit, ~ {
      cs <- .x$Concentration; cs <- cs[!is.na(cs) & cs > 0]
      conc_grid <- exp(seq(log(min(cs)), log(max(cs)), length.out = 200))
      cf <- as.numeric(coef(.y))
      slope <- cf[1]; lower <- cf[2]; upper <- cf[3]; ic50 <- cf[4]
      sig <- lower + (upper - lower) / (1 + (conc_grid / ic50)^slope)
      tibble(Concentration = conc_grid, Signal_norm = sig)
    })
  ) %>%
  select(CompoundID, CompoundName, pred) %>%
  unnest(cols = pred)

# 10. Plot individual smooth curves with viability and per-plot axes
plots <- models %>%
  select(CompoundID, CompoundName) %>%
  distinct() %>%
  mutate(
    p = map2(CompoundID, CompoundName, ~ {
      df_data <- df %>% filter(CompoundID == .x & Concentration > 0)
      df_fit  <- df_pred %>% filter(CompoundID == .x)
      # compute viability
      df_data <- df_data %>% mutate(Viability = 100 - Signal_norm)
      df_fit  <- df_fit  %>% mutate(Viability = 100 - Signal_norm)
      # determine x-axis per plot
      x_min <- min(df_data$Concentration, na.rm = TRUE)
      x_max <- max(df_data$Concentration, na.rm = TRUE)
      ggplot(df_data, aes(x = Concentration, y = Viability)) +
        geom_point() +
        geom_line(data = df_fit, aes(x = Concentration, y = Viability), color = "blue") +
        scale_x_log10(limits = c(x_min, x_max), expand = expansion(mult = c(0.05, 0.05))) +
        scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0.05, 0.05)), oob = scales::oob_keep) +
        labs(
          title = paste("Dose–Response for", .y),
          x = paste0("Concentration [", unique(df_data$Unit), "]"),
          y = "Viability (%)"
        ) +
        theme_minimal()
    })
  )
# Display each plot
walk(plots$p, print)

# 11. Print IC50 table
print(ic50_df)