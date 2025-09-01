# Resazurin Exp30 — per-sample dose–response with drc IC50s + CI ribbons
# Packages ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
  library(drc)
  library(glue)
  library(readr)
  library(scales)
})

# Paths -------------------------------------------------------------------
raw_path <- "Resazurin_IC50/Resazurin_Exp_35_Raw.txt"
cmp_path <- "Resazurin_IC50/Resazurin_Exp_35_Compound.txt"
sample_path <- "Resazurin_IC50/Resazurin_Exp_35_Sample.txt"
out_dir  <- "Resazurin_IC50/resazurin_outputs/Exp34"
dir.create(out_dir, showWarnings = FALSE)
dir.create(file.path(out_dir, "plots"), showWarnings = FALSE)
dir.create(file.path(out_dir, "csv"),   showWarnings = FALSE)

# Helper: parse RAW (Exp34 format: A1/B1 headers; 8 lines; 12 Fluo + 12 Abs + 12 BG per line)
parse_raw_blocks <- function(path) {
  lines <- readLines(path, warn = FALSE)
  hdr_ix <- grep("^([AB]1)_Fluo\\t1\\t2\\t3\\t4\\t5\\t6\\t7\\t8\\t9\\t10\\t11\\t12", lines)
  stopifnot(length(hdr_ix) > 0)

  purrr::map_dfr(hdr_ix, function(i) {
    hdr <- lines[[i]]
    pos <- stringr::str_match(hdr, "^([AB]1)_Fluo")[,2]
    plate <- stringr::str_sub(pos, 1, 1)
    block <- as.integer(stringr::str_sub(pos, 2, 2))  # 1

    block_lines <- lines[(i+1):(i+8)]  # 8 dose levels
    purrr::imap_dfr(block_lines, function(ln, idx) {
      parts <- stringr::str_split(ln, "\\t")[[1]]
      parts <- parts[parts != ""]              # drop blank separators
      lvl   <- suppressWarnings(as.integer(parts[1]))
      if (is.na(lvl)) lvl <- idx               # Exp34 sometimes omits the 1..8 label
      fl    <- suppressWarnings(as.numeric(parts[2:13]))  # 12 Fluo values
      tibble(
        plate, pos, block, level = lvl,
        col = 1:12, meas = "Fluo", value = fl
      )
    })
  })
}

# Map column (1..12) -> sample_id from the Sample file
parse_sample_map <- function(path) {
  # Read everything as character to avoid mixed column types
  smp <- readr::read_tsv(path, col_types = readr::cols(.default = readr::col_character()))

  # Identify grouping columns (everything except ID, if present)
  if ("ID" %in% names(smp)) {
    grp_cols <- setdiff(names(smp), "ID")
  } else {
    smp <- smp %>% dplyr::mutate(ID = dplyr::row_number())
    grp_cols <- setdiff(names(smp), "ID")
  }

  # Preferred: a row whose ID == "Sample" supplies the sample names
  sample_row <- smp %>% dplyr::filter(tolower(ID) == "sample")
  if (nrow(sample_row) == 1) {
    sample_names <- sample_row %>%
      tidyr::pivot_longer(dplyr::all_of(grp_cols), names_to = "grp", values_to = "sample_id") %>%
      dplyr::select(grp, sample_id)
  } else {
    # Fallback: use the column headers themselves as sample IDs
    sample_names <- tibble::tibble(grp = grp_cols, sample_id = grp_cols)
  }

  # Gather flags for columns 1..12 and coerce to logical
  flag_rows <- smp %>%
    dplyr::filter(grepl("^[0-9]+$", ID)) %>%           # rows "1".."12"
    dplyr::mutate(col = as.integer(ID)) %>%
    tidyr::pivot_longer(dplyr::all_of(grp_cols), names_to = "grp", values_to = "flag_chr") %>%
    dplyr::mutate(
      flag_chr = stringr::str_trim(flag_chr),
      flag_num = suppressWarnings(as.numeric(flag_chr)),
      flag = dplyr::case_when(
        is.na(flag_chr)                  ~ FALSE,
        !is.na(flag_num)                 ~ flag_num != 0,
        stringr::str_detect(stringr::str_to_lower(flag_chr), "^(true|t|x|yes|y)$")  ~ TRUE,
        stringr::str_detect(stringr::str_to_lower(flag_chr), "^(false|f|no|n)$")    ~ FALSE,
        TRUE ~ FALSE
      )
    )

  col_map <- flag_rows %>%
    dplyr::filter(flag) %>%
    dplyr::left_join(sample_names, by = "grp") %>%
    dplyr::select(col, sample_id) %>%
    dplyr::mutate(col = as.integer(col), sample_id = as.character(sample_id)) %>%
    dplyr::arrange(col)

  if (!nrow(col_map)) stop("parse_sample_map(): no (col -> sample_id) mapping could be derived from the Sample file.")
  col_map
}

raw_long <- parse_raw_blocks(raw_path)
col_map  <- parse_sample_map(sample_path)

raw_long <- raw_long %>%
  dplyr::left_join(col_map, by = "col")

# sanity checks
stopifnot("sample_id" %in% names(raw_long))
stopifnot(!any(is.na(raw_long$sample_id)))

# Helper: read the compound sheet & build mapping -------------------------
cmp <- read_tsv(cmp_path, show_col_types = FALSE)

# Map position (A1..A4, B1..B4) -> Compound name + Unit
key_map <- cmp %>%
  filter(ID %in% c("Compound", "Unit")) %>%
  pivot_longer(-ID, names_to = "pos", values_to = "val") %>%
  pivot_wider(names_from = ID, values_from = val) %>%
  rename(Unit = Unit)

# Concentrations per level (1..8) by position
dils <- cmp %>%
  filter(ID %in% as.character(1:8)) %>%
  mutate(level = as.integer(ID)) %>%
  dplyr::select(-ID) %>%
  pivot_longer(-level, names_to = "pos", values_to = "dose") %>%
  mutate(dose = as.numeric(dose)) %>%
  left_join(key_map, by = "pos") %>%
  mutate(plate = str_sub(pos, 1, 1),
         block = as.integer(str_sub(pos, 2, 2)))

# Tie raw columns to positions (triplicates per block) --------------------
dat <- raw_long %>%
  dplyr::mutate(rep = ((col - 1L) %% 3L) + 1L) %>%
  dplyr::left_join(dils, by = c("pos", "level", "plate", "block"))

# We only need fluorescence for viability normalisation/fits --------------
fluo <- dat %>%
  filter(meas == "Fluo") %>%
  dplyr::select(sample_id, plate, pos, block, rep, level, dose, Unit, Compound, value)

# Compute top (vehicle, level=1) and bottoms (min signal) -----------------
top_map <- fluo %>%
  filter(level == 1) %>%
  group_by(sample_id, pos) %>%
  summarise(top = mean(value, na.rm = TRUE), .groups = "drop")

min_by_pos <- fluo %>%
  dplyr::filter(level >= 2, is.finite(value)) %>%   # keep only finite values
  dplyr::group_by(sample_id, pos) %>%
  dplyr::summarise(min_signal = min(value), .groups = "drop")

# Bottoms per curve = within-curve minimum (levels 2–8)
bottom_map <- min_by_pos %>%
  dplyr::transmute(sample_id, pos, bottom = min_signal)

# Merge tops/bottoms and compute normalised viability in [0,1] ------------
fluo_n <- fluo %>%
  left_join(top_map,    by = c("sample_id", "pos")) %>%
  left_join(bottom_map, by = c("sample_id", "pos")) %>%
  mutate(viab = (value - bottom) / pmax(top - bottom, .Machine$double.eps),
         viab = pmin(pmax(viab, 0), 1))

# Fitting + prediction helpers --------------------------------------------
# Delta-method 95% CI for LL.4 with fixed lower=0, upper=1
predict_ll4_ci <- function(fit, new_dose) {
  cf <- coef(fit)                   # parameters are c(b, c, d, e), but c and d fixed below
  vc <- vcov(fit)
  # When lower/upper are fixed at 0/1, vcov should be 2x2 for (b, e)
  b  <- unname(cf["b:(Intercept)"])
  e  <- unname(cf["e:(Intercept)"])
  # Safety for name variants
  if (is.na(b)) b <- unname(cf[grep("^b", names(cf))[1]])
  if (is.na(e)) e <- unname(cf[grep("^e", names(cf))[1]])

  x  <- new_dose
  t  <- b * (log(x) - log(e))
  et <- exp(t)
  den <- (1 + et)
  f   <- 1 / den

  K   <- et / (den^2)                      # common factor
  ddb <- -K * (log(x) - log(e))            # df/db
  dde <-  K * ( b / e )                    # df/de

  # Build gradient and variance
  G <- cbind(ddb, dde)                     # N x 2
  vars <- rowSums((G %*% vc) * G)          # diag(G V G^T)
  se   <- sqrt(pmax(vars, 0))
  tibble(dose = x, fit = f,
         lwr = pmax(f - 1.96 * se, 0),
         upr = pmin(f + 1.96 * se, 1))
}

fit_one_curve <- function(df, fix_bounds = TRUE) {
  dfit <- df %>% dplyr::filter(level >= 2, is.finite(dose), dose > 0)
  if (nrow(dfit) < 4L) return(NULL)

  fct <- if (fix_bounds) LL.4(fixed = c(NA, 0, 1, NA)) else LL.4()
  fit <- try(drm(viab ~ dose, data = dfit, fct = fct, robust = "mean"), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)

  # IC50 (with SE & 95% CI)
  ic50 <- try(ED(fit, respLev = 0.5, type = "absolute", interval = "delta"), silent = TRUE)
  if (inherits(ic50, "try-error")) {
    ic50_df <- tibble(ED = NA_real_, SE = NA_real_, Lower = NA_real_, Upper = NA_real_)
  } else {
    ic50_df <- as_tibble(as.data.frame(ic50))
    # Standard ED() naming is (Estimate, Std. Error, Lower, Upper):
    names(ic50_df) <- c("ED", "SE", "Lower", "Upper")
  }

  # Prediction grid & CI ribbon
  rng  <- range(dfit$dose, na.rm = TRUE)
  grid <- tibble(dose = 10^seq(log10(rng[1]) - 0.05, log10(rng[2]) + 0.05, length.out = 200))
  pred <- predict_ll4_ci(fit, grid$dose)

  # Per-fit p-value for slope (dose–response trend)
  p_trend <- tryCatch({
    sm <- summary(fit)
    cm <- if (!is.null(sm$coefficients)) sm$coefficients else coef(sm)
    rn <- rownames(cm); colp <- grep("Pr", colnames(cm), value = TRUE)[1]
    idx <- grep("^b", rn)[1]
    as.numeric(cm[idx, colp])
  }, error = function(e) NA_real_)

  list(fit = fit,
       ic50 = ic50_df,
       pred = pred,
       points = dfit %>% dplyr::select(dose, viab, level, rep),
       p_trend = p_trend)
}

# Run fits for every sample × compound ------------------------------------
res_list <- fluo_n %>%
  dplyr::group_by(sample_id, Compound, Unit, pos) %>%
  dplyr::group_map(~{
    ans <- fit_one_curve(.x)
    tibble(
      sample_id = .y$sample_id, Compound = .y$Compound, Unit = .y$Unit, pos = .y$pos,
      fit_ok = !is.null(ans),
      p_trend = if (!is.null(ans)) ans$p_trend else NA_real_,
      ic50   = list(if (!is.null(ans)) ans$ic50 else tibble(ED=NA_real_, SE=NA_real_, Lower=NA_real_, Upper=NA_real_)),
      pred   = list(if (!is.null(ans)) ans$pred else tibble(dose=numeric(), fit=numeric(), lwr=numeric(), upr=numeric())),
      points = list(if (!is.null(ans)) ans$points else tibble(dose=numeric(), viab=numeric()))
    )
  }) %>% dplyr::bind_rows()
res_ann <- res_list

# IC50 table now includes SE and p_trend
ic50_tbl <- res_list %>%
  dplyr::transmute(sample_id, Compound, Unit, pos, fit_ok, p_trend, ic50) %>%
  tidyr::unnest(ic50) %>%
  dplyr::mutate(across(c(ED, SE, Lower, Upper), as.numeric)) %>%
  dplyr::arrange(sample_id, Compound)

readr::write_csv(ic50_tbl, file.path(out_dir, "csv", "IC50_results.csv"))

# --- Helpers for IC50 comparisons & pretty printing ---------------------------
ic50_compare <- function(e1, se1, e2, se2) {
  ok <- is.finite(e1) && is.finite(se1) && is.finite(e2) && is.finite(se2) && (e1 > 0) && (e2 > 0)
  if (!isTRUE(ok)) {
    return(tibble(ratio = NA_real_, ratio_lwr = NA_real_, ratio_upr = NA_real_, p_value = NA_real_))
  }
  l1 <- log(e1); l2 <- log(e2)
  s1l <- se1 / e1; s2l <- se2 / e2
  diff <- l1 - l2
  se   <- sqrt(s1l^2 + s2l^2)
  z    <- diff / se
  p    <- 2 * pnorm(abs(z), lower.tail = FALSE)
  tibble(ratio = exp(diff),
         ratio_lwr = exp(diff - 1.96*se),
         ratio_upr = exp(diff + 1.96*se),
         p_value = p)
}

pretty_p <- function(p) {
  if (!is.finite(p)) return("p = n/a")
  if (p < 1e-4) paste0("p = ", format(p, digits = 2, scientific = TRUE))
  paste0("p = ", signif(p, 3))
}

fmt_num <- function(x) {
  if (!is.finite(x)) return("n/a")
  if (x >= 0.01) formatC(x, format = "f", digits = 3)
  format(x, digits = 2, scientific = TRUE)
}

# (A) sample vs samplecse for each compound — robust self-join (no list-cols, no char numerics)
ic50_tbl2 <- ic50_tbl %>%
  dplyr::mutate(
    sample_id = as.character(sample_id),
    base      = tolower(sub("cse$", "", sample_id, ignore.case = TRUE)),
    is_cse    = grepl("cse$", sample_id, ignore.case = TRUE)
  )

non <- ic50_tbl2 %>%
  dplyr::filter(!is_cse) %>%
  dplyr::select(base, Compound, Unit, sample_id, ED, SE) %>%
  dplyr::rename(
    IC50_sample = ED,
    SE_sample   = SE
  )

cse <- ic50_tbl2 %>%
  dplyr::filter(is_cse) %>%
  dplyr::select(base, Compound, Unit, sample_id, ED, SE) %>%
  dplyr::rename(
    samplecse_id    = sample_id,
    IC50_sampleCSE  = ED,
    SE_sampleCSE    = SE
  )

cmp_across <- non %>%
  dplyr::inner_join(cse, by = c("base","Compound","Unit")) %>%
  dplyr::rename(sample_id = sample_id, unit = Unit) %>%   # <- add this rename
  dplyr::mutate(dplyr::across(
    c(IC50_sample, SE_sample, IC50_sampleCSE, SE_sampleCSE),
    ~ suppressWarnings(as.numeric(.))
  )) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(ic = ic50_compare(IC50_sampleCSE, SE_sampleCSE, IC50_sample, SE_sample)) %>%
  tidyr::unnest_wider(ic) %>%
  dplyr::arrange(base, Compound)

readr::write_csv(cmp_across, file.path(out_dir, "csv", "IC50_compare_sample_vs_sampleCSE.csv"))

# (B) Within-sample: generic <compound> vs <compound>_CSE
cmp_pairs <- ic50_tbl %>%
  dplyr::group_by(sample_id, comp_base = sub("_CSE$", "", Compound)) %>%
  dplyr::filter(any(grepl("_CSE$", Compound)) & any(!grepl("_CSE$", Compound))) %>%
  dplyr::summarise(
    unit      = dplyr::coalesce(first(Unit[!grepl("_CSE$", Compound)]), first(Unit[grepl("_CSE$", Compound)])),
    e_nonCSE  = ED[!grepl("_CSE$", Compound)][1],
    se_nonCSE = SE[!grepl("_CSE$", Compound)][1],
    e_CSE     = ED[ grepl("_CSE$", Compound)][1],
    se_CSE    = SE[ grepl("_CSE$", Compound)][1],
    .groups = "drop"
  ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(ic = ic50_compare(e_CSE, se_CSE, e_nonCSE, se_nonCSE)) %>%  # ratio = CSE / non-CSE
  tidyr::unnest(ic) %>%
  dplyr::rename(IC50_nonCSE = e_nonCSE, SE_nonCSE = se_nonCSE,
                IC50_CSE    = e_CSE,    SE_CSE    = se_CSE) %>%
  dplyr::arrange(sample_id, comp_base)

readr::write_csv(cmp_pairs, file.path(out_dir, "csv", "IC50_compare_generic_pairs.csv"))

# --- By-sample bundles: base compound + all co-treatments on ONE plot ---------
dir.create(file.path(out_dir, "plots_by_sample_bundles"), recursive = TRUE, showWarnings = FALSE)

# Extract a "base" stem (text before the first underscore); e.g. Rotenone, Rotenone_CSE, Rotenone_NAC -> "Rotenone"
comp_base_of <- function(x) sub("_.*$", "", x)

# Okabe–Ito colourblind-friendly palette (we reserve blue for the base)
cb_others <- c("#E69F00", "#009E73", "#CC79A7", "#D55E00", "#F0E442", "#000000", "#56B4E9")
base_blue <- "#0072B2"

res_bundles <- res_ann %>%
  dplyr::filter(fit_ok) %>%
  dplyr::mutate(comp_base = comp_base_of(Compound),
                is_base   = !grepl("_", Compound)) %>%
  dplyr::group_by(sample_id, comp_base) %>%
  # keep only (sample_id, base) pairs where the pure base exists
  dplyr::filter(any(is_base)) %>%
  dplyr::ungroup()

# One plot per (sample_id × comp_base)
pairs_to_plot <- res_bundles %>%
  dplyr::distinct(sample_id, comp_base)

for (i in seq_len(nrow(pairs_to_plot))) {
  sid <- pairs_to_plot$sample_id[i]
  cb  <- pairs_to_plot$comp_base[i]

  sub <- res_bundles %>%
    dplyr::filter(sample_id == sid, comp_base == cb) %>%
    # one entry per compound
    dplyr::group_by(Compound) %>% dplyr::slice(1) %>% dplyr::ungroup()

  if (!nrow(sub)) next

  # Build data for ribbons/lines/points and IC50 vlines
  preds <- purrr::map_dfr(seq_len(nrow(sub)), function(j) {
    sub$pred[[j]] %>% dplyr::mutate(group = sub$Compound[[j]])
  })
  pts <- purrr::map_dfr(seq_len(nrow(sub)), function(j) {
    sub$points[[j]] %>% dplyr::mutate(group = sub$Compound[[j]])
  })
  vdat <- tibble::tibble(
    dose  = purrr::map_dbl(sub$ic50, ~ suppressWarnings(as.numeric(.x[1, "ED"]))),
    group = sub$Compound
  ) %>% dplyr::filter(is.finite(dose))

  # Colours: base = blue; all others = CB palette (recycled if needed)
  groups <- unique(c(preds$group, pts$group))
  cols   <- setNames(rep_len(cb_others, length(groups)), groups)
  if (cb %in% names(cols)) cols[cb] <- base_blue  # pure base present
  # If the "pure base" appears in the data under its full name (no suffix), ensure it’s blue:
  base_idx <- which(names(cols) == cb)
  if (length(base_idx)) cols[base_idx] <- base_blue

  unit_lbl <- sub$Unit[which(!is.na(sub$Unit))[1]]
  ttl <- glue::glue("{sid} — {cb} (base) ± co-treatments")

  p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = preds, ggplot2::aes(dose, ymin = lwr, ymax = upr, fill = group), alpha = 0.18) +
    ggplot2::geom_line(  data = preds, ggplot2::aes(dose, fit, colour = group), linewidth = 0.9) +
    ggplot2::geom_point( data = pts,   ggplot2::aes(dose, viab, colour = group), size = 2, alpha = 0.85) +
    ggplot2::geom_vline(data = vdat, ggplot2::aes(xintercept = dose, colour = group),
                        linetype = 2, linewidth = 0.7, alpha = 0.9) +
    ggplot2::scale_x_log10(labels = scales::label_number(),
                           expand = ggplot2::expansion(mult = c(0.02, 0.05))) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::scale_colour_manual(values = cols, name = "Treatment") +
    ggplot2::scale_fill_manual(values = scales::alpha(cols, 0.18), guide = "none") +
    ggplot2::labs(title = ttl,
                  x = glue::glue("Concentration ({unit_lbl})"),
                  y = "Normalised viability") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                   panel.grid.minor.x = ggplot2::element_blank())

  fn <- file.path(out_dir, "plots_by_sample_bundles",
                  glue::glue("bundle_{gsub('[^A-Za-z0-9]+','_', sid)}_{gsub('[^A-Za-z0-9]+','_', cb)}.png"))
  ggplot2::ggsave(fn, p, width = 7.2, height = 5.0, dpi = 300)
}

# --- Publication tables: bootstrap medians with 95% CIs -----------------------
set.seed(1234)

boot_ci <- function(x, B = 5000L) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(c(NA_real_, NA_real_))
  bs <- replicate(B, stats::median(sample(x, replace = TRUE), na.rm = TRUE))
  stats::quantile(bs, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
}

# (1) IC50 by compound (pooled across samples)
pub_ic50 <- ic50_tbl %>%
  dplyr::group_by(Compound, Unit) %>%
  dplyr::summarise(
    n = sum(is.finite(ED)),
    median = stats::median(ED, na.rm = TRUE),
    ci = list(boot_ci(ED)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    CI_low = purrr::map_dbl(ci, 1),
    CI_high = purrr::map_dbl(ci, 2)
  ) %>%
  dplyr::select(-ci) %>%
  dplyr::arrange(Compound)

readr::write_csv(pub_ic50, file.path(out_dir, "csv", "pub_IC50_by_compound_bootstrap.csv"))

# (2) Across-sample ratios (sampleCSE / sample) per compound
cmp_across <- readr::read_csv(file.path(out_dir, "csv", "IC50_compare_sample_vs_sampleCSE.csv"), show_col_types = FALSE)

pub_ratio_across <- cmp_across %>%
  dplyr::group_by(Compound, unit) %>%
  dplyr::summarise(
    n = sum(is.finite(ratio)),
    median_ratio = stats::median(ratio, na.rm = TRUE),
    ci = list(boot_ci(ratio)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    CI_low = purrr::map_dbl(ci, 1),
    CI_high = purrr::map_dbl(ci, 2)
  ) %>%
  dplyr::select(-ci) %>%
  dplyr::arrange(Compound)

readr::write_csv(pub_ratio_across, file.path(out_dir, "csv", "pub_Ratio_sampleCSE_over_sample_bootstrap.csv"))

# (3) Within-sample ratios for generic pairs (CSE / non-CSE), pooled across samples
pub_ratio_pairs <- cmp_pairs %>%
  dplyr::group_by(comp_base, unit) %>%
  dplyr::summarise(
    n = sum(is.finite(ratio)),
    median_ratio = stats::median(ratio, na.rm = TRUE),
    ci = list(boot_ci(ratio)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    CI_low  = purrr::map_dbl(ci, 1),
    CI_high = purrr::map_dbl(ci, 2)
  ) %>%
  dplyr::select(-ci) %>%
  dplyr::arrange(comp_base)

readr::write_csv(pub_ratio_pairs, file.path(out_dir, "csv", "pub_Ratio_pairs_generic_bootstrap.csv"))

message("Done. IC50s: ", file.path(out_dir, "csv", "IC50_results.csv"),
        " — Plots in: ", file.path(out_dir, "plots"))