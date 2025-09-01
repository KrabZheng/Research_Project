# Resazurin Exp36 — per-sample dose–response with drc IC50s + CI ribbons
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
raw_path <- "Resazurin_IC50/Resazurin_Exp_36_Raw.txt"
cmp_path <- "Resazurin_IC50/Resazurin_Exp_36_Compound.txt"
sample_path <- "Resazurin_IC50/Resazurin_Exp_36_Sample.txt"
out_dir  <- "Resazurin_IC50/resazurin_outputs/Exp36"

dir.create(out_dir, showWarnings = FALSE)
dir.create(file.path(out_dir, "plots"), showWarnings = FALSE)
dir.create(file.path(out_dir, "csv"),   showWarnings = FALSE)

# Helper: parse RAW (Exp34 format: A1/B1 headers; 8 lines; 12 Fluo + 12 Abs + 12 BG per line)
parse_raw_blocks <- function(path) {
  lines  <- readLines(path, warn = FALSE)
  hdr_ix <- grep("^([AB][1-4])_Fluo\\t1\\t2\\t3\\t4\\t5\\t6\\t7\\t8\\t9\\t10\\t11\\t12", lines)
  stopifnot(length(hdr_ix) > 0)

  purrr::map_dfr(hdr_ix, function(i) {
    hdr   <- lines[[i]]
    pos   <- stringr::str_match(hdr, "^([AB][1-4])_Fluo")[,2]
    plate <- stringr::str_sub(pos, 1, 1)
    block <- as.integer(stringr::str_sub(pos, 2, 2))

    block_lines <- lines[(i+1):(i+8)]  # 8 dose levels
    purrr::imap_dfr(block_lines, function(ln, idx) {
      parts   <- stringr::str_split(ln, "\\t")[[1]]   # KEEP empties
      lvl_raw <- suppressWarnings(as.integer(parts[1]))
      lvl     <- if (!is.na(lvl_raw) && lvl_raw %in% 1:8) lvl_raw else idx

      # Fluo is always the first 12 numeric cells after the first column (which may be empty)
      if (length(parts) < 13) stop("Raw data line has fewer than 13 tab fields:\n", ln)
      fl <- suppressWarnings(as.numeric(parts[2:13]))

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

raw_long %>% dplyr::count(pos, level) %>% print(n = 50)   # should be n == 12 for every pos×level

# sanity checks
stopifnot("sample_id" %in% names(raw_long))
stopifnot(!any(is.na(raw_long$sample_id)))

# Helper: read the compound sheet & build mapping -------------------------
cmp <- readr::read_tsv(cmp_path, col_types = readr::cols(.default = readr::col_character())) %>%
  dplyr::mutate(ID = as.character(ID))

# Map position (A1..A4, B1..B4) -> Compound name + Unit
key_map <- cmp %>%
  dplyr::filter(ID %in% c("Compound", "Unit")) %>%
  tidyr::pivot_longer(-ID, names_to = "pos", values_to = "val") %>%
  tidyr::pivot_wider(names_from = ID, values_from = val) %>%
  dplyr::rename(Unit = Unit)

# Concentrations per level (1..8) by position
dils <- cmp %>%
  dplyr::filter(ID %in% as.character(1:8)) %>%
  dplyr::mutate(level = as.integer(ID)) %>%
  dplyr::select(-ID) %>%
  tidyr::pivot_longer(-level, names_to = "pos", values_to = "dose") %>%
  dplyr::mutate(dose = suppressWarnings(as.numeric(dose))) %>%
  dplyr::left_join(key_map, by = "pos") %>%
  dplyr::mutate(
    plate = stringr::str_sub(pos, 1, 1),
    block = as.integer(stringr::str_sub(pos, 2, 2))
  )

# Sanity checks (fail fast if mapping is broken)
stopifnot(nrow(key_map) > 0, nrow(dils) > 0)

# Tie raw columns to positions (triplicates per block)
dat <- raw_long %>%
  dplyr::mutate(rep = ((col - 1L) %% 3L) + 1L) %>%
  dplyr::left_join(dils, by = c("pos", "level", "plate", "block"))

# quick sanity on the raw layout (should be n==12 for each pos×level)
raw_long %>% dplyr::count(pos, level) %>% print(n = 50)

# --- Compound sheet -> dose mapping -------------------------------------
cmp <- readr::read_tsv(cmp_path, col_types = readr::cols(.default = readr::col_character())) %>%
  dplyr::mutate(ID = as.character(ID))

key_map <- cmp %>%
  dplyr::filter(ID %in% c("Compound", "Unit")) %>%
  tidyr::pivot_longer(-ID, names_to = "pos", values_to = "val") %>%
  tidyr::pivot_wider(names_from = ID, values_from = val) %>%
  dplyr::rename(Unit = Unit)

dils <- cmp %>%
  dplyr::filter(ID %in% as.character(1:8)) %>%
  dplyr::mutate(level = as.integer(ID)) %>%
  dplyr::select(-ID) %>%
  tidyr::pivot_longer(-level, names_to = "pos", values_to = "dose") %>%
  dplyr::mutate(dose = suppressWarnings(as.numeric(dose))) %>%
  dplyr::left_join(key_map, by = "pos") %>%
  dplyr::mutate(
    plate = stringr::str_sub(pos, 1, 1),
    block = as.integer(stringr::str_sub(pos, 2, 2))
  )

stopifnot(nrow(key_map) > 0, nrow(dils) > 0)

# --- Build dat (now it's safe to reference it) ---------------------------
dat <- raw_long %>%
  dplyr::mutate(rep = ((col - 1L) %% 3L) + 1L) %>%
  dplyr::left_join(dils, by = c("pos","level","plate","block"))

# doses must be present for all fluorescence rows
if (any(!is.finite(dat$dose[dat$meas == "Fluo"]))) {
  bad <- dat %>% dplyr::filter(meas == "Fluo", !is.finite(dose)) %>% dplyr::count(pos, level)
  print(bad, n = 50)
  stop("Dose mapping failed for some (pos, level). Check Compound.txt positions/levels.")
}

# Should be 32 pos×level combos, all with n=12
dat %>% dplyr::filter(meas == "Fluo") %>% dplyr::count(pos, level) %>% print(n = Inf)

# --- Now derive fluo and (optionally) inspect a sample -------------------
fluo <- dat %>%
  dplyr::filter(meas == "Fluo") %>%
  dplyr::select(sample_id, plate, pos, block, rep, level, dose, Unit, Compound, value)

# Optional spot-checks AFTER fluo exists
fluo %>% dplyr::filter(sample_id == "ASC399cse") %>% dplyr::count(pos, level)
fluo %>% dplyr::filter(sample_id == "ASC399cse") %>%
  dplyr::summarise(min_val = min(value, na.rm = TRUE),
                   max_val = max(value, na.rm = TRUE))

# We only need fluorescence for viability normalisation/fits --------------
fluo <- dat %>%
  filter(meas == "Fluo") %>%
  dplyr::select(sample_id, plate, pos, block, rep, level, dose, Unit, Compound, value)

# Make sure we have full 12x8 for that sample
fluo %>% filter(sample_id == "ASC399cse") %>% count(pos, level)

# Look at the raw value range (shouldn’t be all ~0)
fluo %>% filter(sample_id == "ASC399cse") %>%
  summarise(min_val = min(value, na.rm = TRUE),
            max_val = max(value, na.rm = TRUE))

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
  left_join(top_map, by = c("sample_id","pos")) %>%
  left_join(bottom_map, by = c("sample_id","pos")) %>%
  mutate(viab = (value - bottom) / pmax(top - bottom, .Machine$double.eps),
         viab = pmin(pmax(viab, 1e-6), 1 - 1e-6))

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

fit_one_curve <- function(df, fix_bounds = TRUE, B_boot = 600L) {
  dfit <- df %>%
    dplyr::filter(level >= 2, is.finite(dose), dose > 0, is.finite(viab))
  if (nrow(dfit) < 4L || dplyr::n_distinct(dfit$dose) < 3L) return(NULL)

  fct <- if (fix_bounds) LL.4(fixed = c(NA, 0, 1, NA)) else LL.4()
  fit <- try(drm(viab ~ dose, data = dfit, fct = fct, robust = "mean"), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)

  # --- IC50 = 'e' ---------------------------------------------------------
  cf <- coef(fit)
  en <- grep("^e", names(cf), value = TRUE)[1]
  e_val <- as.numeric(cf[en])
  se_val <- NA_real_
  lower  <- NA_real_
  upper  <- NA_real_

  # 1) SE from vcov
  vc <- try(vcov(fit), silent = TRUE)
  if (!inherits(vc, "try-error") && !is.null(dim(vc)) && is.finite(e_val) &&
      !is.null(en) && en %in% rownames(vc)) {
    vv <- vc[en, en]
    if (is.finite(vv) && vv >= 0) se_val <- sqrt(vv)
  }

  # 2) SE from summary() if needed
  if (!is.finite(se_val)) {
    sm <- try(summary(fit), silent = TRUE)
    if (!inherits(sm, "try-error") && !is.null(sm$coefficients)) {
      rn <- rownames(sm$coefficients)
      idx <- grep("^e", rn)[1]
      if (!is.na(idx)) {
        se_candidate <- as.numeric(sm$coefficients[idx, "Std. Error"])
        if (is.finite(se_candidate)) se_val <- se_candidate
      }
    }
  }

  # 3) If still no SE, use fiducial limits (ED(..., "fls"))
  if (!is.finite(se_val)) {
    ed_fls <- try(ED(fit, respLev = 0.5, type = "absolute", interval = "fls"), silent = TRUE)
    if (!inherits(ed_fls, "try-error")) {
      df2 <- as.data.frame(ed_fls)
      cn  <- tolower(colnames(df2))
      grab <- function(name) {
        j <- which(cn == tolower(name) | grepl(paste0("^", tolower(name)), cn))
        if (length(j)) as.numeric(df2[1, j[1]]) else NA_real_
      }
      est2   <- grab("estimate")
      lower2 <- grab("lower")
      upper2 <- grab("upper")
      if (!is.finite(e_val) && is.finite(est2)) e_val <- est2
      if (isTRUE(all(is.finite(c(lower2, upper2))))) {
        lower  <- lower2
        upper  <- upper2
        se_val <- (upper2 - lower2) / (2 * 1.96)
      }
    }
  }

  # 4) Last resort: bootstrap ED (parameter 'e')
  if (!is.finite(se_val)) {
    set.seed(1234)
    n <- nrow(dfit)
    boot_vals <- replicate(B_boot, {
      idx <- sample.int(n, n, replace = TRUE)
      dd <- dfit[idx, , drop = FALSE]
      ft <- try(drm(viab ~ dose, data = dd, fct = fct, robust = "mean"), silent = TRUE)
      if (inherits(ft, "try-error")) return(NA_real_)
      cfe <- try(coef(ft), silent = TRUE)
      if (inherits(cfe, "try-error")) return(NA_real_)
      en2 <- grep("^e", names(cfe), value = TRUE)[1]
      as.numeric(cfe[en2])
    })
    boot_vals <- boot_vals[is.finite(boot_vals)]
    if (length(boot_vals) >= 20L) {
      if (!is.finite(e_val)) e_val <- stats::median(boot_vals)
      se_val <- stats::sd(boot_vals)
      lower  <- stats::quantile(boot_vals, 0.025, names = FALSE)
      upper  <- stats::quantile(boot_vals, 0.975, names = FALSE)
    }
  }

  ic50_df <- tibble::tibble(
    ED    = as.numeric(e_val),
    SE    = as.numeric(se_val),
    Lower = as.numeric(if (is.finite(lower)) lower else if (is.finite(se_val)) e_val - 1.96 * se_val else NA_real_),
    Upper = as.numeric(if (is.finite(upper)) upper else if (is.finite(se_val)) e_val + 1.96 * se_val else NA_real_)
  )

  # Predictions (with CI if vcov is well-behaved)
  rng  <- range(dfit$dose, na.rm = TRUE)
  grid <- tibble(dose = 10^seq(log10(rng[1]) - 0.05, log10(rng[2]) + 0.05, length.out = 200))
  pred <- try(predict_ll4_ci(fit, grid$dose), silent = TRUE)
  if (inherits(pred, "try-error")) {
    bpar <- as.numeric(cf[grep("^b", names(cf))[1]])
    epar <- as.numeric(cf[grep("^e", names(cf))[1]])
    f    <- function(x) 1 / (1 + exp(bpar * (log(x) - log(epar))))
    pred <- tibble(dose = grid$dose, fit = f(grid$dose), lwr = NA_real_, upr = NA_real_)
  }

  p_trend <- tryCatch({
    sm <- summary(fit)
    cm <- if (!is.null(sm$coefficients)) sm$coefficients else coef(sm)
    rn <- rownames(cm); colp <- grep("Pr", colnames(cm), value = TRUE)[1]
    idx <- grep("^b", rn)[1]
    as.numeric(cm[idx, colp])
  }, error = function(e) NA_real_)

  list(
    fit    = fit,
    ic50   = ic50_df,
    pred   = pred,
    points = dfit %>% dplyr::select(dose, viab, level, rep),
    p_trend = p_trend
  )
}

# --- Helper: IC50 ratio comparison (delta method on log-IC50) ---------------
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
  tibble(
    ratio     = exp(diff),
    ratio_lwr = exp(diff - 1.96 * se),
    ratio_upr = exp(diff + 1.96 * se),
    p_value   = p
  )
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

# IC50 table now includes SE and p_trend
ic50_tbl <- res_list %>%
  dplyr::transmute(sample_id, Compound, Unit, pos, fit_ok, p_trend, ic50) %>%
  tidyr::unnest(ic50) %>%
  dplyr::mutate(across(c(ED, SE, Lower, Upper), as.numeric)) %>%
  dplyr::arrange(sample_id, Compound)

readr::write_csv(ic50_tbl, file.path(out_dir, "csv", "IC50_results.csv"))

# ==== Baseline vs co-treatment (within-sample) ============================
# Expect compounds like: KBrO3 (baseline) and KBrO3_SFN / KBrO3_NAC / KBrO3_Erastin

ic50_base <- ic50_tbl %>%
  dplyr::mutate(
    comp_base = sub("_.+$", "", Compound),                # "KBrO3" from "KBrO3_SFN"
    comp_var  = dplyr::if_else(grepl("_", Compound), sub("^[^_]+_", "", Compound), "BASELINE"),
    is_cse    = grepl("cse$", sample_id, ignore.case = TRUE)
  )

# Split baseline and variants
base_only <- ic50_base %>%
  dplyr::filter(comp_var == "BASELINE") %>%
  dplyr::select(sample_id, is_cse, comp_base, Unit, IC50_base = ED, SE_base = SE)

vars_only <- ic50_base %>%
  dplyr::filter(comp_var != "BASELINE") %>%
  dplyr::select(sample_id, is_cse, comp_base, Unit, comp_var, IC50_var = ED, SE_var = SE)

# Pair each variant with its baseline in the SAME sample_id (ASC397 and ASC397cse are separate rows)
cotreat_pairs <- vars_only %>%
  dplyr::inner_join(base_only, by = c("sample_id","is_cse","comp_base","Unit")) %>%
  dplyr::mutate(dplyr::across(c(IC50_var, SE_var, IC50_base, SE_base), ~ suppressWarnings(as.numeric(.)))) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(ic = ic50_compare(IC50_var, SE_var, IC50_base, SE_base)) %>%  # ratio = variant / baseline
  tidyr::unnest_wider(ic) %>%
  dplyr::arrange(sample_id, is_cse, comp_base, comp_var)

readr::write_csv(cotreat_pairs, file.path(out_dir, "csv", "IC50_compare_baseline_vs_cotreatments.csv"))

# Publication-style summary: median ratio (variant / baseline) with bootstrap CIs
set.seed(1234)
boot_ci <- function(x, B = 5000L) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(c(NA_real_, NA_real_))
  bs <- replicate(B, stats::median(sample(x, replace = TRUE), na.rm = TRUE))
  stats::quantile(bs, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
}

pub_cotreat <- cotreat_pairs %>%
  dplyr::group_by(comp_base, comp_var, is_cse, Unit) %>%
  dplyr::summarise(
    n = sum(is.finite(ratio)),
    median_ratio = stats::median(ratio, na.rm = TRUE),
    ci = list(boot_ci(ratio)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(CI_low = purrr::map_dbl(ci, 1),
                CI_high = purrr::map_dbl(ci, 2)) %>%
  dplyr::select(-ci) %>%
  dplyr::arrange(comp_base, comp_var, is_cse)

readr::write_csv(pub_cotreat, file.path(out_dir, "csv", "pub_Ratio_baseline_vs_cotreatments_bootstrap.csv"))

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

# (A) sample vs samplecse for each compound — robust self-join (dedup + keep Unit)
ic50_tbl2 <- ic50_tbl %>%
  dplyr::mutate(
    sample_id = as.character(sample_id),
    base      = tolower(sub("cse$", "", sample_id, ignore.case = TRUE)),
    is_cse    = grepl("cse$", sample_id, ignore.case = TRUE)
  )

non <- ic50_tbl2 %>%
  dplyr::filter(!is_cse) %>%
  dplyr::distinct(base, Compound, Unit, .keep_all = TRUE) %>%   # <-- dedupe
  dplyr::select(base, Compound, Unit, sample_id, ED, SE) %>%
  dplyr::rename(IC50_sample = ED, SE_sample = SE)

cse <- ic50_tbl2 %>%
  dplyr::filter(is_cse) %>%
  dplyr::distinct(base, Compound, Unit, .keep_all = TRUE) %>%   # <-- dedupe
  dplyr::select(base, Compound, Unit, sample_id, ED, SE) %>%
  dplyr::rename(samplecse_id = sample_id, IC50_sampleCSE = ED, SE_sampleCSE = SE)

cmp_across <- non %>%
  dplyr::inner_join(cse, by = c("base","Compound","Unit")) %>%
  dplyr::mutate(
    unit = Unit,  # keep Unit and also provide 'unit' for downstream
    dplyr::across(c(IC50_sample, SE_sample, IC50_sampleCSE, SE_sampleCSE),
                  ~ suppressWarnings(as.numeric(.)))
  ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(ic = ic50_compare(IC50_sampleCSE, SE_sampleCSE, IC50_sample, SE_sample)) %>%
  tidyr::unnest_wider(ic) %>%
  dplyr::arrange(base, Compound)

readr::write_csv(cmp_across, file.path(out_dir, "csv", "IC50_compare_sample_vs_sampleCSE.csv"))

# --- Combined plots -----------------------------------------------------------
# Folders for outputs
dir.create(file.path(out_dir, "plots_combined", "across_sample"), recursive = TRUE, showWarnings = FALSE)

# Helper to derive base id (pair sample with its "...cse")
base_id <- function(s) tolower(gsub("cse$", "", s, ignore.case = TRUE))

# Annotate res_list with pairing helpers
res_ann <- res_list %>%
  dplyr::mutate(
    base  = base_id(sample_id),
    is_cse = grepl("cse$", sample_id, ignore.case = TRUE)
  )

valid_bases <- res_ann %>%
  dplyr::distinct(base, is_cse) %>%
  dplyr::count(base, name = "n_groups") %>%
  dplyr::filter(n_groups == 2) %>%
  dplyr::pull(base)

# ---------- (1) Across sample vs samplecse: each compound, two curves ----------
# Only keep bases where both sample and samplecse exist
for (b in valid_bases) {
  sub <- res_ann %>% dplyr::filter(base == b, fit_ok)

  comps <- intersect(
    unique(sub %>% dplyr::filter(!is_cse) %>% dplyr::pull(Compound)),
    unique(sub %>% dplyr::filter( is_cse) %>% dplyr::pull(Compound))
  )

  for (cmp_name in comps) {
    d0 <- sub %>% dplyr::filter(Compound == cmp_name, !is_cse) %>% dplyr::slice(1)
    d1 <- sub %>% dplyr::filter(Compound == cmp_name,  is_cse) %>% dplyr::slice(1)
    if (nrow(d0) == 0 || nrow(d1) == 0) next

    pred <- dplyr::bind_rows(
      d0$pred[[1]]   %>% dplyr::mutate(group = paste0(d0$sample_id[[1]])),
      d1$pred[[1]]   %>% dplyr::mutate(group = paste0(d1$sample_id[[1]]))
    )
    pts  <- dplyr::bind_rows(
      d0$points[[1]] %>% dplyr::mutate(group = paste0(d0$sample_id[[1]])),
      d1$points[[1]] %>% dplyr::mutate(group = paste0(d1$sample_id[[1]]))
    )

    # IC50 vertical lines only
    ic0 <- d0$ic50[[1]][1, "ED", drop = TRUE]
    ic1 <- d1$ic50[[1]][1, "ED", drop = TRUE]
    vdat <- tibble::tibble(dose = c(ic0, ic1),
                           group = c(d0$sample_id[[1]], d1$sample_id[[1]]))

    unit_lbl  <- dplyr::coalesce(d0$Unit[[1]], d1$Unit[[1]])
    title_lbl <- glue::glue("{toupper(b)} — {cmp_name} (sample vs samplecse)")

    p <- ggplot2::ggplot() +
      ggplot2::geom_ribbon(data = pred, ggplot2::aes(dose, ymin = lwr, ymax = upr, fill = group), alpha = 0.18) +
      ggplot2::geom_line(  data = pred, ggplot2::aes(dose, fit, colour = group), linewidth = 0.9) +
      ggplot2::geom_point( data = pts,  ggplot2::aes(dose, viab, colour = group), size = 2, alpha = 0.85) +
      ggplot2::geom_vline(data = vdat, ggplot2::aes(xintercept = dose, colour = group),
                          linetype = 2, linewidth = 0.7, alpha = 0.9) +
      ggplot2::scale_x_log10(labels = scales::label_number(),
                              expand = ggplot2::expansion(mult = c(0.02, 0.05))) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::labs(title = title_lbl,
                    x = glue::glue("Concentration ({unit_lbl})"),
                    y = "Normalised viability", colour = NULL, fill = NULL) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                     panel.grid.minor.x = ggplot2::element_blank())

    fn <- file.path(out_dir, "plots_combined", "across_sample",
                    glue::glue("combined_{b}_{gsub('[^A-Za-z0-9]+','_', cmp_name)}.png"))
    ggplot2::ggsave(fn, p, width = 6.8, height = 4.8, dpi = 300)
  }
}

# ==== Plots: Baseline vs co-treatment curves per sample ===================
dir.create(file.path(out_dir, "plots_combined", "cotreatment_pairs"), recursive = TRUE, showWarnings = FALSE)

# Helper: pull one fitted entry from res_list by sample_id + exact Compound name
get_fit_entry <- function(sid, cmp_name) {
  res_list %>% dplyr::filter(sample_id == sid, Compound == cmp_name, fit_ok) %>% dplyr::slice(1)
}

# Find all (sample_id × comp_base) that have baseline + >=1 variant
have_pairs <- ic50_base %>%
  dplyr::distinct(sample_id, comp_base, Compound) %>%
  dplyr::mutate(is_variant = grepl("_", Compound)) %>%
  dplyr::group_by(sample_id, comp_base) %>%
  dplyr::summarise(
    has_baseline = any(!is_variant),
    variants = list(Compound[is_variant]),
    .groups = "drop"
  ) %>%
  dplyr::filter(has_baseline, lengths(variants) > 0)

for (i in seq_len(nrow(have_pairs))) {
  sid <- have_pairs$sample_id[i]
  base_name <- have_pairs$comp_base[i]
  var_list <- have_pairs$variants[[i]]

  d_base <- get_fit_entry(sid, base_name)
  if (!nrow(d_base)) next
  unit_lbl <- d_base$Unit[[1]]

  for (vfull in var_list) {
    d_var <- get_fit_entry(sid, vfull)
    if (!nrow(d_var)) next

    pred <- dplyr::bind_rows(
      d_base$pred[[1]] %>% dplyr::mutate(group = base_name),
      d_var$pred[[1]]  %>% dplyr::mutate(group = vfull)
    )
    pts <- dplyr::bind_rows(
      d_base$points[[1]] %>% dplyr::mutate(group = base_name),
      d_var$points[[1]]  %>% dplyr::mutate(group = vfull)
    )
    vdat <- tibble::tibble(
      dose  = c(d_base$ic50[[1]][1, "ED", drop=TRUE], d_var$ic50[[1]][1, "ED", drop=TRUE]),
      group = c(base_name, vfull)
    )

    p <- ggplot2::ggplot() +
      ggplot2::geom_ribbon(data = pred, ggplot2::aes(dose, ymin = lwr, ymax = upr, fill = group), alpha = 0.18) +
      ggplot2::geom_line(  data = pred, ggplot2::aes(dose, fit, colour = group), linewidth = 0.9) +
      ggplot2::geom_point( data = pts,  ggplot2::aes(dose, viab, colour = group), size = 2, alpha = 0.85) +
      ggplot2::geom_vline(data = vdat, ggplot2::aes(xintercept = dose, colour = group),
                          linetype = 2, linewidth = 0.7, alpha = 0.9) +
      ggplot2::scale_x_log10(labels = scales::label_number(),
                              expand = ggplot2::expansion(mult = c(0.02, 0.05))) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::labs(title = glue::glue("{sid} — {base_name} vs {vfull}"),
                    x = glue::glue("Concentration ({unit_lbl})"),
                    y = "Normalised viability", colour = NULL, fill = NULL) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                     panel.grid.minor.x = ggplot2::element_blank())

    fn <- file.path(out_dir, "plots_combined", "cotreatment_pairs",
                    glue::glue("pair_{gsub('[^A-Za-z0-9]+','_', sid)}_{gsub('[^A-Za-z0-9]+','_', base_name)}_vs_{gsub('[^A-Za-z0-9]+','_', vfull)}.png"))
    ggplot2::ggsave(fn, p, width = 6.8, height = 4.8, dpi = 300)
  }
}

# --- Publication tables: bootstrap medians with 95% CIs (Exp35) -------------
set.seed(1234)

boot_ci <- function(x, B = 5000L) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(c(NA_real_, NA_real_))
  bs <- replicate(B, stats::median(sample(x, replace = TRUE), na.rm = TRUE))
  stats::quantile(bs, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
}

# (1) IC50 by compound (pooled across samples) ---------------------------------
pub_ic50 <- ic50_tbl %>%
  dplyr::group_by(Compound, Unit) %>%
  dplyr::summarise(
    n = sum(is.finite(ED)),
    median = stats::median(ED, na.rm = TRUE),
    ci = list(boot_ci(ED)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    CI_low  = purrr::map_dbl(ci, 1),
    CI_high = purrr::map_dbl(ci, 2)
  ) %>%
  dplyr::select(-ci) %>%
  dplyr::arrange(Compound)

readr::write_csv(pub_ic50, file.path(out_dir, "csv", "pub_IC50_by_compound_bootstrap.csv"))

# (2) Across-sample ratios (sampleCSE / sample) per compound -------------------
# (works off cmp_across from earlier; normalize unit name)
cmp_across_tbl <- cmp_across   # already contains both Unit and unit

pub_ratio_across <- cmp_across_tbl %>%
  dplyr::group_by(Compound, unit) %>%
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
  dplyr::arrange(Compound)

readr::write_csv(pub_ratio_across, file.path(out_dir, "csv", "pub_Ratio_sampleCSE_over_sample_bootstrap.csv"))

# (3) Baseline vs co-treatment ratios (variant / baseline), within-sample ------

# 3a) Stratified by CSE
pub_cotreat_byCSE <- cotreat_pairs %>%
  dplyr::group_by(comp_base, comp_var, is_cse, Unit) %>%
  dplyr::summarise(
    n = sum(is.finite(ratio)),
    median_ratio = stats::median(ratio, na.rm = TRUE),
    ci = list(boot_ci(ratio)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(CI_low = purrr::map_dbl(ci, 1),
                CI_high = purrr::map_dbl(ci, 2)) %>%
  dplyr::select(-ci) %>%
  dplyr::arrange(comp_base, comp_var, is_cse, Unit)

readr::write_csv(pub_cotreat_byCSE, file.path(out_dir, "csv", "pub_Ratio_baseline_vs_cotreatments_byCSE_bootstrap.csv"))

# 3b) Pooled across CSE
pub_cotreat_pooled <- cotreat_pairs %>%
  dplyr::group_by(comp_base, comp_var, Unit) %>%
  dplyr::summarise(
    n = sum(is.finite(ratio)),
    median_ratio = stats::median(ratio, na.rm = TRUE),
    ci = list(boot_ci(ratio)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(CI_low = purrr::map_dbl(ci, 1),
                CI_high = purrr::map_dbl(ci, 2)) %>%
  dplyr::select(-ci) %>%
  dplyr::arrange(comp_base, comp_var, Unit)

readr::write_csv(pub_cotreat_pooled, file.path(out_dir, "csv", "pub_Ratio_baseline_vs_cotreatments_overall_bootstrap.csv"))

message("Done. IC50s: ", file.path(out_dir, "csv", "IC50_results.csv"),
        " — Plots in: ", file.path(out_dir, "plots"))