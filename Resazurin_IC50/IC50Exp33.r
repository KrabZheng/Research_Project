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
raw_path <- "Resazurin_IC50/Resazurin_Exp_33_Raw.txt"
cmp_path <- "Resazurin_IC50/Resazurin_Exp_33_Compound.txt"
out_dir  <- "Resazurin_IC50/resazurin_outputs"
dir.create(out_dir, showWarnings = FALSE)
dir.create(file.path(out_dir, "plots"), showWarnings = FALSE)
dir.create(file.path(out_dir, "csv"),   showWarnings = FALSE)

# Helper: parse the RAW block-structured txt ------------------------------
parse_raw_blocks <- function(path) {
  lines  <- readLines(path, warn = FALSE)
  hdr_ix <- grep(
    "^[^\\t]+_(A|B)_(Fluo|Abs|Abs_Background)\\t1\\t2\\t3\\t4\\t5\\t6\\t7\\t8\\t9\\t10\\t11\\t12",
    lines
  )
  stopifnot(length(hdr_ix) > 0)

  purrr::map_dfr(hdr_ix, function(i) {
    hdr <- lines[[i]]
    m <- stringr::str_match(hdr, "^([^\\t]+)_(A|B)_(Fluo|Abs|Abs_Background)")
    sample_id <- m[2]; plate <- m[3]; meas <- m[4]

    # the next 8 lines are levels 1..8 (no explicit level column in RAW)
    block <- lines[(i+1):(i+8)]

    # choose which 12 columns to read from each line based on meas
    start <- dplyr::case_when(
      meas == "Fluo"           ~ 2L,   # first 12 numbers
      meas == "Abs"            ~ 15L,  # next 12
      meas == "Abs_Background" ~ 28L,  # last 12
      TRUE ~ 2L
    )

    purrr::imap_dfr(block, function(ln, lvl) {
      parts <- stringr::str_split(ln, "\\t", simplify = TRUE)
      vals  <- suppressWarnings(as.numeric(parts[start:(start+11)]))
      tibble(sample_id, plate, meas, level = as.integer(lvl),
             col = 1:12, value = vals)
    })
  })
}

raw_long <- parse_raw_blocks(raw_path)

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
  mutate(block = ((col - 1L) %/% 3L) + 1L,
         rep   = ((col - 1L) %% 3L) + 1L,
         pos   = paste0(plate, block)) %>%
  left_join(dils, by = c("pos", "level", "plate", "block"))

# We only need fluorescence for viability normalisation/fits --------------
fluo <- dat %>%
  filter(meas == "Fluo") %>%
  dplyr::select(sample_id, plate, pos, block, rep, level, dose, Unit, Compound, value)

# Compute robust top/bottom per sample_id × pos ----------------------------
top_bottom_map <- fluo %>%
  dplyr::group_by(sample_id, pos) %>%
  dplyr::summarise(
    # TOP: take the mean at the *minimum dose* (usually vehicle);
    # if dose is missing, fall back to the maximum signal across all levels.
    top = {
      md <- suppressWarnings(min(dose, na.rm = TRUE))
      t1 <- mean(value[dose == md], na.rm = TRUE)
      if (!is.finite(t1)) suppressWarnings(max(value, na.rm = TRUE)) else t1
    },
    # BOTTOM: prefer min at levels >= 2; if not available, min across all levels.
    bottom = {
      b1 <- suppressWarnings(min(value[level >= 2], na.rm = TRUE))
      if (!is.finite(b1)) suppressWarnings(min(value, na.rm = TRUE)) else b1
    },
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    # Final fallbacks to avoid NAs propagating downstream
    top    = dplyr::coalesce(ifelse(is.finite(top),    top,    NA_real_), bottom),
    bottom = dplyr::coalesce(ifelse(is.finite(bottom), bottom, NA_real_), top)
  )

# Merge and compute normalised viability in [0,1] --------------------------
fluo_n <- fluo %>%
  dplyr::left_join(top_bottom_map, by = c("sample_id", "pos")) %>%
  dplyr::mutate(
    denom = pmax(top - bottom, .Machine$double.eps, na.rm = TRUE),
    viab  = (value - bottom) / denom,
    viab  = pmin(pmax(viab, 0), 1)
  )

# Quick sanity check
print(
  fluo_n %>%
    dplyr::summarise(
      n_rows   = dplyr::n(),
      n_viab_na = sum(!is.finite(viab)),
      n_top_na  = sum(!is.finite(top)),
      n_bot_na  = sum(!is.finite(bottom))
    )
)

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
  dfit <- df %>%
  dplyr::filter(level >= 2, is.finite(dose), dose > 0, is.finite(viab))
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

# IC50 table now includes SE and p_trend
ic50_tbl <- res_list %>%
  dplyr::transmute(sample_id, Compound, Unit, pos, fit_ok, p_trend, ic50) %>%
  tidyr::unnest(ic50) %>%
  dplyr::mutate(across(c(ED, SE, Lower, Upper), as.numeric)) %>%
  dplyr::arrange(sample_id, Compound)

readr::write_csv(ic50_tbl, file.path(out_dir, "csv", "IC50_results.csv"))

# --- Helpers for pairwise IC50 tests (log-scale delta method) ----------------
# --- Helpers for IC50 comparisons & pretty printing ---------------------------
ic50_compare <- function(e1, se1, e2, se2) {
  # Ratio = e1/e2 on original scale; p for log-difference (delta method)
  if (any(!is.finite(c(e1, se1, e2, se2))) || e1 <= 0 || e2 <= 0) {
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

ic50_tbl2 <- ic50_tbl %>%
  dplyr::mutate(base = tolower(gsub("cse$", "", sample_id, ignore.case = TRUE)),
                is_cse = grepl("cse$", sample_id, ignore.case = TRUE))

# (A) sample vs samplecse for each compound (one row per base×compound)
cmp_across <- ic50_tbl2 %>%
  dplyr::group_by(base, Compound) %>%
  dplyr::filter(any(is_cse) & any(!is_cse)) %>%
  dplyr::summarise(
    sample_id      = first(sample_id[!is_cse]),
    samplecse_id   = first(sample_id[is_cse]),
    unit           = dplyr::coalesce(first(Unit[!is_cse]), first(Unit[is_cse])),
    e_sample       = first(ED[!is_cse]),
    se_sample      = first(SE[!is_cse]),
    e_samplecse    = first(ED[is_cse]),
    se_samplecse   = first(SE[is_cse]),
    .groups = "drop"
  ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(ic = ic50_compare(e_samplecse, se_samplecse, e_sample, se_sample)) %>%
  tidyr::unnest(ic) %>%
  dplyr::rename(IC50_sample = e_sample, SE_sample = se_sample,
                IC50_sampleCSE = e_samplecse, SE_sampleCSE = se_samplecse) %>%
  dplyr::arrange(base, Compound)

readr::write_csv(cmp_across, file.path(out_dir, "csv", "IC50_compare_sample_vs_sampleCSE.csv"))

# --- Combined plots: EXACTLY TWO PER COMPOUND (non-CSE & CSE) -----------------

# Folders
dir.create(file.path(out_dir, "plots_by_arm", "nonCSE"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "plots_by_arm", "CSE"),    recursive = TRUE, showWarnings = FALSE)

# Helper to derive base id (pair sample with its "...cse")
base_id <- function(s) tolower(gsub("cse$", "", s, ignore.case = TRUE))

# Annotate res_list with pairing helpers
res_ann <- res_list %>%
  dplyr::mutate(
    base   = base_id(sample_id),
    is_cse = grepl("cse$", sample_id, ignore.case = TRUE)
  )

## --- Palette setup (CB-friendly + explicit overrides) --------------------
## Load your palettes.R if present; otherwise define a safe fallback pal_for()
try(source("palettes.R"), silent = TRUE)
if (!exists("pal_for")) {
  pal_for <- function(ids) {
    ids <- sort(unique(as.character(ids)))
    n   <- length(ids)
    cols <- if (requireNamespace("viridisLite", quietly = TRUE)) {
      viridisLite::viridis(n)
    } else {
      grDevices::hcl.colors(n, palette = "Dark 3")
    }
    stats::setNames(cols, ids)
  }
}

# All groups present in plots
all_groups <- sort(unique(res_ann$sample_id))

# Base palette (colour-blind friendly)
pal_map <- pal_for(all_groups); names(pal_map) <- all_groups
# After pal_map is created (and overrides applied), add:
status_map <- c(
  "P676"   = "former",
  "P679"   = "never",
  "P698" = "current"
)

# Explicit overrides (Paul Tol CB-safe hues)
# P676 = orange, P679 = blue, ASC390 = red, ASC397 = red-ish but distinct
overrides <- c(
  "P676"   = "#EE7733",  # orange
  "P679"   = "#0077BB",  # blue
  "P698" = "#d13714"  # red
)

# Apply same colours to “…cse” variants if they exist in your sample_ids
overrides2 <- c(overrides,
  stats::setNames(unname(overrides), paste0(names(overrides), "cse"))
)

ix <- intersect(names(overrides2), names(pal_map))
pal_map[ix] <- overrides2[ix]

# (Optional) sanity check:
# print(pal_map[c("P676","P676cse","P679","P679cse","ASC390","ASC390cse","ASC397","ASC397cse")])

# Helper to plot a single compound × arm (either non-CSE or CSE)
.plot_compound_arm <- function(sub, cmp_name, arm_label, out_subdir) {
  if (nrow(sub) == 0) return(invisible(NULL))

  # Gather predictions, raw points, and IC50s with sample metadata
  preds <- purrr::map_dfr(seq_len(nrow(sub)), function(i) {
    sub$pred[[i]] %>% dplyr::mutate(sample_id = sub$sample_id[[i]])
  })
  pts <- purrr::map_dfr(seq_len(nrow(sub)), function(i) {
    sub$points[[i]] %>% dplyr::mutate(sample_id = sub$sample_id[[i]])
  })
  vdat <- purrr::map_dfr(seq_len(nrow(sub)), function(i) {
    tibble::tibble(
      dose = suppressWarnings(as.numeric(sub$ic50[[i]][1, "ED", drop = TRUE])),
      sample_id = sub$sample_id[[i]]
    )
  }) %>% dplyr::filter(is.finite(dose))
  ids_present <- sort(unique(preds$sample_id))
  base_names  <- sub("cse$", "", ids_present, ignore.case = TRUE)
  st          <- unname(status_map[base_names])
  legend_labs <- ifelse(is.na(st) | st == "", base_names, paste0(base_names, " (", st, ")"))
  names(legend_labs) <- ids_present


  # Unit label (first non-NA)
  unit_lbl <- {
    uu <- sub$Unit[!is.na(sub$Unit)]
    if (length(uu)) uu[[1]] else NA_character_
  }

  # Build plot (CI ribbons + lines + points + IC50 vlines)
  p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = preds,
      ggplot2::aes(dose, ymin = lwr, ymax = upr, fill = sample_id),
      alpha = 0.15, show.legend = FALSE
    ) +
    ggplot2::geom_line(
      data = preds,
      ggplot2::aes(dose, fit, colour = sample_id),
      linewidth = 0.9
    ) +
    ggplot2::geom_point(
      data = pts,
      ggplot2::aes(dose, viab, colour = sample_id),
      size = 1.9, alpha = 0.75
    ) +
    ggplot2::geom_vline(
      data = vdat,
      ggplot2::aes(xintercept = dose, colour = sample_id),
      linetype = 2, linewidth = 0.6, alpha = 0.7
    ) +
    ggplot2::scale_x_log10(
      labels = scales::label_number(),
      expand = ggplot2::expansion(mult = c(0.02, 0.05))
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)
    ) +
    ggplot2::scale_colour_manual(
      values = pal_map[ids_present],
      breaks = ids_present,
      labels = legend_labs
    ) +
    ggplot2::scale_fill_manual(
      values = scales::alpha(pal_map[ids_present], 0.18),
      breaks = ids_present,
      labels = legend_labs,
      guide  = "none"
    ) +
    ggplot2::labs(
      title  = cmp_name,            # compound only
      x      = glue::glue("Concentration ({unit_lbl})"),
      y      = "Normalised viability",
      colour = "Sample"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title         = ggplot2::element_text(face = "bold"),
      panel.grid.minor.x = ggplot2::element_blank(),
      legend.key.height  = ggplot2::unit(10, "pt")
    )

  fn <- file.path(out_dir, "plots_by_arm", out_subdir,
                  glue::glue("{gsub('[^A-Za-z0-9]+','_', cmp_name)}_{out_subdir}.png"))
  ggplot2::ggsave(fn, p, width = 7.5, height = 5.0, dpi = 300)
  invisible(fn)
}

# Make exactly two plots per compound
compounds <- res_ann %>%
  dplyr::filter(fit_ok) %>%
  dplyr::pull(Compound) %>% unique() %>% sort()

for (cmp_name in compounds) {
  # non-CSE
  sub_non <- res_ann %>% dplyr::filter(Compound == cmp_name, fit_ok, !is_cse)
  .plot_compound_arm(sub_non, cmp_name, "non-CSE samples", "nonCSE")

  # CSE
  sub_cse <- res_ann %>% dplyr::filter(Compound == cmp_name, fit_ok,  is_cse)
  .plot_compound_arm(sub_cse, cmp_name, "CSE samples", "CSE")
}

message("Saved per-compound plots to: ",
        file.path(out_dir, "plots_by_arm", "nonCSE"), " and ",
        file.path(out_dir, "plots_by_arm", "CSE"))

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