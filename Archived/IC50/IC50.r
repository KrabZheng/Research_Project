##############################################################################
#   Resazurin IC50 analysis — final build (Exp30–36)
#   - Robust parsing, paired Base(−CSE) vs CSE(+CSE) analysis
#   - 4PL fits with nls/nlsLM + monotone-spline fallback
#   - Colour-blind friendly, per-patient outputs, Δlog2 forest panels
##############################################################################

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(purrr); library(tibble); library(forcats); library(glue); library(rlang)
  library(ggplot2); library(scales); library(broom); library(drc)
  library(ggbeeswarm)
})

# ----------------------------- CONFIG -------------------------------------- #
CONFIG <- list(
  out_dir = "PROCESSTODAY/IC50_output",
  norm_mode = "per_triplicate",     # "per_triplicate" | "per_plate" | "vehicle_row"
  vehicle_row_idx = 1L,             # only if norm_mode == "vehicle_row"
  min_points_for_fit = 6L,
  base_size = 11,
  bootstrap_B = 200L,
  error_band = "ci95",              # "sd" | "sem" | "ci95"
  seed = 123
)
if (!dir.exists(CONFIG$out_dir)) dir.create(CONFIG$out_dir, recursive = TRUE)

# ----------------------------- Helpers ------------------------------------- #
`%||%` <- function(a, b) if (!is.null(a)) a else b

safe_log10 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[!is.finite(x) | x <= 0] <- NA_real_
  log10(x)
}

strip_cse <- function(x) sub("(?i)cse$", "", x, perl = TRUE)  # "P679cse" -> "P679"
col_block  <- function(col_idx) ceiling(col_idx / 3)          # 1:3, 4:6, 7:9, 10:12

# Okabe–Ito colours (colour-blind friendly)
PALETTES <- list(
  oi = c("#0072B2", "#56B4E9", "#D55E00", "#E69F00",
         "#009E73", "#F0E442", "#CC79A7", "#999999"),
  groups_plot = c(`-CSE` = "#56B4E9", `+CSE` = "#D55E00")
)

mix_hex <- function(col, bg = "#FFFFFF", p = 0.3) {
  c1 <- grDevices::col2rgb(col)/255
  c2 <- grDevices::col2rgb(bg)/255
  out <- (1 - p) * c1 + p * c2
  grDevices::rgb(out[1], out[2], out[3])
}

# Paired colours per patient: "-CSE" darker, "+CSE" lighter
make_patient_pair_palette <- function(patients, levels_group = c("-CSE","+CSE")) {
  patients <- unique(patients)
  base_hues <- setNames(rep(PALETTES$oi, length.out = length(patients)), patients)
  pal <- c()
  for (p in patients) {
    base_col <- base_hues[[p]]
    pal[paste0(p, "|-CSE")] <- mix_hex(base_col, "#000000", p = 0.20)
    pal[paste0(p, "|+CSE")] <- mix_hex(base_col, "#FFFFFF", p = 0.35)
  }
  pal
}

palette_for <- function(levels) setNames(rep(PALETTES$oi, length.out = length(levels)), levels)

get_patients <- function(res) {
  res$wells %>%
    transmute(patient = strip_cse(sample_name)) %>%
    distinct() %>%
    filter(!is.na(patient), patient != "") %>%
    pull(patient) %>% unique() %>% sort()
}

filter_by_patient <- function(tbl, patient) tbl %>% filter(strip_cse(sample_name) == patient)
sanitize_filename <- function(x) gsub("[^A-Za-z0-9_\\-]+", "_", x)

# ------------------------ Read COMPOUND table ------------------------------ #
# Returns: experiment, plate_key(A1..B4), group(A/B), compound_name, unit, dose_idx, dose
read_compound_table <- function(path) {
  raw <- suppressMessages(readr::read_tsv(path, col_types = cols(.default = "?")))
  raw <- raw %>% mutate(IDN = toupper(trimws(as.character(ID))))

  meta_long <- raw %>%
    filter(IDN %in% c("EXPERIMENT","COMPOUND","UNIT","UNITS")) %>%
    pivot_longer(cols = -c(ID, IDN), names_to = "plate_key", values_to = "value") %>%
    mutate(IDN = recode(IDN, "UNITS" = "UNIT"))

  meta_wide <- meta_long %>% pivot_wider(names_from = IDN, values_from = value) %>%
    mutate(COMPOUND = trimws(as.character(COMPOUND))) %>%
    distinct(plate_key, .keep_all = TRUE)

  ladder <- raw %>%
    filter(grepl("^[0-9]+$", ID)) %>%
    mutate(dose_idx = as.integer(ID)) %>%
    pivot_longer(cols = -c(ID, IDN, dose_idx),
                 names_to = "plate_key", values_to = "dose_raw") %>%
    mutate(dose = suppressWarnings(readr::parse_number(as.character(dose_raw)))) %>%
    filter(!is.na(dose)) %>%
    transmute(plate_key, dose_idx, dose) %>%
    distinct(plate_key, dose_idx, .keep_all = TRUE)

  ladder %>%
    left_join(meta_wide %>% rename(experiment = EXPERIMENT,
                                   compound_name = COMPOUND,
                                   unit = UNIT),
              by = "plate_key") %>%
    mutate(
      group = substr(plate_key, 1, 1),
      unit  = if_else(is.na(unit) | unit == "", "µM", as.character(unit)),
      compound_name = if_else(is.na(compound_name) | compound_name == "",
                              plate_key, compound_name)
    ) %>%
    dplyr::select(experiment, plate_key, group, compound_name, unit, dose_idx, dose)
}

# ------------------------- Read SAMPLE layout ------------------------------ #
# Returns: col_idx (1..12), sample_name
read_sample_layout <- function(path) {
  raw <- suppressMessages(readr::read_tsv(path, col_types = cols(.default = "?")))
  cell_line_row <- which(raw$ID == "Compound"); stopifnot(length(cell_line_row) == 1)
  cell_lines <- as.list(raw[cell_line_row, -1]) %>% unlist(use.names = FALSE)
  layout_rows <- raw %>% filter(grepl("^[0-9]+$", ID))
  tibble(col_idx = as.integer(layout_rows$ID)) %>%
    mutate(sample_name = purrr::map_chr(col_idx, function(ci) {
      row_ind <- which(as.integer(layout_rows$ID) == ci)
      which_one <- which(as.integer(layout_rows[row_ind, -1]) == 1)
      if (length(which_one) == 1) cell_lines[[which_one]] else paste0("Block", col_block(ci))
    }))
}

# ---------------------------- Parse RAW plates ----------------------------- #
# Returns: file, plate_label, mode, row_idx, col_idx, readout
parse_raw_files <- function(paths, modes = c("Fluo")) {
  is_header <- function(line) grepl("\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\t12$", line)
  parse_one <- function(path) {
    lines <- readr::read_lines(path)
    header_idx <- which(vapply(lines, is_header, logical(1)))
    purrr::map_dfr(header_idx, function(hi) {
      header_line <- lines[hi]
      label <- strsplit(header_line, "\t", fixed = TRUE)[[1]][1]
      mode  <- sub(".*_([A-Za-z]+)$", "\\1", label)
      if (!(mode %in% modes)) return(tibble())
      label_core <- sub("_([A-Za-z]+)$", "", label)

      cand <- lines[(hi + 1):min(hi + 30, length(lines))]
      cand <- cand[grepl("[0-9]", cand)]
      if (!length(cand)) return(tibble())
      cand <- cand[seq_len(min(8, length(cand)))]

      df <- purrr::imap_dfr(cand, function(r, i) {
        bits <- strsplit(r, "\t", fixed = TRUE)[[1]]
        nums <- suppressWarnings(as.numeric(bits[-1]))
        nums <- nums[is.finite(nums)]
        if (length(nums) < 12) return(tibble())
        tibble(row_idx = i, vals = list(nums[1:12]))
      })
      if (!nrow(df)) return(tibble())

      mat <- do.call(rbind, lapply(df$vals, function(x) matrix(x, nrow = 1)))
      colnames(mat) <- as.character(seq_len(ncol(mat)))
      as_tibble(mat) |>
        mutate(row_idx = df$row_idx) |>
        pivot_longer(-row_idx, names_to = "col_idx", values_to = "readout") |>
        mutate(col_idx = as.integer(col_idx),
               file = path, plate_label = label_core, mode = mode)
    })
  }
  res <- purrr::map_dfr(paths, parse_one)
  if (!nrow(res)) stop("No data parsed from RAW files; check headers.", call. = FALSE)
  res
}

# -------------------- Annotate wells with layout & compounds --------------- #
annotate_layout <- function(raw_long, comp_tbl, sample_map = NULL) {
  req <- c("row_idx","col_idx","readout","plate_label")
  missing <- setdiff(req, names(raw_long))
  if (length(missing)) stop(glue("annotate_layout(): missing {toString(missing)}"), call. = FALSE)

  raw_long <- raw_long %>%
    mutate(
      group_key = case_when(
        grepl("_[AB]$", plate_label) ~ sub(".*_([AB])$", "\\1", plate_label),             # "..._A"
        grepl("^(A|B)[0-9]+$", plate_label) ~ sub("^((A|B)[0-9]+)$", "\\1", plate_label), # "A1"
        TRUE ~ NA_character_
      ),
      block = col_block(col_idx)
    )

  # Exp30/33 style (_A/_B; compounds via A1..A4/B1..B4)
  ab_style <- raw_long %>% filter(!is.na(group_key) & group_key %in% c("A","B") & grepl("_[AB]$", plate_label))
  ab_annot <- NULL
  if (nrow(ab_style)) {
    ab_annot <- ab_style %>%
      mutate(plate_key = paste0(group_key, block)) %>%
      left_join(comp_tbl, by = c("plate_key", "row_idx" = "dose_idx")) %>%
      mutate(sample_name = sub("_([AB])$", "", plate_label),
             dose_idx = as.integer(row_idx))
  }

  # Exp34–36 style (plate_label itself A1/B1...; Sample map supplies sample_name)
  ax_style <- raw_long %>% filter(!is.na(group_key) & grepl("^(A|B)[0-9]+$", group_key))
  ax_annot <- NULL
  if (nrow(ax_style)) {
    ax_annot <- ax_style %>%
      mutate(plate_key = group_key) %>%
      left_join(comp_tbl, by = c("plate_key", "row_idx" = "dose_idx")) %>%
      mutate(dose_idx = as.integer(row_idx))
    if (!is.null(sample_map) && nrow(sample_map)) {
      ax_annot <- ax_annot %>% left_join(sample_map, by = "col_idx")
    } else {
      ax_annot <- ax_annot %>% mutate(sample_name = paste0("Block", block))
    }
  }

  bind_rows(ab_annot, ax_annot) %>%
    mutate(
      dose = as.numeric(dose),
      unit = as.character(unit),
      sample_name  = sample_name %||% sub("_([AB])$", "", plate_label),
      compound_name = coalesce(na_if(compound_name, ""), plate_key)
    )
}

# ----------------------------- Normalisation ------------------------------- #
normalise_readouts <- function(df, norm_mode = CONFIG$norm_mode, vehicle_row_idx = CONFIG$vehicle_row_idx) {
  df <- df %>% filter(!is.na(readout))
  by_trip  <- c("file","plate_label","plate_key","sample_name","block")
  by_plate <- c("file","plate_label")

  scale_minmax <- function(x) {
    ymin <- suppressWarnings(min(x, na.rm = TRUE))
    ymax <- suppressWarnings(max(x, na.rm = TRUE))
    if (!is.finite(ymin) || !is.finite(ymax) || ymax <= ymin) rep(NA_real_, length(x)) else (x - ymin) / (ymax - ymin)
  }

  if (norm_mode == "per_triplicate") {
    df <- df %>% group_by(across(all_of(by_trip)))  %>% mutate(y_try = scale_minmax(readout)) %>% ungroup()
  } else if (norm_mode == "per_plate") {
    df <- df %>% group_by(across(all_of(by_plate))) %>% mutate(y_try = scale_minmax(readout)) %>% ungroup()
  } else {
    df <- df %>% group_by(across(all_of(by_trip)))  %>% mutate({
      base_vals <- readout[row_idx == vehicle_row_idx]
      base <- if (length(base_vals)) mean(base_vals, na.rm = TRUE) else NA_real_
      y_try <- if (!is.finite(base) || base == 0) NA_real_ else readout / base
    }) %>% ungroup()
  }

  df <- df %>% group_by(across(all_of(by_trip))) %>% mutate(all_na_group = all(!is.finite(y_try))) %>% ungroup()
  if (any(df$all_na_group)) {
    df <- df %>% group_by(across(all_of(by_plate))) %>% mutate(y_plate = scale_minmax(readout)) %>% ungroup() %>%
      mutate(y = ifelse(all_na_group, y_plate, y_try))
  } else df <- df %>% mutate(y = y_try)

  df %>% dplyr::select(-any_of(c("y_try","y_plate","all_na_group"))) %>% mutate(y = pmin(pmax(y, 0), 1))
}

# ---- Exp30 tweak: B2 uses B1's min/max; B4 uses B3's min/max (no swap) ---- #
override_minmax_B_pairs <- function(df) {
  by_keys <- c("file","plate_label","sample_name","plate_key")
  mm <- df %>%
    group_by(across(all_of(by_keys))) %>%
    summarise(min_r = suppressWarnings(min(readout, na.rm = TRUE)),
              max_r = suppressWarnings(max(readout, na.rm = TRUE)), .groups = "drop")

  donor_map <- tibble(plate_key = c("B2","B4"),
                      donor     = c("B1","B3"))

  ref <- mm %>%
    inner_join(donor_map, by = "plate_key") %>%
    transmute(file, plate_label, sample_name,
              plate_key = plate_key,
              donor_key = donor)

  donor_vals <- mm %>%
    transmute(file, plate_label, sample_name, donor_key = plate_key,
              ref_min = min_r, ref_max = max_r)

  df %>%
    left_join(ref,      by = c("file","plate_label","sample_name","plate_key")) %>%
    left_join(donor_vals, by = c("file","plate_label","sample_name","donor_key")) %>%
    left_join(mm %>% transmute(file, plate_label, sample_name, plate_key,
                               own_min = min_r, own_max = max_r),
              by = c("file","plate_label","sample_name","plate_key")) %>%
    mutate(min_use = if_else(plate_key %in% c("B2","B4"), ref_min, own_min),
           max_use = if_else(plate_key %in% c("B2","B4"), ref_max, own_max),
           y = (readout - min_use) / pmax(max_use - min_use, 1e-12),
           y = pmin(pmax(y, 0), 1)) %>%
    dplyr::select(-donor_key,-ref_min,-ref_max,-own_min,-own_max,-min_use,-max_use)
}

# --------------------------- Zero-dose jitter ------------------------------- #
jitter_zero_dose <- function(df) {
  overall_min <- suppressWarnings(min(df$dose[df$dose > 0 & is.finite(df$dose)], na.rm = TRUE))
  if (!is.finite(overall_min)) overall_min <- 1e-12
  df %>%
    group_by(plate_key) %>%
    mutate(
      min_pos = suppressWarnings(min(dose[dose > 0 & is.finite(dose)], na.rm = TRUE)),
      min_pos = ifelse(is.finite(min_pos), min_pos, overall_min),
      dose_adj = ifelse(!is.finite(dose) | dose <= 0, min_pos / 2, dose)
    ) %>% ungroup() %>% dplyr::select(-min_pos) %>% rename(dose_orig = dose, dose = dose_adj)
}

# --------------------------- Model fitting --------------------------------- #
fit_ll4 <- function(dat) {
  dat2 <- dat %>% filter(is.finite(dose), is.finite(y)) %>% arrange(dose)
  if (nrow(dat2) < 4L) return(NULL)
  if (n_distinct(dat2$dose) < 3L) return(NULL)
  eps <- 1e-6
  dat2$y <- pmin(pmax(dat2$y, eps), 1 - eps)
  safe_try <- function(expr) suppressWarnings(try(expr, silent = TRUE))

  fit <- safe_try(
    drc::drm(y ~ dose, data = dat2,
             fct = drc::LL.4(fixed = c(NA, 0, 1, NA)),
             logDose = 10,
             control = drc::drmc(constr = FALSE, errorm = FALSE, noMessage = TRUE),
             na.action = na.omit)
  )
  if (!inherits(fit, "try-error")) return(fit)

  dmin <- min(dat2$dose, na.rm = TRUE); dmax <- max(dat2$dose, na.rm = TRUE)
  start_ic50 <- 10^mean(log10(c(dmin, dmax))); start_slope <- 1

  if (requireNamespace("minpack.lm", quietly = TRUE)) {
    fit2 <- safe_try(
      minpack.lm::nlsLM(y ~ 1/(1 + (dose/ic50)^slope), data = dat2,
                        start = list(ic50 = start_ic50, slope = start_slope),
                        lower = c(ic50 = dmin/10, slope = 0.01),
                        upper = c(ic50 = dmax*10, slope = 10))
    )
    if (!inherits(fit2, "try-error")) return(fit2)
  }

  fit3 <- safe_try(
    nls(y ~ 1/(1 + (dose/ic50)^slope), data = dat2,
        start = list(ic50 = start_ic50, slope = start_slope),
        algorithm = "port",
        lower = c(ic50 = dmin/10, slope = 0.01),
        upper = c(ic50 = dmax*10, slope = 10),
        control = list(warnOnly = TRUE))
  )
  if (!inherits(fit3, "try-error")) return(fit3)

  NULL
}

predict_model <- function(model, newdose) {
  if (is.null(model)) return(rep(NA_real_, length(newdose)))
  if (inherits(model, "drc")) {
    as.numeric(predict(model, newdata = data.frame(dose = newdose)))
  } else if (inherits(model, "nls") || inherits(model, "nlsLM")) {
    cf <- coef(model); slope <- unname(cf[["slope"]]); ic50 <- unname(cf[["ic50"]])
    1/(1 + (newdose/ic50)^slope)
  } else rep(NA_real_, length(newdose))
}

bounded_ic50 <- function(fit, dmin, dmax, target = 0.5) {
  if (is.null(fit) || !is.finite(dmin) || !is.finite(dmax) || dmax <= dmin) return(NA_real_)
  lower <- max(dmin/10, 1e-12); upper <- dmax*10
  f <- function(logd) predict_model(fit, 10^logd) - target
  a <- log10(lower); b <- log10(upper)
  fa <- f(a); fb <- f(b)
  if (!is.finite(fa) || !is.finite(fb) || fa*fb > 0) return(NA_real_)
  root <- tryCatch(uniroot(f, c(a, b))$root, error = function(e) NA_real_)
  if (!is.finite(root)) NA_real_ else 10^root
}

# --- Monotone spline fallback & IC50 from spline ---
spline_model <- function(dat) {
  dd <- dat %>% filter(is.finite(dose), is.finite(y)) %>% arrange(dose) %>% distinct(dose, .keep_all = TRUE)
  if (nrow(dd) < 3) return(NULL)
  x <- log10(dd$dose); y <- pmin(pmax(dd$y, 1e-6), 1 - 1e-6)
  sf <- stats::splinefun(x, y, method = "monoH.FC")
  list(type = "spline", fun = sf, xrange = range(x), dmin = min(dd$dose), dmax = max(dd$dose))
}
is_spline <- function(x) is.list(x) && !is.null(x$type) && identical(x$type, "spline")
predict_spline <- function(sm, newdose) if (is.null(sm)) rep(NA_real_, length(newdose)) else sm$fun(log10(newdose))
ic50_from_spline <- function(sm, target = 0.5) {
  if (is.null(sm)) return(tibble(ic50 = NA_real_, ic50_flag = "nofit"))
  a <- sm$xrange[1]; b <- sm$xrange[2]
  fa <- sm$fun(a) - target; fb <- sm$fun(b) - target
  if (!is.finite(fa) || !is.finite(fb)) return(tibble(ic50 = NA_real_, ic50_flag = "nofit"))
  if (fa*fb > 0) {
    if (fa < 0 && fb < 0) return(tibble(ic50 = sm$dmax, ic50_flag = "right_censored"))
    if (fa > 0 && fb > 0) return(tibble(ic50 = sm$dmin, ic50_flag = "left_censored"))
    return(tibble(ic50 = NA_real_, ic50_flag = "nocross"))
  }
  root <- tryCatch(uniroot(function(z) sm$fun(z) - target, c(a,b))$root, error = function(e) NA_real_)
  tibble(ic50 = if (is.finite(root)) 10^root else NA_real_, ic50_flag = "ok")
}

ed50_ci <- function(fit) {
  out <- try(drc::ED(fit, 50, type = "absolute", interval = "delta"), silent = TRUE)
  if (inherits(out, "try-error")) return(c(NA_real_, NA_real_, NA_real_))
  c(out[1,"Estimate"], out[1,"Lower"], out[1,"Upper"])
}

# ------------------------------ Pipeline ----------------------------------- #
run_pipeline <- function(raw_paths, compound_path, sample_path = NULL, experiment_label = NA_character_) {
  message("→ Output dir: ", normalizePath(CONFIG$out_dir, winslash = "/"))
  message("Parsing RAW files …")
  raw_long <- parse_raw_files(raw_paths, modes = c("Fluo"))

  message("Reading compound table …")
  comp_tbl <- read_compound_table(compound_path)

  sample_map <- NULL
  if (!is.null(sample_path) && file.exists(sample_path)) {
    message("Reading sample layout …")
    sample_map <- read_sample_layout(sample_path)
  }

  message("Annotating wells …")
  wells <- annotate_layout(raw_long, comp_tbl, sample_map)

  message("Normalising readouts (", CONFIG$norm_mode, ") …")
  wells <- normalise_readouts(wells, norm_mode = CONFIG$norm_mode, vehicle_row_idx = CONFIG$vehicle_row_idx)

  # EXP30: apply B2←B1 and B4←B3 min–max override
  if (!is.null(experiment_label) && grepl("Exp30", experiment_label, ignore.case = TRUE)) {
    wells <- override_minmax_B_pairs(wells)
  }

  wells <- jitter_zero_dose(wells)
  wells <- wells %>% mutate(rep_id = paste0("r", ((col_idx - 1) %% 3) + 1))

  # Per-dose summaries (mean ± err across triplicates)
  dose_summary <- wells %>%
    group_by(experiment = experiment_label %||% first(experiment),
             plate_label, plate_key, sample_name, compound_name, unit,
             dose_idx, dose) %>%
    summarise(y_mean = mean(y, na.rm = TRUE),
              y_sd   = sd(y, na.rm = TRUE),
              n      = sum(is.finite(y)), .groups = "drop") %>%
    mutate(y_err = case_when(
      CONFIG$error_band == "sem"  ~ y_sd / sqrt(pmax(n, 1)),
      CONFIG$error_band == "ci95" ~ qt(0.975, pmax(n - 1, 1)) * y_sd / sqrt(pmax(n, 1)),
      TRUE                        ~ y_sd))

  message("Fitting dose–response models …")

  fits <- wells %>%
    group_by(experiment = experiment_label %||% first(experiment),
             plate_label, plate_key, sample_name, compound_name, unit) %>%
    group_modify(~ {
      replic <- .x %>% dplyr::select(dose_idx, dose, y)
      dat_for_fit <- replic %>% group_by(dose_idx, dose) %>%
        summarise(y = mean(y, na.rm = TRUE), .groups = "drop") %>% arrange(dose)
      dmin <- suppressWarnings(min(dat_for_fit$dose, na.rm = TRUE))
      dmax <- suppressWarnings(max(dat_for_fit$dose, na.rm = TRUE))

      mod <- fit_ll4(dat_for_fit)
      ic50_est <- bounded_ic50(mod, dmin, dmax, 0.5)
      ci_l <- NA_real_; ci_u <- NA_real_; flag <- "ok"; mtype <- "other"

      if (!is.null(mod)) {
        if (inherits(mod,"drc")) {
          mtype <- "drc"
          ci <- ed50_ci(mod)
          if (is.finite(ci[1])) { ic50_est <- ci[1]; ci_l <- ci[2]; ci_u <- ci[3] }
        } else if (inherits(mod,"nlsLM")) mtype <- "nlsLM" else if (inherits(mod,"nls")) mtype <- "nls"
      }

      sm_obj <- NULL
      if (!is.finite(ic50_est)) {
        sm_obj <- spline_model(dat_for_fit)
        ic <- ic50_from_spline(sm_obj)
        ic50_est <- ic$ic50; flag <- ic$ic50_flag; mtype <- "spline"
      }

      need_boot <- !is.finite(ci_l) || !is.finite(ci_u)
      if (need_boot && is.finite(ic50_est) && nrow(dat_for_fit) >= 4) {
        if (!is.null(CONFIG$seed)) set.seed(CONFIG$seed)
        B <- CONFIG$bootstrap_B %||% 200L
        boot_vals <- numeric(B)
        for (b in seq_len(B)) {
          boot_summary <- replic %>% group_by(dose_idx, dose) %>%
            summarise(y = mean(sample(y[is.finite(y)], size = n(), replace = TRUE)), .groups="drop") %>%
            arrange(dose)
          mb <- fit_ll4(boot_summary)
          v <- bounded_ic50(mb, dmin, dmax, 0.5)
          if (!is.finite(v)) {
            smb <- spline_model(boot_summary)
            vb <- ic50_from_spline(smb)$ic50
            boot_vals[b] <- if (is.finite(vb)) vb else NA_real_
          } else boot_vals[b] <- v
        }
        boot_vals <- boot_vals[is.finite(boot_vals)]
        if (length(boot_vals) >= 20) {
          qs <- stats::quantile(boot_vals, c(.025,.975), na.rm = TRUE, names = FALSE)
          ci_l <- qs[1]; ci_u <- qs[2]
        }
      }

      model_to_store <- if (mtype == "spline") sm_obj else mod
      tibble(model = list(model_to_store), model_type = mtype,
             ic50 = ic50_est, ic50_l = ci_l, ic50_u = ci_u, ic50_flag = flag)
    }) %>% ungroup() %>%
    mutate(logIC50 = safe_log10(ic50),
           logIC50_l = safe_log10(ic50_l),
           logIC50_u = safe_log10(ic50_u))

  # Replicate-specific fits (to show spread in boxplots)
  replicate_fits <- wells %>%
    group_by(experiment = experiment_label %||% first(experiment),
             plate_label, plate_key, sample_name, rep_id, compound_name, unit) %>%
    group_modify(~{
      dat <- .x %>% dplyr::select(dose_idx, dose, y) %>% arrange(dose)
      dmin <- suppressWarnings(min(dat$dose, na.rm = TRUE))
      dmax <- suppressWarnings(max(dat$dose, na.rm = TRUE))
      mod <- fit_ll4(dat)
      ic  <- bounded_ic50(mod, dmin, dmax, 0.5)
      mty <- if (!is.null(mod)) (if (inherits(mod,"drc")) "drc" else if (inherits(mod,"nlsLM")) "nlsLM" else "nls") else "none"
      if (!is.finite(ic)) { sm <- spline_model(dat); ic <- ic50_from_spline(sm)$ic50; mod <- sm; mty <- "spline" }
      tibble(model = list(mod), model_type = mty, ic50 = ic)
    }) %>% ungroup() %>%
    mutate(logIC50 = safe_log10(ic50), logIC50_l = NA_real_, logIC50_u = NA_real_)

  list(wells = wells, dose_summary = dose_summary, fits = fits, replicate_fits = replicate_fits)
}

# ------------------------------ Statistics --------------------------------- #
boot_ic50_log10 <- function(df_repl, B = CONFIG$bootstrap_B %||% 400L) {
  if (!is.null(CONFIG$seed)) set.seed(CONFIG$seed)
  make_summary <- function(d) d %>% group_by(dose_idx, dose) %>% summarise(y = mean(y, na.rm = TRUE), .groups="drop") %>% arrange(dose)
  summ <- make_summary(df_repl)
  dmin <- suppressWarnings(min(summ$dose, na.rm = TRUE))
  dmax <- suppressWarnings(max(summ$dose, na.rm = TRUE))
  m0   <- fit_ll4(summ)
  ic0  <- bounded_ic50(m0, dmin, dmax, 0.5)
  if (!is.finite(ic0)) { sm0 <- spline_model(summ); ic0 <- ic50_from_spline(sm0)$ic50 }
  est <- safe_log10(ic0)

  draws <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    bs <- df_repl %>% group_by(dose_idx, dose) %>%
      summarise(y = mean(sample(y[is.finite(y)], size = n(), replace = TRUE)), .groups="drop") %>% arrange(dose)
    mb <- fit_ll4(bs)
    icb <- bounded_ic50(mb, dmin, dmax, 0.5)
    if (!is.finite(icb)) { smb <- spline_model(bs); icb <- ic50_from_spline(smb)$ic50 }
    draws[b] <- safe_log10(icb)
  }
  ci <- if (sum(is.finite(draws)) >= 40) stats::quantile(draws, c(.025,.975), na.rm = TRUE, names = FALSE) else c(NA_real_, NA_real_)
  list(est = est, l = ci[1], u = ci[2], draws = draws)
}

paired_delta_boot <- function(df_pair, B = CONFIG$bootstrap_B %||% 400L) {
  if (!is.null(CONFIG$seed)) set.seed(CONFIG$seed)
  d_base <- df_pair %>% filter(Group == "Base")
  d_cse  <- df_pair %>% filter(Group == "CSE")
  a <- boot_ic50_log10(d_base, B); c <- boot_ic50_log10(d_cse, B)
  ok <- is.finite(a$draws) & is.finite(c$draws)
  d_draws <- c$draws[ok] - a$draws[ok]
  d_ci <- if (length(d_draws) >= 40) stats::quantile(d_draws, c(.025,.975), na.rm = TRUE, names = FALSE) else c(NA_real_, NA_real_)
  list(est = c$est - a$est, l = d_ci[1], u = d_ci[2], fold = 10^(c$est - a$est))   # log10 scale
}

compare_ic50_within_patient_by_compound <- function(fits) {
  dat <- fits %>%
    transmute(experiment, compound_name,
              patient_id = strip_cse(sample_name),
              cond = if_else(grepl("cse$", tolower(sample_name)), "CSE", "Base"),
              ic50, ic50_l, ic50_u, logIC50) %>%
    filter(!is.na(patient_id))

  wide <- dat %>% distinct(experiment, compound_name, patient_id, cond, .keep_all = TRUE) %>%
    pivot_wider(names_from = cond, values_from = c(ic50, ic50_l, ic50_u, logIC50), names_sep = ".")

  wide %>%
    mutate(
      has_pair = is.finite(ic50.Base) & is.finite(ic50.CSE),
      FC       = if_else(has_pair, ic50.CSE / ic50.Base, NA_real_),
      log10FC  = if_else(has_pair, log10(FC), NA_real_),

      se_log_base = if_else(is.finite(ic50_u.Base) & is.finite(ic50_l.Base),
                            (log(ic50_u.Base) - log(ic50_l.Base))/(2*1.96), NA_real_),
      se_log_cse  = if_else(is.finite(ic50_u.CSE) & is.finite(ic50_l.CSE),
                            (log(ic50_u.CSE) - log(ic50_l.CSE))/(2*1.96), NA_real_),
      se_log_ratio = sqrt(se_log_base^2 + se_log_cse^2),
      logFC        = ifelse(has_pair, log(ic50.CSE/ic50.Base), NA_real_),
      logFC_l      = ifelse(has_pair & is.finite(se_log_ratio), logFC - 1.96*se_log_ratio, NA_real_),
      logFC_u      = ifelse(has_pair & is.finite(se_log_ratio), logFC + 1.96*se_log_ratio, NA_real_),
      FC_l         = ifelse(is.finite(logFC_l), exp(logFC_l), NA_real_),
      FC_u         = ifelse(is.finite(logFC_u), exp(logFC_u), NA_real_)
    ) %>%
    dplyr::select(experiment, compound_name, patient_id, ic50.Base, ic50.CSE, FC, FC_l, FC_u, log10FC) %>%
    arrange(experiment, compound_name, patient_id)
}

# ------------------------------ Plotting ----------------------------------- #
plot_dose_curves <- function(dose_summary, fits, facet = c("by_compound","by_sample")) {
  facet <- match.arg(facet)

  # Predictions grid from stored models
  pred <- purrr::pmap_dfr(
    list(fits$model, fits$model_type, fits$plate_key, fits$sample_name, fits$compound_name, fits$unit),
    function(m, mtype, plate_key, sample_name, compound_name, unit) {
      obs <- dose_summary %>% filter(plate_key == !!plate_key, sample_name == !!sample_name)
      dmin <- suppressWarnings(min(obs$dose, na.rm = TRUE))
      dmax <- suppressWarnings(max(obs$dose, na.rm = TRUE))
      if (!is.finite(dmin) || !is.finite(dmax) || dmax <= dmin) return(tibble())
      grid <- tibble(dose = 10^seq(log10(dmin), log10(dmax), length.out = 200))
      yhat <- tryCatch({
        if (is.null(m)) rep(NA_real_, nrow(grid))
        else if (identical(mtype,"spline") && is_spline(m)) as.numeric(predict_spline(m, grid$dose))
        else if (inherits(m,"drc") || inherits(m,"nls") || inherits(m,"nlsLM")) predict_model(m, grid$dose)
        else rep(NA_real_, nrow(grid))
      }, error = function(e) rep(NA_real_, nrow(grid)))
      tibble(plate_key, sample_name, compound_name, unit, dose = grid$dose, yhat = as.numeric(yhat))
    }
  ) %>% filter(is.finite(yhat))

  ds <- dose_summary %>%
    mutate(patient = strip_cse(sample_name),
           Group   = ifelse(grepl("cse$", tolower(sample_name)), "+CSE", "-CSE"),
           pair_key = paste(patient, Group, sep = "|"))
  pred <- pred %>%
    mutate(patient = strip_cse(sample_name),
           Group   = ifelse(grepl("cse$", tolower(sample_name)), "+CSE", "-CSE"),
           pair_key = paste(patient, Group, sep = "|"))

  pal_pairs <- make_patient_pair_palette(unique(ds$patient))

  g <- ggplot(ds, aes(dose, y_mean,
                      colour = pair_key, fill = pair_key, linetype = Group)) +
    geom_ribbon(aes(ymin = pmax(y_mean - y_err, 0), ymax = pmin(y_mean + y_err, 1)),
                alpha = 0.15, colour = NA, show.legend = FALSE) +
    geom_point(size = 1.7, alpha = 0.95) +
    geom_line(data = pred, aes(dose, yhat), linewidth = 0.7) +
    scale_x_log10(labels = label_number(), breaks = log_breaks(n = 4)) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    scale_colour_manual(values = pal_pairs, guide = guide_legend(title = "Patient · Condition")) +
    scale_fill_manual(values = pal_pairs, guide = "none") +
    scale_linetype_manual(values = c(`-CSE` = "solid", `+CSE` = "dashed"), guide = "none") +
    labs(x = "Dose", y = "Normalised viability") +
    theme_bw(base_size = CONFIG$base_size)

  if (facet == "by_compound") {
    g <- g + facet_wrap(~ compound_name + unit, scales = "free_x")
  } else {
    g <- g + facet_wrap(~ sample_name, scales = "free_x")
  }
  g
}

plot_ic50_box <- function(rep_fits) {
  dat <- rep_fits %>%
    filter(is.finite(logIC50)) %>%
    mutate(patient = strip_cse(sample_name),
           Group   = ifelse(grepl("cse$", tolower(sample_name)), "+CSE", "-CSE"))

  ggplot(dat, aes(Group, logIC50, fill = Group, colour = Group)) +
    ggbeeswarm::geom_quasirandom(width = 0.15, alpha = 0.9, size = 2, show.legend = FALSE) +
    geom_boxplot(width = 0.60, alpha = 0.25, outlier.shape = NA) +
    scale_colour_manual(values = PALETTES$groups_plot, guide = "none") +
    scale_fill_manual(values = PALETTES$groups_plot, guide = "none") +
    labs(x = NULL, y = expression(log[10](IC[50]))) +
    facet_grid(patient ~ compound_name, scales = "free_y", switch = "y") +
    theme_bw(base_size = CONFIG$base_size)
}

plot_ic50_pairs <- function(fits) {
  dat <- fits %>%
    filter(is.finite(logIC50)) %>%
    mutate(patient = strip_cse(sample_name),
           Group   = ifelse(grepl("cse$", tolower(sample_name)), "+CSE", "-CSE"),
           pair_key = paste(patient, Group, sep = "|"))

  keep <- dat %>% count(compound_name, patient, Group) %>%
    pivot_wider(names_from = Group, values_from = n, values_fill = 0) %>%
    filter(`-CSE` > 0, `+CSE` > 0) %>% dplyr::select(compound_name, patient)

  dat2 <- dat %>% inner_join(keep, by = c("compound_name","patient"))
  if (nrow(dat2) == 0) return(ggplot() + theme_bw() + labs(caption = "No paired −CSE/+CSE"))

  pal_pairs <- make_patient_pair_palette(unique(dat2$patient))

  ggplot(dat2, aes(Group, logIC50, group = patient)) +
    geom_line(alpha = 0.35, colour = "grey60") +
    geom_linerange(aes(ymin = logIC50_l, ymax = logIC50_u, colour = pair_key),
                   alpha = 0.9, linewidth = 0.6, na.rm = TRUE) +
    geom_point(aes(colour = pair_key), size = 2) +
    scale_colour_manual(values = pal_pairs, guide = guide_legend(title = "Patient · Condition")) +
    facet_wrap(~ compound_name, scales = "free_y", drop = FALSE) +
    labs(x = NULL, y = expression(log[10](IC[50]))) +
    theme_bw(base_size = CONFIG$base_size)
}

# Single-panel forest (per patient): Δlog2(IC50), coloured by compound
forest_one_patient_panel <- function(res, B = CONFIG$bootstrap_B %||% 400L, patient) {
  wells <- res$wells %>%
    mutate(Group = ifelse(grepl("cse$", tolower(sample_name)), "CSE", "Base"),
           patient = strip_cse(sample_name)) %>%
    filter(patient == !!patient)

  keep <- wells %>% distinct(compound_name, Group) %>% count(compound_name) %>% filter(n == 2L) %>% dplyr::select(compound_name)
  pairs <- wells %>% semi_join(keep, by = "compound_name") %>%
    dplyr::select(compound_name, dose_idx, dose, y, Group)
  if (nrow(pairs) == 0) return(list(plot = ggplot() + theme_void(), stats = tibble()))

  stats <- pairs %>%
    group_by(compound_name) %>%
    group_modify(~{
      out <- paired_delta_boot(.x, B = B)  # log10 scale
      conv <- 1 / log10(2)
      tibble(delta_log2 = out$est * conv, ci_l2 = out$l * conv, ci_u2 = out$u * conv)
    }) %>% ungroup() %>%
    filter(is.finite(delta_log2)) %>%
    mutate(compound_name = as.character(compound_name)) %>%
    arrange(delta_log2) %>%
    mutate(compound_name = factor(compound_name, levels = unique(compound_name)))

  if (nrow(stats) == 0) return(list(plot = ggplot() + theme_void(), stats = tibble()))
  pal_comp <- palette_for(levels(stats$compound_name))

  p <- ggplot(stats, aes(x = delta_log2, y = compound_name, colour = compound_name)) +
    geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey70") +
    geom_errorbarh(aes(xmin = ci_l2, xmax = ci_u2), height = 0, linewidth = 0.9, alpha = 0.95) +
    geom_point(size = 2.2) +
    scale_colour_manual(values = pal_comp, guide = guide_legend(title = "Compound")) +
    labs(x = expression(Delta*log[2](IC[50])*" (+CSE − −CSE)"),
         y = "Compound",
         subtitle = "Right = higher IC50 under CSE (resistance); left = sensitised") +
    theme_bw(base_size = CONFIG$base_size)

  list(plot = p, stats = stats)
}

# ------------------------------ I/O helpers -------------------------------- #
write_csv_safe <- function(x, path) { readr::write_csv(x, path); message("✓ Written: ", path) }
ggsave_safe <- function(plot, filename, width = 8, height = 5.5) {
  fp <- file.path(CONFIG$out_dir, filename)
  ggplot2::ggsave(fp, plot = plot, width = width, height = height, dpi = 300)
  message("✓ Saved plot: ", fp)
}

# ------------------------------ Paths -------------------------------------- #
paths <- list(
  exp30 = list(
    raw = c("PROCESSTODAY/Resazurin_Exp_30_CSE_H2O2_KBrO3_Menadione_RSL3_RSL3CSE_Erastin_ErastinCSE.txt"),
    compound = "PROCESSTODAY/Resazurin_Exp_30_Compound.txt",
    sample = NULL,
    label = "Resazurin_Exp30"
  ),
  exp33 = list(
    raw = c("PROCESSTODAY/Resazurin_Exp_33_BSO_BSOCSE_Rotenone_RotenoneCSE.txt"),
    compound = "PROCESSTODAY/Resazurin_Exp_33_Compound.txt",
    sample = NULL,
    label = "Resazurin_Exp33"
  ),
  exp34 = list(
    raw = c("PROCESSTODAY/Resazurin_Exp_34_ASC397 ASC397cse ASC399 ASC399cse Rotenone H2O2.txt"),
    compound = "PROCESSTODAY/Resazurin_Exp_34_Compound.txt",
    sample   = "PROCESSTODAY/Resazurin_Exp_34_Sample.txt",
    label = "Resazurin_Exp34"
  ),
  exp35 = list(
    raw = c("PROCESSTODAY/Resazurin_Exp_35_ASC397 ASC397cse ASC399 ASC399cse KBrO3.txt"),
    compound = "PROCESSTODAY/Resazurin_Exp_35_Compound.txt",
    sample   = "PROCESSTODAY/Resazurin_Exp_35_Sample.txt",
    label = "Resazurin_Exp35"
  ),
  exp36 = list(
    raw = c("PROCESSTODAY/Resazurin_Exp_36_ASC397 ASC397cse ASC399 ASC399cse CSE.txt"),
    compound = "PROCESSTODAY/Resazurin_Exp_36_Compound.txt",
    sample   = "PROCESSTODAY/Resazurin_Exp_36_Sample.txt",
    label = "Resazurin_Exp36"
  )
)

# ------------------------------ Runner ------------------------------------- #
run_and_save <- function(p) {
  res <- run_pipeline(raw_paths        = p$raw,
                      compound_path    = p$compound,
                      sample_path      = p$sample,
                      experiment_label = p$label)

  # Tables
  write_csv_safe(res$fits,         file.path(CONFIG$out_dir, paste0(p$label, "_IC50_fits.csv")))
  write_csv_safe(res$dose_summary, file.path(CONFIG$out_dir, paste0(p$label, "_DoseSummary.csv")))

  # Whole-experiment plots
  g1 <- plot_dose_curves(res$dose_summary, res$fits, facet = "by_compound")
  ggsave_safe(g1, paste0(p$label, "_DoseCurves_byCompound.png"), width = 9.5, height = 6.2)

  g2 <- plot_ic50_pairs(res$fits)
  ggsave_safe(g2, paste0(p$label, "_IC50_Paired.png"))

  g3 <- plot_ic50_box(res$replicate_fits)
  ggsave_safe(g3, paste0(p$label, "_IC50_Box_byCompound.png"), width = 9.5, height = 6.5)

  # Overall Δlog2 forest (all patients, facetted by compound names)
  # (Keep the earlier multi-facet version if you like; below we do per-patient panels.)
  # fw_all <- forest_within_patient(res, B = CONFIG$bootstrap_B)  # optional

  # Per-patient set
  pats <- get_patients(res)
  for (pat in pats) {
    pat_tag <- sanitize_filename(toupper(pat))

    fits_p <- res$fits         %>% filter_by_patient(patient = pat)
    ds_p   <- res$dose_summary %>% filter_by_patient(patient = pat)

    gdc <- plot_dose_curves(ds_p, fits_p, facet = "by_compound")
    ggsave_safe(gdc, paste0(p$label, "_PAT_", pat_tag, "_DoseCurves_byCompound.png"),
                width = 9.5, height = 6.2)

    gpairs <- plot_ic50_pairs(fits_p)
    ggsave_safe(gpairs, paste0(p$label, "_PAT_", pat_tag, "_IC50_Paired.png"))

    rep_fits_p <- res$replicate_fits %>% filter_by_patient(patient = pat)
    gbox <- plot_ic50_box(rep_fits_p)
    ggsave_safe(gbox, paste0(p$label, "_PAT_", pat_tag, "_IC50_Box_byCompound.png"),
                width = 7.5, height = 5.8)

    fw1 <- forest_one_patient_panel(res, B = CONFIG$bootstrap_B, patient = pat)
    if (nrow(fw1$stats)) {
      ggsave_safe(fw1$plot, paste0(p$label, "_PAT_", pat_tag, "_Forest_DeltaIC50.png"),
                  width = 8, height = 6)
      write_csv_safe(fw1$stats, file.path(CONFIG$out_dir, paste0(p$label, "_PAT_", pat_tag, "_Forest_DeltaIC50.csv")))
    } else {
      message("  ↳ ", pat, ": no −CSE/+CSE pair → forest skipped")
    }

    stats_p <- compare_ic50_within_patient_by_compound(fits_p)
    if (nrow(stats_p)) {
      write_csv_safe(stats_p, file.path(CONFIG$out_dir, paste0(p$label, "_PAT_", pat_tag, "_Stats_within_patient_byCompound.csv")))
    }
  }

  # Whole-experiment paired stats
  stats <- compare_ic50_within_patient_by_compound(res$fits)
  write_csv_safe(stats, file.path(CONFIG$out_dir, paste0(p$label, "_Stats_within_patient_byCompound.csv")))

  invisible(res)
}

# ------------------------------ Run ---------------------------------------- #
# Example:
# source("IC50_resazurin_final.R")
# res30 <- run_and_save(paths$exp30)
# res33 <- run_and_save(paths$exp33)
# res34 <- run_and_save(paths$exp34)
# res35 <- run_and_save(paths$exp35)
# res36 <- run_and_save(paths$exp36)
