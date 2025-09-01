##############################################################################
#   IC50 analysis — multi-plate, robust (v17 math, CI ribbons)     2025-08-29
##############################################################################

## ---- Config & paths ---------------------------------------------------------
dir_in           <- "Resazurin_IC50"
out_dir          <- file.path(dir_in, "ic50_outputs")
plots_dir        <- file.path(out_dir, "plots")
plots_comb_dir   <- file.path(out_dir, "plots_combined")

cfg <- list(
  seed               = 1L,
  min_unique_doses   = 4,            # require ≥4 unique positive doses for a fit
  slope_bounds       = c(0.2, 20),   # acceptable |b| range for QC
  aicc_tie_tol       = 2,            # ΔAICc ≤ 2 => prefer simpler model
  B_ci               = 1200,         # bootstrap draws for CI (keep reasonable for multi-plate)
  b_draw_bounds      = c(0.1, 50),   # sanity for bootstrap draws
  e_bounds_factor    = c(0.25, 4),   # e draws confined inside tested range × factors
  loess_span         = 0.8,
  unit_default       = "µM"
)

channel_mode     <- "Fluo"           # "Fluo" or "AbsMinusBG"
exclude_compounds<- c("Rotenone")    # case-insensitive
font_family      <- "Arial"
save_svg         <- TRUE
also_png         <- TRUE

set.seed(cfg$seed)

## ---- Libraries --------------------------------------------------------------
suppressPackageStartupMessages({
  library(MASS)        # keep MASS before dplyr to avoid select()/filter() masks
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(purrr); library(tibble); library(drc); library(ggplot2); library(glue)
  library(scales)
})

## ---- Utils ------------------------------------------------------------------
`%ni%` <- Negate(`%in%`)
ensure_dir <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE)
has_svglite <- requireNamespace("svglite", quietly = TRUE)

base_theme <- theme_minimal(base_size = 12) +
  theme(
    text = element_text(family = font_family, size = 12),
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank() # hide x gridlines
  )

save_plot <- function(p, path_base, w = 6.0, h = 4.5){
  if (is.null(p)) return(invisible(FALSE))
  if (save_svg) {
    if (has_svglite) try(ggsave(paste0(path_base, ".svg"), p, width = w, height = h, device = svglite::svglite), silent = TRUE)
    else             try(ggsave(paste0(path_base, ".svg"), p, width = w, height = h), silent = TRUE)
  }
  if (isTRUE(also_png)) try(ggsave(paste0(path_base, ".png"), p, width = w, height = h, dpi = 220), silent = TRUE)
  invisible(TRUE)
}

make_pd <- function(S, eps = 1e-8){
  S <- as.matrix(S); S <- (S + t(S))/2
  ev <- eigen(S, symmetric = TRUE)
  ev$values[ev$values < eps] <- eps
  ev$vectors %*% diag(ev$values, nrow = length(ev$values)) %*% t(ev$vectors)
}

quiet_ED_delta <- function(fit, p = 50){
  ed <- NULL
  ok <- try(suppressWarnings(utils::capture.output(
    ed <- drc::ED(fit, p, interval = "delta")
  )), silent = TRUE)
  if (inherits(ok, "try-error")) return(NULL)
  ed
}

## ---- v17 concentration standardiser ----------------------------------------
.standardise_conc <- function(x, unit){
  u <- tolower(gsub("\\s+", "", unit %||% ""))
  u <- gsub("μ","u",u)
  res <- list(value = x, unit_std = unit, ok = TRUE, convertible = FALSE)
  if (u %in% c("m","mol/l"))                 { res$value <- x*1e6; res$unit_std <- "µM"; res$convertible <- TRUE; return(res) }
  if (u %in% c("mm","mmol/l"))               { res$value <- x*1e3; res$unit_std <- "µM"; res$convertible <- TRUE; return(res) }
  if (u %in% c("um","µm","umol/l","µmol/l")) { res$value <- x;     res$unit_std <- "µM"; res$convertible <- TRUE; return(res) }
  if (u %in% c("nm","nmol/l"))               { res$value <- x*1e-3;res$unit_std <- "µM"; res$convertible <- TRUE; return(res) }
  if (u %in% c("pm","pmol/l"))               { res$value <- x*1e-6;res$unit_std <- "µM"; res$convertible <- TRUE; return(res) }

  if (u %in% c("%","percent","%v/v","v/v%","%w/v","w/v%")) { res$value <- x; res$unit_std <- "%"; return(res) }
  if (u %in% c("mg/ml","mgperml"))                         { res$value <- x; res$unit_std <- "mg/mL"; return(res) }
  if (u %in% c("ug/ml","µg/ml","mcg/ml","ugperml","µgperml")) { res$value <- x; res$unit_std <- "µg/mL"; return(res) }
  if (u %in% c("ng/ml","ngperml"))                         { res$value <- x; res$unit_std <- "ng/mL"; return(res) }

  res$ok <- FALSE
  res
}

convert_conc_to_std <- function(df){
  tmp <- purrr::pmap(list(df$Concentration, df$Unit), .standardise_conc)
  df$Dose              <- vapply(tmp, function(z) z$value, numeric(1))
  df$UnitStd           <- vapply(tmp, function(z) z$unit_std, character(1))
  df$UnitOK            <- vapply(tmp, function(z) z$ok, logical(1))
  df$ConvertibleTo_uM  <- vapply(tmp, function(z) isTRUE(z$convertible), logical(1))
  df
}

## ---- Plate readers (layout + sample map + raw) ------------------------------
build_layout <- function(path){
  raw <- readr::read_delim(path, "\t", col_names = FALSE, trim_ws = TRUE, show_col_types = FALSE)
  ids   <- as.character(raw[raw$X1 == "ID",       -1])
  names <- as.character(raw[raw$X1 == "Compound", -1])
  units <- as.character(raw[raw$X1 == "Unit",     -1])
  dose_rows <- raw %>% dplyr::filter(stringr::str_detect(X1, "^[0-9]+$"))
  purrr::map_dfr(seq_along(ids), function(j){
    purrr::map_dfr(seq_len(nrow(dose_rows)), function(i){
      tibble::tibble(CompoundID=ids[j], CompoundName=names[j],
                     Row=LETTERS[as.integer(dose_rows$X1[i])],
                     Concentration=as.numeric(dose_rows[i, j+1][[1]]), Unit=units[j])
    })
  })
}

get_group_names <- function(path){
  sraw <- readr::read_delim(path, "\t", col_names = FALSE, trim_ws = TRUE, show_col_types = FALSE)
  as.character(sraw[sraw$X1 == "Compound", -1])
}

get_viability_one_plate <- function(raw_path){
  base      <- sub("_Raw\\d*\\.txt$", "", raw_path)
  cmp_path  <- paste0(base, "_Compound.txt")
  smp_path  <- paste0(base, "_Sample.txt")
  if (!file.exists(cmp_path)) stop("Compound file missing for ", basename(raw_path))

  # --- read compound layout, compute slots-per-set for column→slot mapping
  layout_long  <- build_layout(cmp_path)
  slots_info <- layout_long %>%
    dplyr::distinct(CompoundID) %>%
    dplyr::mutate(
      SetLetter  = substr(CompoundID, 1, 1),
      SlotIndex  = suppressWarnings(as.integer(gsub("^[A-Z]", "", CompoundID)))
    ) %>%
    dplyr::group_by(SetLetter) %>%
    dplyr::summarise(n_slots = dplyr::n_distinct(SlotIndex), .groups = "drop")

  # if we have a Sample.txt, use it; otherwise we'll parse names from headers
  have_sample_map <- file.exists(smp_path)
  if (have_sample_map) {
    group_names  <- get_group_names(smp_path)
    reps_per_grp <- 12 / length(group_names)
  }

  all_lines <- readr::read_lines(raw_path)
  rows <- which(stringr::str_detect(all_lines, stringr::regex("_fluo", ignore_case = TRUE)))
  rows <- rows[rows + 8 <= length(all_lines)]
  if (!length(rows)) stop("No \"_Fluo\" blocks in ", basename(raw_path))

  one_block <- function(i){
    plate <- suppressMessages(readr::read_delim(
      I(paste(all_lines[i:(i+8)], collapse = "\n")),
      "\t", trim_ws = TRUE, show_col_types = FALSE
    ))
    fluo_start <- which(stringr::str_detect(names(plate), stringr::regex("_fluo$", ignore_case = TRUE)))[1]
    if (is.na(fluo_start)) return(NULL)

    # header token (e.g., "ASC390_A_Fluo" OR "A1_Fluo")
    header_token <- stringr::str_remove(names(plate)[fluo_start], "_[Ff]luo$")
    # parse set/sample from token
    # 1) "SampleID_SET" (e.g. ASC390_A)
    m1 <- stringr::str_match(header_token, "^([A-Za-z0-9]+)_([A-H])$")
    # 2) "A1" / "B3" (set + slot index shown, but we will compute slot by column anyway)
    m2 <- stringr::str_match(header_token, "^([A-H])(\\d+)$")

    if (!is.na(m1[1,1])) {
      sample_id  <- m1[1,2]
      set_letter <- m1[1,3]
    } else if (!is.na(m2[1,1])) {
      sample_id  <- NA_character_                     # no sample in header
      set_letter <- m2[1,2]
    } else {
      # last resort: try to find trailing "_A" etc.
      set_letter <- stringr::str_match(header_token, "([A-H])$")[,2]
      if (is.na(set_letter)) set_letter <- "A"
      sample_id <- NA_character_
    }

    sig <- plate[, (fluo_start + 1):(fluo_start + 12)]
    names(sig) <- sprintf("%02d", 1:12)

    df <- sig %>%
      dplyr::mutate(Row = LETTERS[1:8]) %>%
      tidyr::pivot_longer(-Row, names_to = "Col", values_to = "Signal",
                          names_transform = list(Col = as.integer)) %>%
      dplyr::mutate(
        SetLetter = set_letter,
        # map column → slot index using number of slots for this set
        n_slots   = dplyr::coalesce(slots_info$n_slots[match(SetLetter, slots_info$SetLetter)], 1L),
        rep_block = pmax(1L, floor(12L / n_slots)),
        SlotIndex = pmin(pmax(1L, ((Col - 1L) %/% rep_block) + 1L), n_slots),
        CompoundID= paste0(SetLetter, SlotIndex)
      ) %>%
      dplyr::left_join(layout_long, by = c("CompoundID", "Row"))

    # attach Group/Patient
    if (have_sample_map) {
      df <- df %>%
        dplyr::mutate(
          ColGroup = ceiling(Col / reps_per_grp),
          Group    = group_names[ColGroup]
        )
    } else {
      # no Sample.txt: try to use sample_id from header; fallback to set label
      df <- df %>%
        dplyr::mutate(Group = ifelse(is.na(sample_id), paste0("Set", SetLetter), sample_id))
    }

    df %>%
      dplyr::mutate(
        Patient   = stringr::str_remove(Group, "(_)?cse$"),
        Condition = ifelse(stringr::str_detect(Group, "cse$|_cse$"), "CSE", "Control"),
        Experiment= basename(base)
      )
  }

  viab <- purrr::map_dfr(rows, one_block)

  # standardise concentrations
  viab <- convert_conc_to_std(viab)
  if (any(!viab$UnitOK, na.rm = TRUE)) {
    bad <- viab %>% dplyr::filter(!UnitOK) %>% dplyr::distinct(CompoundID, Unit) %>% dplyr::arrange(CompoundID)
    warning("Unrecognised units detected; proceeding with raw values:\n",
            paste0(" - ", bad$CompoundID, ": '", bad$Unit, "'", collapse = "\n"))
  }

  # robust anchors (v17): U = median at Dose==0 (if present); L = median at highest dose
  viab <- viab %>%
    dplyr::group_by(Experiment, CompoundName, Group) %>%
    dplyr::mutate(
      has_zero = any(is.finite(Signal) & is.finite(Dose) & Dose == 0),
      U = ifelse(has_zero,
                 stats::median(Signal[is.finite(Signal) & Dose == 0], na.rm = TRUE),
                 NA_real_),
      top_dose = suppressWarnings(max(Dose[is.finite(Dose)], na.rm = TRUE)),
      L = stats::median(Signal[is.finite(Signal) & Dose >= top_dose], na.rm = TRUE),
      U = ifelse(is.finite(U), U, stats::median(Signal[is.finite(Signal)], na.rm = TRUE)),
      L = ifelse(is.finite(L), L, stats::quantile(Signal[is.finite(Signal)], probs = 0.1, na.rm = TRUE)),
      Denom = pmax(U - L, 1e-6),
      Viability = pmin(pmax((Signal - L)/Denom * 100, 0), 100)
    ) %>%
    dplyr::ungroup()

  viab
}

## ---- Read all plates ---------------------------------------------------------
ensure_dir(out_dir); ensure_dir(plots_dir); ensure_dir(plots_comb_dir)

raw_files <- list.files(dir_in, pattern = "^Resazurin_Exp_\\d+_Raw\\d*\\.txt$", full.names = TRUE)
if (!length(raw_files)) stop("No 'Resazurin_Exp_*_Raw*.txt' files found in ", dir_in)
message(glue("Found {length(raw_files)} raw file(s):\n- {paste(basename(raw_files), collapse = '\n- ')}"))

all_viability <- purrr::map_dfr(raw_files, get_viability_one_plate) %>%
  dplyr::filter(is.finite(Viability), is.finite(Dose)) %>%
  dplyr::filter(tolower(CompoundName) %ni% tolower(exclude_compounds)) %>%
  dplyr::mutate(
    Sample = Group,                                  # keep your “sample” naming
    CompGroup = sub("_[^_]+$", "", CompoundName)     # base compound (before last _suffix)
  )

## ---- Robust fitting helpers (v17) -------------------------------------------
# Dose-aggregated + MAD trimming + weighted LL.4 candidates; choose by AICc
resample_fit_once <- function(df_repl){
  df_trim <- df_repl %>%
    dplyr::group_by(Dose) %>%
    dplyr::mutate(med = median(Viability, na.rm=TRUE),
                  madv = mad(Viability, constant = 1.4826, na.rm=TRUE),
                  keep = ifelse(is.finite(madv) & madv > 0, abs(Viability - med) <= 3*madv, TRUE)) %>%
    dplyr::ungroup() %>% dplyr::filter(keep)
  df_ag <- df_trim %>% dplyr::group_by(Dose) %>%
    dplyr::summarise(Viability = mean(Viability, na.rm=TRUE), wt = dplyr::n(), .groups="drop")

  if (sum(df_ag$wt) == 0 || nrow(df_ag) < cfg$min_unique_doses) return(NULL)

  fit_variant <- function(fix_c, fix_d){
    suppressWarnings(tryCatch(
      drc::drm(Viability ~ Dose, data = df_ag, weights = wt,
               fct = drc::LL.4(fixed = c(NA, fix_c, fix_d, NA)),
               control = drc::drmc(maxIt = 1000)),
      error = function(e) NULL))
  }
  cand <- list(
    free   = fit_variant(NA,   NA),
    fix_c  = fit_variant(0,    NA),
    fix_d  = fit_variant(NA,   100),
    fix_cd = fit_variant(0,    100)
  )
  keep <- vapply(cand, Negate(is.null), logical(1))
  if (!any(keep)) return(NULL)
  cand <- cand[keep]

  n <- sum(df_ag$wt)
  aicc <- function(f){ k <- length(coef(f)); den <- n - k - 1; AIC(f) + if (den>0) 2*k*(k+1)/den else Inf }
  aics <- vapply(cand, aicc, numeric(1))
  ks   <- vapply(cand, function(f) length(coef(f)), numeric(1))
  best <- which.min(aics); ties <- which(aics - aics[best] <= cfg$aicc_tie_tol)
  if (length(ties) > 1) best <- ties[ which.min(ks[ties]) ]
  best_fit <- cand[[best]]
  attr(best_fit, "variant") <- names(cand)[best]
  attr(best_fit, "data_ag") <- df_ag
  best_fit
}

safe_fit <- function(df){
  df <- df[df$Dose > 0 & is.finite(df$Dose) & is.finite(df$Viability), ]
  if (nrow(df) < cfg$min_unique_doses) return(NULL)
  resample_fit_once(df)
}

# Logistic fallback curve if model still fails
fallback_logistic <- function(df, xgrid){
  dfx <- df %>% dplyr::arrange(Dose)
  # crude IC50 via interpolation around 50
  ic50_guess <- {
    g <- dfx %>% dplyr::filter(Dose > 0, is.finite(Dose), is.finite(Viability))
    if (nrow(g) < 2) median(g$Dose, na.rm = TRUE) else {
      idx <- which(diff(sign(g$Viability - 50)) != 0)
      if (length(idx)) {
        i <- idx[1]
        x1 <- log10(g$Dose[i]); x2 <- log10(g$Dose[i+1])
        y1 <- g$Viability[i];   y2 <- g$Viability[i+1]
        10^approx(x = c(y1,y2), y = c(x1,x2), xout = 50)$y
      } else median(g$Dose, na.rm = TRUE)
    }
  }
  hill <- 1
  mu <- 100/(1 + (xgrid/ic50_guess)^hill)
  tibble::tibble(Dose = xgrid, Viability = mu, Lower = NA_real_, Upper = NA_real_, source = "fallback")
}

## ---- Predictions with bootstrap CI (percentile) ------------------------------
.ll4_mean <- function(x, b, c, d, e) c + (d - c) / (1 + (x / e)^b)

bootstrap_mean_matrix <- function(fit, df, grid, B, b_bounds = cfg$b_draw_bounds, e_bounds = NULL){
  V <- try(vcov(fit), silent = TRUE)
  if (!inherits(V, "try-error") && is.matrix(V) && all(is.finite(V))) {
    Vpd <- try(make_pd(V), silent = TRUE)
    if (!inherits(Vpd, "try-error")) {
      cf <- coef(fit)
      draws <- try(MASS::mvrnorm(B, mu = cf, Sigma = Vpd), silent = TRUE)
      if (!inherits(draws, "try-error")) {
        nb <- grep("^b:", names(cf)); ne <- grep("^e:", names(cf))
        ok <- rep(TRUE, nrow(draws))
        if (length(nb)) { b <- draws[, nb, drop=TRUE]; ok <- ok & is.finite(b) & b>b_bounds[1] & b<b_bounds[2] }
        if (length(ne)) { e <- draws[, ne, drop=TRUE]; ok <- ok & is.finite(e) & e>0
                          if (!is.null(e_bounds)) ok <- ok & e>e_bounds[1] & e<e_bounds[2] }
        draws <- draws[ok, , drop = FALSE]
        if (nrow(draws) > 0) {
          c0 <- if ("c:(Intercept)" %in% names(cf)) cf["c:(Intercept)"] else 0
          d0 <- if ("d:(Intercept)" %in% names(cf)) cf["d:(Intercept)"] else 100
          M <- sapply(seq_len(nrow(draws)), function(i){
            b <- draws[i, startsWith(names(cf), "b:")]
            e <- draws[i, startsWith(names(cf), "e:")]
            c <- if ("c:(Intercept)" %in% names(cf)) draws[i, "c:(Intercept)"] else c0
            d <- if ("d:(Intercept)" %in% names(cf)) draws[i, "d:(Intercept)"] else d0
            .ll4_mean(grid, b, c, d, e)
          })
          M <- as.matrix(M); if (ncol(M) == 1) M <- cbind(M)
          return(M)
        }
      }
    }
  }
  # fall back to nonparametric replicate bootstrap (lighter for multi-plate)
  Bnp <- min(B, 500)
  preds <- matrix(NA_real_, nrow = length(grid), ncol = Bnp)
  for (i in seq_len(Bnp)) {
    dfb <- df %>% dplyr::group_by(Dose) %>% dplyr::slice_sample(n = dplyr::n(), replace = TRUE) %>% dplyr::ungroup()
    fb  <- resample_fit_once(dfb); if (is.null(fb)) next
    preds[, i] <- as.numeric(predict(fb, newdata = data.frame(Dose = grid)))
  }
  preds
}

predict_curves <- function(viab, n = 200, level = 0.95){
  a <- (1 - level) / 2; probs <- c(a, 1 - a)

  viab %>%
    dplyr::group_by(Sample, Patient, CompoundName, UnitStd) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      fit = purrr::map(data, safe_fit),
      pred = purrr::map2(fit, data, ~{
        fit <- .x; df <- .y
        concs <- sort(unique(df$Dose[df$Dose > 0 & is.finite(df$Dose)]))
        if (!length(concs)) return(NULL)
        grid <- exp(seq(log(min(concs)), log(max(concs)), length.out = n))

        if (!is.null(fit)) {
          e_bounds <- range(concs) * cfg$e_bounds_factor
          M <- bootstrap_mean_matrix(fit, df, grid, B = cfg$B_ci,
                                     b_bounds = cfg$b_draw_bounds, e_bounds = e_bounds)
          if (!is.null(M) && ncol(M) > 0) {
            mu <- rowMeans(M, na.rm = TRUE)
            lo <- apply(M, 1, quantile, probs = probs[1], na.rm = TRUE)
            hi <- apply(M, 1, quantile, probs = probs[2], na.rm = TRUE)
            return(tibble::tibble(Dose = grid,
                                  Viability = pmin(pmax(mu, 0), 100),
                                  Lower = pmax(pmin(lo, hi), 0),
                                  Upper = pmin(pmax(lo, hi), 100),
                                  source = "model"))
          }
        }
        # fallback logistic
        fallback_logistic(df, grid)
      })
    ) %>%
    dplyr::select(-data, -fit) %>%
    tidyr::unnest(pred) %>%
    dplyr::ungroup()
}

## ---- IC50 + QC --------------------------------------------------------------
qc_metrics <- function(df, fit){
  res <- try(residuals(fit), silent = TRUE); if (inherits(res, "try-error")) res <- NA_real_
  rmse <- if (all(is.finite(res))) sqrt(mean(res^2)) else NA_real_
  sse <- if (all(is.finite(res))) sum(res^2) else NA_real_
  sst <- if (all(is.finite(df$Viability))) sum((df$Viability - mean(df$Viability, na.rm = TRUE))^2) else NA_real_
  r2  <- if (is.finite(sse) && is.finite(sst) && sst > 0) 1 - sse/sst else NA_real_
  dsum <- df %>%
    dplyr::filter(Dose > 0, is.finite(Dose), is.finite(Viability)) %>%
    dplyr::summarise(mv = mean(Viability, na.rm = TRUE), .by = Dose) %>%
    dplyr::arrange(Dose)
  rho <- if (nrow(dsum) >= 4) suppressWarnings(cor(dsum$mv, log(dsum$Dose), method = "spearman")) else NA_real_
  tibble::tibble(RMSE = rmse, R2 = r2, Spearman = rho)
}

compute_ic50 <- function(viab){
  viab %>%
    dplyr::group_by(Sample, Patient, CompoundName, UnitStd) %>%
    dplyr::group_modify(~{
      f <- safe_fit(.x)
      if (is.null(f)) return(tibble::tibble(
        IC50 = NA, IC50_L = NA, IC50_U = NA, in_range = NA,
        b = NA, c = NA, d = NA, aicc = NA, variant = NA, converged = FALSE,
        RMSE = NA, R2 = NA, Spearman = NA, n_doses = dplyr::n_distinct(.x$Dose[.x$Dose > 0])
      ))

      cf <- coef(f); nm <- names(cf)
      getp <- function(k) if (k %in% nm) unname(cf[k]) else NA_real_
      ic50 <- getp(grep("^e:", nm, value = TRUE)[1])

      ed <- quiet_ED_delta(f)
      if (!is.null(ed) && all(is.finite(ed[1, ]))) {
        ic50_l <- unname(ed[1, 2]); ic50_u <- unname(ed[1, 3])
      } else {
        ic50_l <- ic50_u <- NA_real_
        V <- try(vcov(f), silent = TRUE)
        if (!inherits(V, "try-error") && is.matrix(V) && all(is.finite(V))) {
          Vpd <- try(make_pd(V), silent = TRUE)
          if (!inherits(Vpd, "try-error")) {
            draws <- try(MASS::mvrnorm(1500, mu = cf, Sigma = Vpd), silent = TRUE)
            if (!inherits(draws, "try-error")) {
              ee <- draws[, grep("^e:", names(cf)), drop = TRUE]
              ee <- ee[is.finite(ee) & ee > 0]
              if (length(ee) >= 50) {
                qs <- quantile(ee, c(0.025, 0.975), na.rm = TRUE)
                ic50_l <- qs[[1]]; ic50_u <- qs[[2]]
              }
            }
          }
        }
      }

      concs <- sort(unique(.x$Dose[.x$Dose > 0 & is.finite(.x$Dose)]))
      aicc <- tryCatch({ k <- length(cf); n <- nrow(.x); AIC(f) + 2*k*(k+1)/max(n - k - 1, 1) }, error = function(e) NA_real_)
      converged <- tryCatch(isTRUE(f$fit$convergence == 0), error = function(e) TRUE)
      qc <- qc_metrics(.x, f)

      tibble::tibble(
        IC50 = ic50, IC50_L = ic50_l, IC50_U = ic50_u,
        in_range = is.finite(ic50) && ic50 > min(concs) && ic50 < max(concs),
        b = getp("b:(Intercept)"),
        c = getp("c:(Intercept)"),
        d = getp("d:(Intercept)"),
        aicc = aicc,
        variant = attr(f, "variant"),
        converged = converged,
        RMSE = qc$RMSE, R2 = qc$R2, Spearman = qc$Spearman,
        n_doses = length(concs)
      )
    }) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      slope_ok   = is.finite(b) & abs(b) >= cfg$slope_bounds[1] & abs(b) <= cfg$slope_bounds[2],
      mono_ok    = is.finite(Spearman) & Spearman <= -0.5,
      support_ok = is.finite(n_doses) & n_doses >= cfg$min_unique_doses
    )
}

qc_report <- function(ic50_df){
  ic50_df %>%
    dplyr::mutate(
      converged  = dplyr::coalesce(converged,  FALSE),
      support_ok = dplyr::coalesce(support_ok, FALSE),
      slope_ok   = dplyr::coalesce(slope_ok,   FALSE),
      mono_ok    = dplyr::coalesce(mono_ok,    FALSE),
      in_range   = dplyr::coalesce(in_range,   FALSE),
      ic50_ok    = is.finite(IC50),
      pass       = converged & support_ok & slope_ok & mono_ok & in_range & ic50_ok
    )
}

## ---- Fit/predict all ---------------------------------------------------------
# predictions (with CI ribbons)
all_preds_ci <- predict_curves(all_viability, level = 0.95)

# IC50 table with QC
all_ic50 <- compute_ic50(all_viability)
qc_df    <- qc_report(all_ic50)

## ---- Write CSVs --------------------------------------------------------------
ensure_dir(out_dir)
readr::write_csv(
  all_ic50 %>% dplyr::arrange(CompoundName, Sample, Patient),
  file.path(out_dir, "IC50_results.csv")
)
readr::write_csv(
  qc_df %>% dplyr::arrange(CompoundName, Sample, Patient),
  file.path(out_dir, "fit_diagnostics.csv")
)
# simple “suggested combos” ranking
suggested <- qc_df %>%
  dplyr::mutate(score = 0.5*as.numeric(in_range) + 0.2*pmax(0,R2) + 0.3*pmax(0, -Spearman)) %>%
  dplyr::arrange(dplyr::desc(score))
readr::write_csv(suggested, file.path(out_dir, "suggested_combinations.csv"))

## ---- Plot helpers ------------------------------------------------------------
unit_label_for <- function(df_units){
  u <- unique(df_units)
  if (length(u) == 1) sprintf("Concentration (%s)", u) else "Concentration [per panel]"
}

plot_single <- function(sample, compound){
  dline <- dplyr::filter(all_preds_ci, Sample == sample, CompoundName == compound)
  dpts  <- dplyr::filter(all_viability, Sample == sample, CompoundName == compound, Dose > 0, is.finite(Dose))
  if (!nrow(dpts)) return(NULL)
  xlab <- unit_label_for(dline$UnitStd)

  ggplot() +
    { if (any(is.finite(dline$Lower) & is.finite(dline$Upper)))
        geom_ribbon(data = dline, aes(Dose, ymin = Lower, ymax = Upper), alpha = .15, colour = NA) else NULL } +
    geom_point(data = dpts, aes(Dose, Viability), size = 1.7, alpha = .8) +
    geom_line (data = dline, aes(Dose, Viability), linewidth = .8) +
    scale_x_log10(labels = label_number()) +
    coord_cartesian(ylim = c(0, 100)) +
    labs(title = glue("{compound} — {sample}"), x = xlab, y = "Viability (%)") +
    base_theme
}

# across samples: base compound + co-treatments (facet by Sample)
plot_compound_across_samples <- function(comp_group){
  d <- all_viability %>% dplyr::filter(CompGroup == comp_group, Dose > 0)
  if (!nrow(d)) return(NULL)
  dline <- all_preds_ci %>% dplyr::filter(CompoundName %in% unique(d$CompoundName))
  xlab  <- unit_label_for(dline$UnitStd)

  ggplot() +
    { if (nrow(dline) && any(is.finite(dline$Lower) & is.finite(dline$Upper)))
        geom_ribbon(data = dline, aes(Dose, ymin = Lower, ymax = Upper, fill = CompoundName),
                    alpha = .12, colour = NA, show.legend = FALSE) else NULL } +
    geom_point(data = d, aes(Dose, Viability, colour = CompoundName), size = 1.6, alpha = .75) +
    geom_line (data = dline, aes(Dose, Viability, colour = CompoundName), linewidth = .8) +
    scale_x_log10(labels = label_number()) +
    coord_cartesian(ylim = c(0, 100)) +
    facet_wrap(~ Sample, scales = "free_x") +
    labs(title = glue("{comp_group} and co-treatments — across samples"),
         x = xlab, y = "Viability (%)", colour = "Treatment") +
    base_theme + theme(legend.position = "bottom")
}

# within a sample: base vs co-treatments for that sample
plot_sample_vs_cotreats <- function(sample, comp_group){
  d <- all_viability %>% dplyr::filter(Sample == sample, CompGroup == comp_group, Dose > 0)
  if (!nrow(d)) return(NULL)
  dline <- all_preds_ci %>% dplyr::filter(Sample == sample, CompoundName %in% unique(d$CompoundName))
  xlab  <- unit_label_for(dline$UnitStd)

  ggplot() +
    { if (nrow(dline) && any(is.finite(dline$Lower) & is.finite(dline$Upper)))
        geom_ribbon(data = dline, aes(Dose, ymin = Lower, ymax = Upper, fill = CompoundName),
                    alpha = .12, colour = NA, show.legend = FALSE) else NULL } +
    geom_point(data = d, aes(Dose, Viability, colour = CompoundName), size = 1.7, alpha = .8) +
    geom_line (data = dline, aes(Dose, Viability, colour = CompoundName), linewidth = .8) +
    scale_x_log10(labels = label_number()) +
    coord_cartesian(ylim = c(0, 100)) +
    labs(title = glue("{comp_group} — {sample} (base vs co-treatments)"),
         x = xlab, y = "Viability (%)", colour = "Treatment") +
    base_theme + theme(legend.position = "bottom")
}

## ---- Save plots --------------------------------------------------------------
# singles
all_viability %>% dplyr::distinct(Sample, CompoundName) %>%
  purrr::pwalk(function(Sample, CompoundName){
    gp <- plot_single(Sample, CompoundName)
    save_plot(gp, file.path(plots_dir, glue("{CompoundName}__{Sample}")))
  })

# across-sample combined (group base + *_suffix)
for (grp in sort(unique(all_viability$CompGroup))) {
  gp <- plot_compound_across_samples(grp)
  save_plot(gp, file.path(plots_comb_dir, glue("{grp}__across_samples")))
}

# within-sample combined
for (smp in sort(unique(all_viability$Sample))) {
  for (grp in sort(unique(all_viability$CompGroup))) {
    gp <- plot_sample_vs_cotreats(smp, grp)
    if (!is.null(gp)) save_plot(gp, file.path(plots_comb_dir, glue("{grp}__{smp}")))
  }
}

message("Outputs written to:\n",
        "- ", file.path(out_dir, "IC50_results.csv"), "\n",
        "- ", file.path(out_dir, "fit_diagnostics.csv"), "\n",
        "- ", file.path(out_dir, "suggested_combinations.csv"), "\n",
        "- Single plots: ", normalizePath(plots_dir, winslash = "/"), "\n",
        "- Combined plots: ", normalizePath(plots_comb_dir, winslash = "/"))
