##############################################################################
#   IC50 analysis — v17 (publication-grade robustness)          2025-08-26
##############################################################################

## ---- Paths ------------------------------------------------------------------
raw_path      <- "IC50/Resazurin_Exp_36/20250819_Resazurin ASC397 ASC397cse ASC399 ASC399cse CSE.txt"
compound_path <- "IC50/Resazurin_Exp_36/Resazurin_Exp_36_Compound.txt"
sample_path   <- "IC50/Resazurin_Exp_36/Resazurin_Exp_36_Sample.txt"

## ---- Config -----------------------------------------------------------------
cfg <- list(
  seed               = 1L,
  min_unique_doses   = 4,            # require ≥4 unique positive doses for a fit
  slope_bounds       = c(0.2, 20),   # acceptable |b| range for QC
  aicc_tie_tol       = 2,            # ΔAICc ≤ 2 => prefer simpler model
  B_ci               = 3000,         # bootstrap draws for CI (percentile)
  B_pi               = 3000,         # bootstrap draws for PI (percentile)
  b_draw_bounds      = c(0.1, 50),   # sanity for bootstrap draws
  e_bounds_factor    = c(0.25, 4),   # e must lie within [0.25*min, 4*max] doses during draws
  loess_span         = 0.8,          # smoothness for σ(x) loess
  unit_default       = "µM"          # unified unit for plotting/fitting
)

set.seed(cfg$seed)

## ---- Libraries ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(MASS)        # put MASS first
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(purrr); library(tibble); library(drc); library(ggplot2)
})

## ---- PD repair + quiet ED + dose resample -----------------------------------
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

# Nonparametric bootstrap: resample replicates within each dose, refit, predict/ED
resample_fit_once <- function(df_repl){
  # MAD-based outlier trimming, then aggregate to dose means with weights
  df_trim <- df_repl %>%
    group_by(Dose) %>%
    mutate(med = median(Viability, na.rm=TRUE),
           madv = mad(Viability, constant = 1.4826, na.rm=TRUE),
           keep = ifelse(is.finite(madv) & madv > 0, abs(Viability - med) <= 3*madv, TRUE)) %>%
    ungroup() %>% filter(keep)
  df_ag <- df_trim %>% group_by(Dose) %>%
    summarise(Viability = mean(Viability, na.rm=TRUE), wt = n(), .groups="drop")
  if (nrow(df_ag) < cfg$min_unique_doses) return(NULL)

  fit_variant <- function(fix_c, fix_d){
    suppressWarnings(tryCatch(
      drm(Viability ~ Dose, data = df_ag, weights = wt,
          fct = LL.4(fixed = c(NA, fix_c, fix_d, NA)),
          control = drmc(maxIt = 1000)),
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
  attr(best_fit, "data_ag") <- df_ag      # keep aggregated data for later
  best_fit
}

## ---- Concentration standardiser (multi-unit) ----
.standardise_conc <- function(x, unit){
  u <- tolower(gsub("\\s+", "", unit %||% ""))
  u <- gsub("μ","u",u)  # normalise micro symbol

  res <- list(value = x, unit_std = unit, ok = TRUE, convertible = FALSE)

  # molarity → µM
  if (u %in% c("m","mol/l"))              { res$value <- x*1e6; res$unit_std <- "µM"; res$convertible <- TRUE; return(res) }
  if (u %in% c("mm","mmol/l"))            { res$value <- x*1e3; res$unit_std <- "µM"; res$convertible <- TRUE; return(res) }
  if (u %in% c("um","µm","umol/l","µmol/l")) { res$value <- x;     res$unit_std <- "µM"; res$convertible <- TRUE; return(res) }
  if (u %in% c("nm","nmol/l"))            { res$value <- x*1e-3; res$unit_std <- "µM"; res$convertible <- TRUE; return(res) }
  if (u %in% c("pm","pmol/l"))            { res$value <- x*1e-6; res$unit_std <- "µM"; res$convertible <- TRUE; return(res) }

  # percentage
  if (u %in% c("%","percent","%v/v","v/v%","%w/v","w/v%")) { res$value <- x; res$unit_std <- "%";     return(res) }

  # common mass/volume units (left as-is; cannot convert to µM without MW)
  if (u %in% c("mg/ml","mgperml"))        { res$value <- x; res$unit_std <- "mg/mL"; return(res) }
  if (u %in% c("ug/ml","µg/ml","mcg/ml","ugperml","µgperml")) { res$value <- x; res$unit_std <- "µg/mL"; return(res) }
  if (u %in% c("ng/ml","ngperml"))        { res$value <- x; res$unit_std <- "ng/mL"; return(res) }

  # fallback: unknown unit — keep numeric value, mark not-ok
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

## ---- Layout readers ----------------------------------------------------------
build_layout <- function(path){
  raw <- read_delim(path, "\t", col_names = FALSE, trim_ws = TRUE, show_col_types = FALSE)
  ids   <- as.character(raw[raw$X1 == "ID",       -1])
  names <- as.character(raw[raw$X1 == "Compound", -1])
  units <- as.character(raw[raw$X1 == "Unit",     -1])
  dose_rows <- raw %>% filter(str_detect(X1, "^[0-9]+$"))
  map_dfr(seq_along(ids), function(j){
    map_dfr(seq_len(nrow(dose_rows)), function(i){
      tibble(CompoundID=ids[j], CompoundName=names[j],
             Row=LETTERS[as.integer(dose_rows$X1[i])],
             Concentration=as.numeric(dose_rows[i, j+1][[1]]), Unit=units[j])
    })
  })
}

get_group_names <- function(path){
  sraw <- read_delim(path, "\t", col_names = FALSE, trim_ws = TRUE, show_col_types = FALSE)
  as.character(sraw[sraw$X1 == "Compound", -1])
}

## ---- Use robust aggregated fit for stability --------------------------------
safe_fit <- function(df){
  df <- df[is.finite(df$Dose) & is.finite(df$Viability) & df$Dose > 0, ]
  if (nrow(df) < cfg$min_unique_doses) return(NULL)
  # bootstrap-compatible single fit
  resample_fit_once(df)
}

## ---- Extract + normalise viability (robust anchors) --------------------------
get_viability_data <- function(rows, all_lines, layout, group_names, reps_per_group){
  viab <- map_dfr(rows, function(i){
    plate <- suppressMessages(read_delim(I(paste(all_lines[i:(i+8)], collapse = "\n")),
                                         "\t", trim_ws = TRUE, show_col_types = FALSE))
    fluo_start <- which(str_detect(names(plate), regex("_fluo$", ignore_case = TRUE)))[1]
    sig <- plate[, (fluo_start + 1):(fluo_start + 12)] |> setNames(sprintf("%02d", 1:12))
    PlateID <- str_remove(names(plate)[fluo_start], "_[Ff]luo$")

    sig %>%
      mutate(Row = LETTERS[1:8]) %>%
      pivot_longer(-Row, names_to = "Col", values_to = "Signal",
                   names_transform = list(Col = as.integer)) %>%
      mutate(
        PlateID    = PlateID,
        ColGroup   = ceiling(Col / reps_per_group),   # 1..(12/rep)
        Group      = group_names[ColGroup],
        Patient    = ifelse(str_detect(Group, "_cse$"), str_remove(Group, "_cse$"),
                            ifelse(str_detect(Group, "cse$"), str_remove(Group, "cse$"), Group)),
        Condition  = ifelse(str_detect(Group, "cse$|_cse$"), "CSE", "Control"),
        CompoundID = PlateID
      ) %>%
      left_join(layout, by = c("CompoundID", "Row"))
  })

  # Standardise concentrations (supports µM, %, mg/mL, etc.)
  viab <- convert_conc_to_std(viab)
  if (any(!viab$UnitOK, na.rm = TRUE)) {
    bad <- viab %>% filter(!UnitOK) %>% distinct(CompoundID, Unit) %>% arrange(CompoundID)
    warning("Unrecognised units detected; proceeding with raw values:\n",
            paste0(" - ", bad$CompoundID, ": '", bad$Unit, "'", collapse = "\n"))
  }

  # Robust per-plate × group normalisation (U from 0-dose; L from top dose)
  viab <- viab %>%
    group_by(PlateID, CompoundID, Group) %>%
    mutate(
      has_zero = any(is.finite(Signal) & is.finite(Dose) & Dose == 0),
      U = ifelse(has_zero,
                 median(Signal[is.finite(Signal) & Dose == 0], na.rm = TRUE),
                 NA_real_),
      top_dose = suppressWarnings(max(Dose[is.finite(Dose)], na.rm = TRUE)),
      L = median(Signal[is.finite(Signal) & Dose >= top_dose], na.rm = TRUE),
      # Fallbacks if anchors are NA/degenerate
      U = ifelse(is.finite(U), U, median(Signal[is.finite(Signal)], na.rm = TRUE)),
      L = ifelse(is.finite(L), L, quantile(Signal[is.finite(Signal)], probs = 0.1, na.rm = TRUE)),
      Denom = pmax(U - L, 1e-6),
      Viability = pmin(pmax((Signal - L) / Denom * 100, 0), 100)
    ) %>%
    ungroup()

  viab
}

## ---- Replicate validation ----------------------------------------------------
validate_replicates <- function(viab){
  rep_chk <- viab %>%
    filter(is.finite(Dose)) %>%
    count(PlateID, Group, Dose, name = "n_rep") %>%
    summarise(min_rep = min(n_rep), max_rep = max(n_rep), .by = c(PlateID, Group)) %>%
    mutate(ok = min_rep == max_rep & min_rep %in% c(2, 3, 4))
  if (any(!rep_chk$ok)) {
    message("⚠ Replicate inconsistency detected:")
    print(rep_chk %>% filter(!ok))
  }
  invisible(rep_chk)
}

## ---- QC metrics --------------------------------------------------------------
qc_metrics <- function(df, fit){
  # Basic residual metrics
  res <- try(residuals(fit), silent = TRUE); if (inherits(res, "try-error")) res <- NA_real_
  rmse <- if (all(is.finite(res))) sqrt(mean(res^2)) else NA_real_

  # Pseudo-R2 (vs. mean of Viability)
  sse <- if (all(is.finite(res))) sum(res^2) else NA_real_
  sst <- if (all(is.finite(df$Viability))) sum((df$Viability - mean(df$Viability, na.rm = TRUE))^2) else NA_real_
  r2  <- if (is.finite(sse) && is.finite(sst) && sst > 0) 1 - sse/sst else NA_real_

  # Monotonicity: Spearman rho of mean(Y) vs log-dose
  dsum <- df %>%
    filter(Dose > 0, is.finite(Dose), is.finite(Viability)) %>%
    summarise(mv = mean(Viability, na.rm = TRUE), .by = Dose) %>%
    arrange(Dose)
  rho <- if (nrow(dsum) >= 4) suppressWarnings(cor(dsum$mv, log(dsum$Dose), method = "spearman")) else NA_real_

  tibble(RMSE = rmse, R2 = r2, Spearman = rho)
}

## ---- Bootstrap utilities for bands -------------------------------------------
.ll4_mean <- function(x, b, c, d, e) c + (d - c) / (1 + (x / e)^b)

bootstrap_mean_matrix <- function(fit, df, grid, B, b_bounds = cfg$b_draw_bounds, e_bounds = NULL){
  # Try parameter bootstrap from vcov; repair Sigma if needed
  V <- try(vcov(fit), silent = TRUE)
  if (!inherits(V, "try-error") && is.matrix(V) && all(is.finite(V))) {
    Vpd <- try(make_pd(V), silent = TRUE)
    if (!inherits(Vpd, "try-error")) {
      cf <- coef(fit)
      draws <- try(MASS::mvrnorm(B, mu = cf, Sigma = Vpd), silent = TRUE)
      if (!inherits(draws, "try-error")) {
        nm <- colnames(draws); has_c <- "c:(Intercept)" %in% names(cf); has_d <- "d:(Intercept)" %in% names(cf)
        nb <- grep("^b:", names(cf)); ne <- grep("^e:", names(cf))
        ok <- rep(TRUE, nrow(draws))
        if (length(nb)) { b <- draws[, nb, drop=TRUE]; ok <- ok & is.finite(b) & b>b_bounds[1] & b<b_bounds[2] }
        if (length(ne)) { e <- draws[, ne, drop=TRUE]; ok <- ok & is.finite(e) & e>0
                          if (!is.null(e_bounds)) ok <- ok & e>e_bounds[1] & e<e_bounds[2] }
        draws <- draws[ok, , drop = FALSE]
        if (nrow(draws) > 0) {
          c0 <- if (has_c) cf["c:(Intercept)"] else 0
          d0 <- if (has_d) cf["d:(Intercept)"] else 100
          M <- sapply(seq_len(nrow(draws)), function(i){
            b <- draws[i, startsWith(names(cf), "b:")]
            e <- draws[i, startsWith(names(cf), "e:")]
            c <- if (has_c) draws[i, "c:(Intercept)"] else c0
            d <- if (has_d) draws[i, "d:(Intercept)"] else d0
            .ll4_mean(grid, b, c, d, e)
          })
          M <- as.matrix(M); if (ncol(M) == 1) M <- cbind(M)
          return(M)
        }
      }
    }
  }

  # Fallback: nonparametric bootstrap by resampling replicates within dose
  Bnp <- min(B, 800)  # keep it reasonable if we fell back
  preds <- matrix(NA_real_, nrow = length(grid), ncol = Bnp)
  # Use original replicate-level df (group-wise nest passed this df)
  for (i in seq_len(Bnp)) {
    dfb <- df %>% group_by(Dose) %>% slice_sample(n = n(), replace = TRUE) %>% ungroup()
    fb  <- resample_fit_once(dfb); if (is.null(fb)) next
    preds[, i] <- as.numeric(predict(fb, newdata = data.frame(Dose = grid)))
  }
  preds
}

sigma_x_from_loess <- function(fit, grid, span = cfg$loess_span){
  res <- residuals(fit); sigma <- sd(res, na.rm = TRUE); if (!is.finite(sigma)) sigma <- 0
  if (length(res) >= 6) {
    tmp <- data.frame(fitv = fitted(fit), absr = abs(res))
    lo  <- stats::loess(absr ~ fitv, data = tmp, span = span, degree = 1)
    sig_x <- predict(lo, newdata = data.frame(fitv = .ll4_mean(grid,
                                                               coef(fit)[startsWith(names(coef(fit)), "b:")],
                                                               ifelse("c:(Intercept)" %in% names(coef(fit)), coef(fit)["c:(Intercept)"], 0),
                                                               ifelse("d:(Intercept)" %in% names(coef(fit)), coef(fit)["d:(Intercept)"], 100),
                                                               coef(fit)[startsWith(names(coef(fit)), "e:")])))
    # Extrapolate if needed
    if (anyNA(sig_x)) {
      xs <- sort(unique(tmp$fitv)); ys <- predict(lo, newdata = data.frame(fitv = xs))
      sig_x <- approx(xs, ys, xout = .ll4_mean(grid,
                                               coef(fit)[startsWith(names(coef(fit)), "b:")],
                                               ifelse("c:(Intercept)" %in% names(coef(fit)), coef(fit)["c:(Intercept)"], 0),
                                               ifelse("d:(Intercept)" %in% names(coef(fit)), coef(fit)["d:(Intercept)"], 100),
                                               coef(fit)[startsWith(names(coef(fit)), "e:")]),
                      rule = 2)$y
    }
    sig_x <- pmax(sig_x, 1e-6)
  } else sig_x <- rep(sigma, length(grid))
  sig_x
}

## ---- Predictions with percentile bands ---------------------------------------
predict_curves <- function(viab, n = 200, band = c("ci", "pi"), level = 0.95){
  band <- match.arg(band); a <- (1 - level) / 2; probs <- c(a, 1 - a)

  viab %>%
    group_by(Patient, Group, CompoundID, CompoundName, Unit, UnitStd) %>%
    nest() %>%
    mutate(
      fit = map(data, safe_fit),
      pred = map2(fit, data, ~{
        fit <- .x; df <- .y
        if (is.null(fit)) return(NULL)

        concs <- sort(unique(df$Dose[df$Dose > 0 & is.finite(df$Dose)]))
        if (length(concs) < cfg$min_unique_doses) return(NULL)
        grid <- exp(seq(log(min(concs)), log(max(concs)), length.out = n))

        # Bootstrap (parameter with PD repair, else nonparam)
        e_bounds <- range(concs) * cfg$e_bounds_factor
        M <- bootstrap_mean_matrix(fit, df, grid,
                                   B = if (band == "ci") cfg$B_ci else cfg$B_pi,
                                   b_bounds = cfg$b_draw_bounds, e_bounds = e_bounds)
        if (is.null(M) || !ncol(M)) {
          mu <- as.numeric(predict(fit, newdata = data.frame(Dose = grid)))
          return(tibble(Dose = grid, Viability = pmin(pmax(mu, 0), 100),
                        Lower = pmin(pmax(mu, 0), 100), Upper = pmin(pmax(mu, 0), 100)))
        }

        mu <- rowMeans(M, na.rm = TRUE)

        if (band == "ci") {
          lo <- apply(M, 1, quantile, probs = probs[1], na.rm = TRUE)
          hi <- apply(M, 1, quantile, probs = probs[2], na.rm = TRUE)
        } else {
          # Heteroscedastic noise via LOESS (use aggregated data if available)
          sig_x <- sigma_x_from_loess(fit, grid)
          MP    <- M + replicate(ncol(M), rnorm(length(grid), 0, sig_x))
          lo <- apply(MP, 1, quantile, probs = probs[1], na.rm = TRUE)
          hi <- apply(MP, 1, quantile, probs = probs[2], na.rm = TRUE)
        }

        tibble(Dose = grid,
               Viability = pmin(pmax(mu, 0), 100),
               Lower     = pmax(pmin(lo, hi), 0),
               Upper     = pmin(pmax(lo, hi), 100))
      })
    ) %>%
    select(-data, -fit) %>%
    unnest(pred)
}

## ---- IC50 + CIs + flags ------------------------------------------------------
compute_ic50 <- function(viab){
  viab %>%
    group_by(Patient, Group, CompoundID, CompoundName, Unit, UnitStd) %>%
    group_modify(~{
      f <- safe_fit(.x)
      if (is.null(f)) return(tibble(
        IC50 = NA, IC50_L = NA, IC50_U = NA, in_range = NA,
        b = NA, c = NA, d = NA, aicc = NA, variant = NA, converged = FALSE,
        RMSE = NA, R2 = NA, Spearman = NA, n_doses = n_distinct(.x$Dose[.x$Dose > 0])
      ))

      cf <- coef(f); nm <- names(cf)
      getp <- function(k) if (k %in% nm) unname(cf[k]) else NA_real_

      # Point estimate (parameter 'e')
      ic50 <- getp(grep("^e:", nm, value = TRUE)[1])

      # Try delta-method CI quietly
      ed <- quiet_ED_delta(f)
      if (!is.null(ed) && all(is.finite(ed[1, ]))) {
        ic50_l <- unname(ed[1, 2]); ic50_u <- unname(ed[1, 3])
      } else {
        # Bootstrap CI: parameter-draw if vcov ok, else nonparam resampling
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
        if (!is.finite(ic50_l) || !is.finite(ic50_u)) {
          # Nonparametric bootstrap on replicates
          Bnp <- 800
          ee  <- rep(NA_real_, Bnp)
          for (i in seq_len(Bnp)) {
            dfb <- .x %>% group_by(Dose) %>% slice_sample(n = n(), replace = TRUE) %>% ungroup()
            fb  <- resample_fit_once(dfb); if (is.null(fb)) next
            cf2 <- coef(fb); nm2 <- names(cf2)
            ei  <- if ("e:(Intercept)" %in% nm2) unname(cf2["e:(Intercept)"]) else NA_real_
            ee[i] <- ei
          }
          ee <- ee[is.finite(ee) & ee > 0]
          if (length(ee) >= 50) {
            qs <- quantile(ee, c(0.025, 0.975), na.rm = TRUE)
            ic50_l <- qs[[1]]; ic50_u <- qs[[2]]
          }
        }
      }

      concs <- sort(unique(.x$Dose[.x$Dose > 0 & is.finite(.x$Dose)]))
      aicc <- tryCatch({ k <- length(cf); n <- nrow(.x); AIC(f) + 2*k*(k+1)/max(n - k - 1, 1) }, error = function(e) NA_real_)
      converged <- tryCatch(isTRUE(f$fit$convergence == 0), error = function(e) TRUE)  # consider TRUE if object exists

      qc <- qc_metrics(.x, f)

      tibble(
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
    ungroup() %>%
    mutate(
      slope_ok   = is.finite(b) & abs(b) >= cfg$slope_bounds[1] & abs(b) <= cfg$slope_bounds[2],
      mono_ok    = is.finite(Spearman) & Spearman <= -0.5,
      support_ok = is.finite(n_doses) & n_doses >= cfg$min_unique_doses
    )
}

## ---- QC report utilities -----------------------------------------------------
qc_report <- function(ic50_df){
  ic50_df %>%
    mutate(
      converged  = coalesce(converged,  FALSE),
      support_ok = coalesce(support_ok, FALSE),
      slope_ok   = coalesce(slope_ok,   FALSE),
      mono_ok    = coalesce(mono_ok,    FALSE),
      in_range   = coalesce(in_range,   FALSE),
      ic50_ok    = is.finite(IC50),
      pass       = converged & support_ok & slope_ok & mono_ok & in_range & ic50_ok,
      flag_reasons = pmap_chr(
        list(converged, support_ok, slope_ok, mono_ok, in_range, ic50_ok),
        function(con, sup, sl, mono, inr, icok){
          bad <- c(
            if (!isTRUE(con))  "non-convergence",
            if (!isTRUE(sup))  "<4 unique doses",
            if (!isTRUE(sl))   "slope |b| out of bounds",
            if (!isTRUE(mono)) "non-monotonic",
            if (!isTRUE(inr))  "ED50 out of range",
            if (!isTRUE(icok)) "IC50 NA"
          )
          if (length(bad) == 0) "OK" else paste(bad, collapse = "; ")
        }
      )
    ) %>%
    arrange(Patient, Group, CompoundName, CompoundID)
}

qc_summary <- function(ic50_df){
  df <- qc_report(ic50_df)
  cat("\n================  IC50 QC SUMMARY  ================\n")
  summ <- df %>% summarise(
    n_total = n(),
    n_pass  = sum(pass, na.rm = TRUE),
    n_fail  = sum(!pass, na.rm = TRUE)
  )
  print(summ)

  if (any(!df$pass)) {
    cat("\nTop failure reasons:\n")
    print(df %>% filter(!pass) %>% count(flag_reasons, sort = TRUE) %>% head(10))
    cat("\nFirst few failed curves:\n")
    print(df %>% filter(!pass) %>%
            select(Patient, Group, CompoundName, CompoundID, IC50, IC50_L, IC50_U, flag_reasons) %>%
            head(12))
  } else {
    cat("\nAll curves passed QC.\n")
  }
  invisible(df)
}

## ---- Plot helpers ------------------------------------------------------------
.plot_base <- function(){
  theme_minimal(base_size = 11) %+replace% theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
}

plot_across <- function(preds, viab, cmpd){
  df_line <- filter(preds, CompoundID == cmpd)
  df_pts  <- filter(viab,  CompoundID == cmpd)
  if (!nrow(df_line)) stop("No predictions for CompoundID = ", cmpd)

  u <- unique(df_line$UnitStd)
  xlab <- if (length(u) == 1) sprintf("Dose [%s]", u) else "Dose [mixed units]"
  df_pts <- df_pts %>% filter(is.finite(Dose) & Dose > 0)

  ggplot() +
    geom_ribbon(data = df_line,
                aes(Dose, ymin = Lower, ymax = Upper,
                    fill = interaction(Patient, Group)),
                alpha = .20, colour = NA, show.legend = FALSE) +
    geom_point(data = df_pts,
               aes(Dose, Viability, colour = interaction(Patient, Group)),
               size = 1.6, alpha = .7) +
    geom_line(data = df_line,
              aes(Dose, Viability, colour = interaction(Patient, Group)),
              linewidth = .8) +
    scale_x_log10() +
    coord_cartesian(ylim = c(0, 100)) +
    .plot_base() +
    labs(title  = sprintf("%s across patients & conditions", unique(df_line$CompoundName)),
         x      = xlab,
         y      = "Viability (%)",
         colour = "Patient / Group")
}

plot_within <- function(preds, viab, cmpd){
  df_line <- filter(preds, CompoundID == cmpd)
  df_pts  <- filter(viab,  CompoundID == cmpd)
  if (!nrow(df_line)) stop("No predictions for CompoundID = ", cmpd)

  u <- unique(df_line$UnitStd)
  xlab <- if (length(u) == 1) sprintf("Dose [%s]", u) else "Dose [mixed units]"
  df_pts <- df_pts %>% filter(is.finite(Dose) & Dose > 0)

  ggplot() +
    geom_ribbon(data = df_line,
                aes(Dose, ymin = Lower, ymax = Upper, fill = Group),
                alpha = .20, colour = NA, show.legend = FALSE) +
    geom_point(data = df_pts,
               aes(Dose, Viability, colour = Group),
               size = 1.6, alpha = .7) +
    geom_line(data = df_line,
              aes(Dose, Viability, colour = Group),
              linewidth = .8) +
    scale_x_log10() +
    coord_cartesian(ylim = c(0, 100)) +
    facet_wrap(~ Patient) +
    .plot_base() +
    labs(title = sprintf("%s by patient conditions", unique(df_line$CompoundName)),
         x = xlab, 
         y = "Viability (%)")
}

# Compare compounds within a single patient (optionally per group)
plot_patient_compounds <- function(preds, viab, patient, group = NULL, facet_by_unit = TRUE) {
  df_line <- dplyr::filter(preds, Patient == patient)
  df_pts  <- dplyr::filter(viab,  Patient == patient)

  if (!is.null(group)) {
    df_line <- dplyr::filter(df_line, Group == group)
    df_pts  <- dplyr::filter(df_pts,  Group == group)
  }
  if (nrow(df_line) == 0 || nrow(df_pts) == 0) {
    stop(sprintf("No data found for patient '%s'%s.",
                 patient, if (!is.null(group)) paste0(" and group '", group, "'") else ""))
  }

  df_line <- dplyr::mutate(df_line, CompoundLabel = paste0(CompoundName, " (", CompoundID, ")"))
  df_pts  <- dplyr::mutate(df_pts,  CompoundLabel = paste0(CompoundName, " (", CompoundID, ")"))
  df_pts  <- df_pts %>% dplyr::filter(is.finite(Dose) & Dose > 0)

  units_present <- unique(df_line$UnitStd)
  xlab <- if (length(units_present) == 1) sprintf("Dose [%s]", units_present) else "Dose [per unit panel]"

  p <- ggplot() +
    geom_ribbon(data = df_line,
                aes(Dose, ymin = Lower, ymax = Upper, fill = CompoundLabel),
                alpha = .15, colour = NA, show.legend = FALSE) +
    geom_point(data = df_pts,
               aes(Dose, Viability, colour = CompoundLabel),
               size = 1.6, alpha = .7) +
    geom_line(data = df_line,
              aes(Dose, Viability, colour = CompoundLabel),
              linewidth = .8) +
    scale_x_log10() +
    coord_cartesian(ylim = c(0, 100)) +
    .plot_base() +
    labs(
      title  = sprintf("Patient %s — dose–response across compounds%s",
                       patient, if (is.null(group)) "" else paste0(" (", group, ")")),
      x      = xlab, y = "Viability (%)",
      colour = "Compound"
    )

  if (facet_by_unit && length(units_present) > 1) {
    if (is.null(group)) {
      p <- p + facet_grid(UnitStd ~ Group, scales = "free_x")
    } else {
      p <- p + facet_wrap(~ UnitStd, nrow = 1, scales = "free_x")
    }
  } else if (is.null(group)) {
    p <- p + facet_wrap(~ Group, nrow = 1)
  }

  p
}

## ---- Pipeline ----------------------------------------------------------------
# Read inputs
layout_long  <- build_layout(compound_path)
group_names  <- get_group_names(sample_path)
reps_per_grp <- 12 / length(group_names)

all_lines <- read_lines(raw_path)
rows <- which(str_detect(all_lines, regex("_fluo", ignore_case = TRUE)))
rows <- rows[rows + 8 <= length(all_lines)]
if (length(rows) == 0) stop("No \"_Fluo\" blocks found.")

# Viability with robust normalisation + µM harmonisation
all_viability <- get_viability_data(rows, all_lines, layout_long, group_names, reps_per_grp)

# Safety: each compound must have a single UnitStd (no mixing within a compound)
per_cmpd <- all_viability %>% distinct(CompoundID, UnitStd) %>% count(CompoundID, name = "n_units")
if (any(per_cmpd$n_units > 1)) {
  offenders <- per_cmpd %>% filter(n_units > 1) %>% pull(CompoundID)
  stop("Mixed units found within the same compound: ", paste(offenders, collapse = ", "),
       ". Please fix the input layout/units for these plates.")
}

# Replicate validation
validate_replicates(all_viability)

# IC50 table with CIs and QC flags
all_ic50 <- compute_ic50(all_viability)

# QC report (prints summary; returns full annotated table)
qc_df <- qc_summary(all_ic50)

# Optional: persist full QC table for your records / supplement
# write_csv(qc_df, "IC50_QC_report.csv")

# Predictions (percentile bands)
all_preds_ci <- predict_curves(all_viability, band = "ci", level = .95)
all_preds_pi <- predict_curves(all_viability, band = "pi", level = .95)

## ---- Methods blurb -----------------------------------------------------------
methods_blurb <- function(cfg){
  N  <- cfg$min_unique_doses
  bL <- cfg$slope_bounds[1]; bU <- cfg$slope_bounds[2]
  B  <- cfg$B_ci; span <- cfg$loess_span
  paste(
"Dose–response modelling and IC50 estimation.",
"",
"Signals were normalised per plate and group using a two-point anchor. ",
"The upper anchor (U) was the median signal at 0-dose; the lower anchor (L) ",
"was the median signal at the highest tested dose. Viability was calculated as ",
"100 × (S − L)/(U − L) and truncated to [0, 100]. Robust fallbacks (overall medians ",
"or quantiles) were used if U or L were unavailable.",
"",
"Concentrations were converted to µM prior to modelling. For each patient × group × ",
"compound, a four-parameter log-logistic model (LL.4) was fitted to positive doses ",
"(log-scaled x). The selected variant (free; c fixed to 0; d fixed to 100; or both) ",
"was chosen by AICc, preferring the simpler model when ΔAICc ≤ 2. Fits required at ",
least(N, " unique positive doses and used increased iteration limits (drc::drm)."),
"",
"Uncertainty in the mean curve was quantified via parametric bootstrap draws ",
"from the asymptotic parameter covariance (", B, " draws). Pointwise 95% confidence ",
"bands were the 2.5th–97.5th percentiles of the bootstrap distribution. Predictive ",
"95% bands incorporated heteroscedastic residual variance modelled by LOESS of ",
"|residuals| versus fitted values (span = ", sprintf("%.2f", span),
"), with noise added inside each bootstrap path before computing percentiles.",
"",
"The IC50 (ED50) was extracted from LL.4 with 95% confidence intervals using the ",
"delta method. An in-range flag indicates whether the ED50 lies within the tested ",
"dose range.",
"",
"Quality control comprised replicate-count checks, convergence status, dose support ",
"(≥ ", N, " unique positive doses), slope bounds (|b| ∈ [", sprintf("%.1f", bL), ", ",
sprintf("%.1f", bU), "]), pseudo-R² and RMSE, and monotonicity (Spearman correlation ",
"between mean viability and log dose). Non-conforming curves were flagged and excluded ",
"from summary statements unless explicitly stated.",
"",
"Analyses were conducted in R with packages drc and MASS; code and configuration ",
"are archived with session information for reproducibility.",
sep = ""
  )
}

print(plot_patient_compounds(all_preds_ci, all_viability, "ASC397", group = "ASC397"))
print(plot_patient_compounds(all_preds_ci, all_viability, "ASC397", group = "ASC397cse"))
print(plot_patient_compounds(all_preds_ci, all_viability, "ASC399", group = "ASC399"))
print(plot_patient_compounds(all_preds_ci, all_viability, "ASC399", group = "ASC399cse"))
print(plot_across(all_preds_ci, all_viability, "A1"))
print(plot_across(all_preds_ci, all_viability, "A2"))
print(plot_across(all_preds_ci, all_viability, "A3"))
print(plot_across(all_preds_ci, all_viability, "A4"))
## ---- Example usages ----------------------------------------------------------
# print(plot_across(all_preds_ci, all_viability, "A1"))
# print(plot_across(all_preds_pi, all_viability, "A1"))
# print(plot_within(all_preds_ci, all_viability, "A1"))
# print(plot_within(all_preds_pi, all_viability, "A1"))
# print(plot_patient_compounds(all_preds_pi, all_viability, "ASC397", group = "ASC397"))
# print(plot_patient_compounds(all_preds_ci, all_viability, "ASC397"))

# Methods blurb (console + optional file)
cat("\n================  METHODS BLURB  ==================\n")
cat(methods_blurb(cfg), "\n")

# Optional: write to a text file for manuscript methods/supplement
# writeLines(methods_blurb(cfg), "Methods_IC50_modelling.txt")

## ---- Export optional ---------------------------------------------------------
# write_csv(all_ic50, "IC50_results.csv")
# saveRDS(list(ic50 = all_ic50, preds_ci = all_preds_ci, preds_pi = all_preds_pi,
#              viability = all_viability, config = cfg, session = sessionInfo()),
#         "IC50_analysis_bundle.rds")
