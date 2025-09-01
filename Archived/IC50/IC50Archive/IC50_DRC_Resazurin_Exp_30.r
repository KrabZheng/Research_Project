##############################################################################
#   IC50 analysis — v17 (publication-grade, 4 compounds/plate)   2025-08-26
##############################################################################

## ---- Paths ------------------------------------------------------------------
raw_path      <- "IC50/Resazurin_Exp_30/20250704_Resazurin exp 30_CSE_H2O2_KBrO3_Menadione_RSL3_RSL3CSE_Erastin_ErastinCSE.txt"
compound_path <- "IC50/Resazurin_Exp_30/Resazurin_Exp_30_Compound.txt"

# Optional donor maps: tell a compound to borrow its 0-dose ("U") and/or top-dose ("L") anchors
# Keys and values are CompoundIDs like "A1","A2","A3","A4","B1","B2","B3","B4"
ref_map <- c()                                   # legacy (use same donor for both U and L)
min_map <- c("B2" = "B1", "B4" = "B3")           # U (0-dose) donor per compound
max_map <- c("B1" = "B2", "B3" = "B4")           # L (top-dose) donor per compound

validate_donor_maps <- function(..., compounds = NULL, allow_self = FALSE) {
  # ... maps: any number of named character vectors (compound -> donor)
  # compounds: optional character vector of valid compound IDs
  # allow_self: if FALSE, self-references are not allowed
  maps <- list(...)
  all_edges <- purrr::map_dfr(maps, ~tibble(from = names(.x), to = as.character(.x)))
  if (nrow(all_edges) == 0) return(invisible(TRUE))
  # Use all unique compounds if not provided
  if (is.null(compounds)) compounds <- unique(c(all_edges$from, all_edges$to))
  # Check all donors exist
  missing <- setdiff(all_edges$to, compounds)
  if (length(missing)) {
    stop("Donor map error: missing donor(s): ", paste(missing, collapse = ", "))
  }
  # Check for self-references if not allowed
  if (!allow_self) {
    self_refs <- all_edges$from[all_edges$from == all_edges$to]
    if (length(self_refs)) {
      stop("Donor map error: self-reference(s) not allowed: ", paste(self_refs, collapse = ", "))
    }
  }
  # Cycle detection (DFS)
  g <- split(all_edges$to, all_edges$from)
  visited <- setNames(rep(0L, length(compounds)), compounds) # 0=unvisited, 1=visiting, 2=done
  cycle_path <- character()
  has_cycle <- FALSE
  dfs <- function(node, stack = character()) {
    if (visited[[node]] == 1L) {
      cycle_path <<- c(stack, node)
      has_cycle <<- TRUE
      return()
    }
    if (visited[[node]] == 2L) return()
    visited[[node]] <<- 1L
    for (nbr in g[[node]] %||% character()) {
      dfs(nbr, c(stack, node))
      if (has_cycle) return()
    }
    visited[[node]] <<- 2L
  }
  for (v in compounds) {
    if (visited[[v]] == 0L) dfs(v)
    if (has_cycle) break
  }
  if (has_cycle) {
    stop("Donor map error: cycle detected among compounds: ",
         paste(cycle_path, collapse = " -> "))
  }
  invisible(TRUE)
}

# Validate donor maps before proceeding
all_compounds <- unique(c(names(min_map), min_map, names(max_map), max_map, names(ref_map), ref_map))
validate_donor_maps(min_map, max_map, ref_map, compounds = all_compounds, allow_self = FALSE)

# ...existing code...
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
  loess_span         = 0.8           # smoothness for σ(x) loess
)
cfg$bands_method <- "parametric"   # "parametric" or "auto" (allow nonparam fallback)
cfg$model_plateau <- "fix_d"  # "fix_d" (recommended) or "fix_cd" (strict) or "free"
cfg$y_eps        <- 0.25   # clamp means to [0.25, 99.75] % before fitting
# (Lines 30-31 removed—duplicate configuration entries)set.seed(cfg$seed)

## ---- Libraries ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(MASS)
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(purrr); library(tibble); library(drc); library(ggplot2)
})

## ---- Concentration standardiser (µM, %, mass/vol) ----------------------------
`%||%` <- function(a,b) if (is.null(a)) b else a

.standardise_conc <- function(x, unit){
  u <- tolower(gsub("\\s+", "", unit %||% ""))
  u <- gsub("μ","u",u)  # normalise micro symbol

  res <- list(value = x, unit_std = unit, ok = TRUE, convertible = FALSE)

  # molarity → µM
  if (u %in% c("m","mol/l"))                    { res$value <- x*1e6; res$unit_std <- "µM"; res$convertible <- TRUE; return(res) }
  if (u %in% c("mm","mmol/l"))                  { res$value <- x*1e3; res$unit_std <- "µM"; res$convertible <- TRUE; return(res) }
  if (u %in% c("um","µm","umol/l","µmol/l"))    { res$value <- x;     res$unit_std <- "µM"; res$convertible <- TRUE; return(res) }
  if (u %in% c("nm","nmol/l"))                  { res$value <- x*1e-3; res$unit_std <- "µM"; res$convertible <- TRUE; return(res) }
  if (u %in% c("pm","pmol/l"))                  { res$value <- x*1e-6; res$unit_std <- "µM"; res$convertible <- TRUE; return(res) }

  # percentage & common mass/volume (left as-is; cannot convert to µM without MW)
  if (u %in% c("%","percent","%v/v","v/v%","%w/v","w/v%"))  { res$value <- x; res$unit_std <- "%";      return(res) }
  if (u %in% c("mg/ml","mgperml"))                           { res$value <- x; res$unit_std <- "mg/mL";  return(res) }
  if (u %in% c("ug/ml","µg/ml","mcg/ml","ugperml","µgperml")){ res$value <- x; res$unit_std <- "µg/mL";  return(res) }
  if (u %in% c("ng/ml","ngperml"))                           { res$value <- x; res$unit_std <- "ng/mL";  return(res) }

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

## ---- Layout helpers ----------------------------------------------------------
build_layout <- function(path) {
  raw <- read_delim(path, "\t", col_names = FALSE, trim_ws = TRUE, show_col_types = FALSE)
  ids   <- as.character(raw[raw$X1 == "ID",       -1])
  names <- as.character(raw[raw$X1 == "Compound", -1])
  units <- as.character(raw[raw$X1 == "Unit",     -1])
  dose_rows <- raw %>% filter(str_detect(X1, "^[0-9]+$"))
  map_dfr(seq_along(ids), function(j) {
    map_dfr(seq_len(nrow(dose_rows)), function(i) {
      tibble(CompoundID = ids[j], CompoundName = names[j],
             Row = LETTERS[as.integer(dose_rows$X1[i])],
             Concentration = as.numeric(dose_rows[i, j + 1][[1]]), Unit = units[j])
    })
  })
}

diagnose_anchors <- function(v){
  v %>%
    dplyr::summarise(
      n0 = sum(Dose == 0, na.rm = TRUE),
      U = median(U_used, na.rm = TRUE),
      L = median(L_used, na.rm = TRUE),
      Denom = median(Denom_used, na.rm = TRUE),
      y0 = median(Viability[Dose == 0], na.rm = TRUE),
      ytop = median(Viability[Dose == max(Dose, na.rm = TRUE)], na.rm = TRUE),
      .by = c(PlateID, CompoundID, UnitStd)
    ) %>%
    dplyr::arrange(PlateID, CompoundID)
}
# Map column groups 1..4 → "A1..A4" or "B1..B4" based on plate suffix
make_group_map <- function(plate_id) {
  suffix <- stringr::str_sub(plate_id, -1)           # "A" or "B"
  codes  <- if (suffix == "A") paste0("A", 1:4) else paste0("B", 1:4)
  setNames(codes, 1:4)
}
# Vectorised mapping: prefer 'specific', then 'legacy'; default to self
map_compound <- function(x, specific = NULL, legacy = NULL){
  x <- as.character(x)
  out <- rep(NA_character_, length(x))
  if (!is.null(specific) && length(specific) > 0) {
    idx <- match(x, names(specific))
    pick <- unname(specific[idx])
    out <- ifelse(!is.na(idx), pick, NA_character_)
  }
  if (!is.null(legacy) && length(legacy) > 0) {
    idx <- match(x, names(legacy))
    pick <- unname(legacy[idx])
    out <- ifelse(is.na(out) & !is.na(idx), pick, out)
  }
  ifelse(is.na(out), x, out)
}

## ---- Robust viability extraction (with optional donor maps) ------------------
get_viability_data <- function(rows, all_lines, layout,
                               ref_map = NULL,   # legacy: both U & L borrow from here
                               min_map = NULL,   # optional: per-compound U donor
                               max_map = NULL) { # optional: per-compound L donor

  raw <- purrr::map_dfr(rows, function(i){
    plate <- readr::read_delim(I(paste(all_lines[i:(i+8)], collapse = "\n")),
                               "\t", trim_ws = TRUE, show_col_types = FALSE)

    fluo_start <- which(stringr::str_detect(names(plate), "_[Ff]luo$"))[1]
    sig <- plate[, (fluo_start + 1):(fluo_start + 12)] |>
      setNames(sprintf("%02d", 1:12))

    PlateID   <- stringr::str_remove(names(plate)[fluo_start], "_[Ff]luo$")
    BaseGroup <- stringr::str_remove(PlateID, "_?[AB]$")   # patient/condition label
    grp_map   <- make_group_map(PlateID)

    sig |>
      dplyr::mutate(Row = LETTERS[1:8]) |>
      tidyr::pivot_longer(-Row, names_to = "Col", values_to = "Signal",
                          names_transform = list(Col = as.integer)) |>
      dplyr::mutate(
        PlateID     = PlateID,
        ColGroup    = ceiling(Col / 3),
        CompoundID  = grp_map[as.character(ColGroup)],
        Half        = substr(CompoundID, 1, 1),        # "A" or "B"
        Group       = BaseGroup,                       # single patient/condition
        Patient     = BaseGroup                        # keep naming consistent upstream
      ) |>
      dplyr::left_join(layout, by = c("CompoundID", "Row"))
  })

  # Standardise doses (µM, %, mg/mL, ...)
  raw <- convert_conc_to_std(raw)
  if (any(!raw$UnitOK, na.rm = TRUE)) {
    bad <- raw %>% filter(!UnitOK) %>% distinct(CompoundID, Unit) %>% arrange(CompoundID)
    warning("Unrecognised units detected; proceeding with raw values:\n",
            paste0(" - ", bad$CompoundID, ": '", bad$Unit, "'", collapse = "\n"))
  }

  # Each compound must use only ONE UnitStd
  per_cmpd <- raw %>% distinct(CompoundID, UnitStd) %>% count(CompoundID, name = "n_units")
  if (any(per_cmpd$n_units > 1)) {
    offenders <- per_cmpd %>% filter(n_units > 1) %>% pull(CompoundID)
    stop("Mixed units found within the same compound: ", paste(offenders, collapse = ", "),
         ". Please fix the input layout/units for these plates.")
  }

  # Choose donor IDs for U (0-dose) and L (top-dose)
  choose_map <- function(x, map_specific, map_legacy) {
    out <- if (!is.null(map_specific) && length(map_specific) && x %in% names(map_specific))
      as.character(map_specific[x]) else NA_character_
    ifelse(!is.na(out),
           out,
           ifelse(!is.null(map_legacy) && length(map_legacy) && x %in% names(map_legacy),
                  as.character(map_legacy[x]), x))
  }

  raw <- raw %>%
    mutate(
      BaseU_ID = map_compound(CompoundID, min_map, ref_map),
      BaseL_ID = map_compound(CompoundID, max_map, ref_map)
    )
  # Compute anchors on donors (median for robustness)
  U_tbl <- raw %>%
    group_by(PlateID, BaseU_ID) %>%
    summarise(U = median(Signal[is.finite(Signal) & Dose == 0], na.rm = TRUE),
              .groups = "drop")

  L_tbl <- raw %>%
    group_by(PlateID, BaseL_ID) %>%
    summarise(
      top_dose = suppressWarnings(max(Dose[is.finite(Dose)], na.rm = TRUE)),
      L = median(Signal[is.finite(Signal) & Dose >= top_dose], na.rm = TRUE),
      .groups = "drop"
    )

  # Self fallbacks if donor wells are missing
  self_tbl <- raw %>%
    group_by(PlateID, CompoundID) %>%
    summarise(
      U_self = median(Signal[is.finite(Signal) & Dose == 0], na.rm = TRUE),
      top_dose = suppressWarnings(max(Dose[is.finite(Dose)], na.rm = TRUE)),
      L_self = median(Signal[is.finite(Signal) & Dose >= top_dose], na.rm = TRUE),
      .groups = "drop"
    )

  viab <- raw %>%
    left_join(U_tbl,    by = c("PlateID" = "PlateID", "BaseU_ID" = "BaseU_ID")) %>%
    left_join(L_tbl,    by = c("PlateID" = "PlateID", "BaseL_ID" = "BaseL_ID")) %>%
    left_join(self_tbl, by = c("PlateID", "CompoundID")) %>%
    mutate(
      U_use = ifelse(is.finite(U), U, U_self),
      L_use = ifelse(is.finite(L), L, L_self),
      Denom = pmax(U_use - L_use, 1e-6),
      Viability = pmin(pmax((Signal - L_use) / Denom * 100, 0), 100)
    ) %>%
    filter(is.finite(Dose) & Dose > 0) %>%
    select(-matches("^top_dose($|\\.)")) %>%
    mutate(U_used = U_use, L_used = L_use, Denom_used = pmax(U_use - L_use, 1e-6))

  viab
}

## ---- Replicate validation ----------------------------------------------------
validate_replicates <- function(viab){
  rep_chk <- viab %>%
    filter(is.finite(Dose)) %>%
    count(PlateID, CompoundID, Dose, name = "n_rep") %>%
    summarise(min_rep = min(n_rep), max_rep = max(n_rep), .by = c(PlateID, CompoundID)) %>%
    mutate(ok = min_rep == max_rep & min_rep %in% c(2,3,4))
  if (any(!rep_chk$ok)) {
    message("⚠ Replicate inconsistency detected:")
    print(rep_chk %>% filter(!ok))
  }
  invisible(rep_chk)
}

## ---- Fit utilities: robust aggregated LL.4 + PD/NP bootstraps ---------------
make_pd <- function(S, eps = 1e-8){
  S <- as.matrix(S); S <- (S + t(S))/2
  ev <- eigen(S, symmetric = TRUE)
  ev$values[ev$values < eps] <- eps
  ev$vectors %*% diag(ev$values, nrow = length(ev$values)) %*% t(ev$vectors)
}

# Single robust fit: trim 3*MAD per dose, aggregate to means (weights = n)
resample_fit_once <- function(df_repl){
  df <- df_repl[is.finite(df_repl$Dose) & is.finite(df_repl$Viability) & df_repl$Dose > 0, ]
  if (!nrow(df)) return(NULL)

  # Trim outliers per dose, then aggregate
  df_trim <- df %>%
    dplyr::group_by(Dose) %>%
    dplyr::mutate(med = median(Viability, na.rm = TRUE),
                  madv = mad(Viability, constant = 1.4826, na.rm = TRUE),
                  keep = ifelse(is.finite(madv) & madv > 0, abs(Viability - med) <= 3*madv, TRUE)) %>%
    dplyr::ungroup() %>% dplyr::filter(keep)

  df_ag <- df_trim %>%
    dplyr::group_by(Dose) %>%
    dplyr::summarise(Viability = mean(Viability, na.rm=TRUE), wt = dplyr::n(), .groups = "drop")

  doses <- sort(unique(df_ag$Dose))
  if (length(doses) < cfg$min_unique_doses) return(NULL)

  # Clamp Y just off the boundaries to stabilise drm’s derivatives
  df_ag$Viability <- pmin(pmax(df_ag$Viability, cfg$y_eps), 100 - cfg$y_eps)

  # Heuristic initial values for LL.4
  ll4_starts <- function(dat, fix_c, fix_d){
    y <- dat$Viability; x <- dat$Dose
    # rough plateaus
    c0 <- if (is.na(fix_c)) max(0, min(y, na.rm=TRUE)) else NA
    d0 <- if (is.na(fix_d)) min(100, max(y, na.rm=TRUE)) else NA
    # ED50 guess: dose where y crosses mid between plateaus
    ylo <- if (is.na(fix_c)) c0 else 0
    yhi <- if (is.na(fix_d)) d0 else 100
    ymid <- (ylo + yhi)/2
    ord <- order(x)
    xc <- try(approx(x[ord], y[ord], xout = ymid, ties = "ordered")$x, silent = TRUE)
    e0 <- if (!inherits(xc,"try-error") && is.finite(xc)) xc else exp(mean(log(range(x, na.rm=TRUE))))
    b0 <- 1.2
    # build vector in drm parameter order (b, c?, d?, e)
    if (!is.na(fix_c) && !is.na(fix_d)) return(c(b0, e0))
    if (!is.na(fix_c) &&  is.na(fix_d)) return(c(b0, d0, e0))
    if ( is.na(fix_c) && !is.na(fix_d)) return(c(b0, c0, e0))
    c(b0, c0, d0, e0)
  }

  fit_variant <- function(fix_c, fix_d){
    st <- ll4_starts(df_ag, fix_c, fix_d)
    suppressWarnings(tryCatch(
      drm(Viability ~ Dose, data = df_ag, weights = wt,
          fct = LL.4(fixed = c(NA, fix_c, fix_d, NA)),
          start = st,
          control = drmc(maxIt = 2000, relTol = 1e-8)),
      error = function(e) NULL))
  }

  # Candidate set (respect cfg$model_plateau)
  cand <- switch(cfg$model_plateau,
    "fix_cd" = list(fix_cd = fit_variant(0, 100)),
    "fix_d"  = list(fix_d  = fit_variant(NA, 100),
                    fix_cd = fit_variant(0, 100)),
    "free"   = list(free   = fit_variant(NA, NA),
                    fix_c  = fit_variant(0,  NA),
                    fix_d  = fit_variant(NA, 100),
                    fix_cd = fit_variant(0,  100))
  )
  keep <- vapply(cand, Negate(is.null), logical(1))
  if (!any(keep)) return(NULL)
  cand <- cand[keep]

  n <- sum(df_ag$wt)
  aicc <- function(f){ k <- length(coef(f)); den <- n - k - 1; AIC(f) + if (den>0) 2*k*(k+1)/den else Inf }
  aics <- vapply(cand, aicc, numeric(1))
  ks   <- vapply(cand, function(f) length(coef(f)), numeric(1))
  best <- which.min(aics); ties <- which(aics - aics[best] <= cfg$aicc_tie_tol)
  if (length(ties) > 1) best <- ties[which.min(ks[ties])]
  best_fit <- cand[[best]]
  attr(best_fit, "variant") <- names(cand)[best]
  attr(best_fit, "data_ag") <- df_ag
  best_fit
}

safe_fit <- function(df){ resample_fit_once(df) }

.ll4_mean <- function(x,b,c,d,e) c + (d - c) / (1 + (x / e)^b)

sigma_x_from_loess <- function(fit, grid, span = cfg$loess_span){
  res <- residuals(fit); sigma <- sd(res, na.rm = TRUE); if (!is.finite(sigma)) sigma <- 0
  if (length(res) >= 6) {
    tmp <- data.frame(fitv = fitted(fit), absr = abs(res))
    lo  <- stats::loess(absr ~ fitv, data = tmp, span = span, degree = 1)
    mu  <- .ll4_mean(grid,
                     coef(fit)[startsWith(names(coef(fit)), "b:")],
                     ifelse("c:(Intercept)" %in% names(coef(fit)), coef(fit)["c:(Intercept)"], 0),
                     ifelse("d:(Intercept)" %in% names(coef(fit)), coef(fit)["d:(Intercept)"], 100),
                     coef(fit)[startsWith(names(coef(fit)), "e:")])
    sig_x <- predict(lo, newdata = data.frame(fitv = mu))
    if (anyNA(sig_x)) {
      xs <- sort(unique(tmp$fitv)); ys <- predict(lo, newdata = data.frame(fitv = xs))
      sig_x <- approx(xs, ys, xout = mu, rule = 2)$y
    }
    pmax(sig_x, 1e-6)
  } else rep(sigma, length(grid))
}

bootstrap_mean_matrix <- function(fit, df, grid, B,
                                  b_bounds = cfg$b_draw_bounds,
                                  e_bounds = NULL,
                                  method = cfg$bands_method){
  # Parametric (preferred)
  V <- try(vcov(fit), silent = TRUE)
  if (!inherits(V, "try-error") && is.matrix(V) && all(is.finite(V))) {
    V <- (V + t(V))/2
    ev <- eigen(V, symmetric = TRUE)
    lam <- pmax(ev$values, 1e-10)            # ridge to avoid zero/negative eigenvalues
    Vpd <- ev$vectors %*% diag(lam, length(lam)) %*% t(ev$vectors)

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

  if (method == "parametric") return(matrix(numeric(0), nrow = length(grid), ncol = 0))

  # Optional: nonparametric fallback
  Bnp <- min(B, 800)
  preds <- matrix(NA_real_, nrow = length(grid), ncol = Bnp)
  for (i in seq_len(Bnp)) {
    dfb <- df %>% dplyr::group_by(Dose) %>% dplyr::slice_sample(n = dplyr::n(), replace = TRUE) %>%
           dplyr::ungroup()
    fb  <- resample_fit_once(dfb); if (is.null(fb)) next
    preds[, i] <- as.numeric(predict(fb, newdata = data.frame(Dose = grid)))
  }
  preds
}

quiet_ED_delta <- function(fit, p = 50){
  ed <- NULL
  ok <- try(suppressWarnings(utils::capture.output(
    ed <- drc::ED(fit, p, interval = "delta")
  )), silent = TRUE)
  if (inherits(ok, "try-error")) return(NULL)
  ed
}

## ---- Predictions with percentile bands ---------------------------------------
predict_curves <- function(viab, n = 200, band = c("ci", "pi"), level = 0.95){
  band <- match.arg(band); a <- (1 - level) / 2; probs <- c(a, 1 - a)

  viab %>%
    group_by(PlateID, Patient, Group, CompoundID, CompoundName, Unit, UnitStd, Half) %>%
    nest() %>%
    mutate(
      fit = map(data, safe_fit),
      pred = map2(fit, data, ~{
        fit <- .x; df <- .y
        if (is.null(fit)) return(NULL)

        concs <- sort(unique(df$Dose[df$Dose > 0 & is.finite(df$Dose)]))
        if (length(concs) < cfg$min_unique_doses) return(NULL)
        grid <- exp(seq(log(min(concs)), log(max(concs)), length.out = n))

        e_bounds <- range(concs) * cfg$e_bounds_factor
        M <- bootstrap_mean_matrix(fit, df, grid,
                           B = if (band == "ci") cfg$B_ci else cfg$B_pi,
                           b_bounds = cfg$b_draw_bounds, e_bounds = e_bounds,
                           method = cfg$bands_method)
        if (is.null(M) || length(M) == 0 || is.null(dim(M)) || ncol(M) == 0) {
  mu <- as.numeric(predict(fit, newdata = data.frame(Dose = grid)))
  return(tibble(Dose = grid, Viability = pmin(pmax(mu, 0), 100),
                Lower = pmin(pmax(mu, 0), 100), Upper = pmin(pmax(mu, 0), 100)))
}

        mu <- rowMeans(M, na.rm = TRUE)
        if (band == "ci") {
          lo <- apply(M, 1, quantile, probs = probs[1], na.rm = TRUE)
          hi <- apply(M, 1, quantile, probs = probs[2], na.rm = TRUE)
        } else {
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
qc_metrics <- function(df, fit){
  res <- try(residuals(fit), silent = TRUE); if (inherits(res, "try-error")) res <- NA_real_
  rmse <- if (all(is.finite(res))) sqrt(mean(res^2)) else NA_real_
  sse <- if (all(is.finite(res))) sum(res^2) else NA_real_
  sst <- if (all(is.finite(df$Viability))) sum((df$Viability - mean(df$Viability, na.rm=TRUE))^2) else NA_real_
  r2  <- if (is.finite(sse) && is.finite(sst) && sst > 0) 1 - sse/sst else NA_real_
  dsum <- df %>% filter(Dose > 0, is.finite(Dose), is.finite(Viability)) %>%
    summarise(mv = mean(Viability, na.rm=TRUE), .by = Dose) %>% arrange(Dose)
  rho <- if (nrow(dsum) >= 4) suppressWarnings(cor(dsum$mv, log(dsum$Dose), method="spearman")) else NA_real_
  tibble(RMSE = rmse, R2 = r2, Spearman = rho)
}

compute_ic50 <- function(viab){
  viab %>%
    group_by(PlateID, Patient, Group, CompoundID, CompoundName, Unit, UnitStd, Half) %>%
    group_modify(~{
      f <- safe_fit(.x)
      if (is.null(f)) return(tibble(
        IC50 = NA, IC50_L = NA, IC50_U = NA, in_range = NA,
        b = NA, c = NA, d = NA, aicc = NA, variant = NA, converged = FALSE,
        RMSE = NA, R2 = NA, Spearman = NA, n_doses = n_distinct(.x$Dose[.x$Dose > 0])
      ))
      cf <- coef(f); nm <- names(cf); getp <- function(k) if (k %in% nm) unname(cf[k]) else NA_real_
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
        if (!is.finite(ic50_l) || !is.finite(ic50_u)) {
          Bnp <- 800
          ee  <- rep(NA_real_, Bnp)
          for (i in seq_len(Bnp)) {
            dfb <- .x %>% group_by(Dose) %>% slice_sample(n = n(), replace = TRUE) %>% ungroup()
            fb  <- resample_fit_once(dfb); if (is.null(fb)) next
            cf2 <- coef(fb); nm2 <- names(cf2)
            ee[i] <- if ("e:(Intercept)" %in% nm2) unname(cf2["e:(Intercept)"]) else NA_real_
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
      converged <- tryCatch(isTRUE(f$fit$convergence == 0), error = function(e) TRUE)
      qc <- qc_metrics(.x, f)

      tibble(
        IC50 = ic50, IC50_L = ic50_l, IC50_U = ic50_u,
        in_range = is.finite(ic50) && ic50 > min(concs) && ic50 < max(concs),
        b = getp("b:(Intercept)"),
        c = getp("c:(Intercept)"),
        d = getp("d:(Intercept)"),
        aicc = aicc, variant = attr(f, "variant"), converged = converged,
        RMSE = qc$RMSE, R2 = qc$R2, Spearman = qc$Spearman, n_doses = length(concs)
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
      converged = coalesce(converged, FALSE),
      support_ok = coalesce(support_ok, FALSE),
      slope_ok = coalesce(slope_ok, FALSE),
      mono_ok = coalesce(mono_ok, FALSE),
      in_range = coalesce(in_range, FALSE),
      ic50_ok = is.finite(IC50),
      pass = converged & support_ok & slope_ok & mono_ok & in_range & ic50_ok,
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
    arrange(PlateID, CompoundID, Patient, Group)
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

plot_compound <- function(preds, viab, cmpd,
                          colour_by = c("PlateID","Group"),
                          plate = NULL, facet_by_plate = FALSE,
                          show_legend = TRUE) {
  colour_by <- match.arg(colour_by)

  df_line <- dplyr::filter(preds, CompoundID == cmpd)
  df_pts  <- dplyr::filter(viab,  CompoundID == cmpd, Dose > 0)

  if (!is.null(plate)) {
    df_line <- dplyr::filter(df_line, PlateID == plate)
    df_pts  <- dplyr::filter(df_pts,  PlateID == plate)
  }
  if (!nrow(df_line)) stop("No predictions for ", cmpd,
                           if (!is.null(plate)) paste0(" on plate ", plate))

  # Map a single variable to both colour and fill (like Exp-36)
  if (colour_by == "Group") {
    df_line <- dplyr::mutate(df_line, Colour = interaction(Patient, Group, sep=" / ", drop=TRUE))
    df_pts  <- dplyr::mutate(df_pts,  Colour = interaction(Patient, Group, sep=" / ", drop=TRUE))
    legend_name <- "Patient / Group"
  } else {
    df_line <- dplyr::mutate(df_line, Colour = PlateID)
    df_pts  <- dplyr::mutate(df_pts,  Colour = PlateID)
    legend_name <- "Plate"
  }

  u <- unique(df_line$UnitStd)
  xlab <- if (length(u) == 1) sprintf("Dose [%s]", u) else "Dose [mixed units]"

  p <- ggplot() +
    geom_ribbon(data = df_line,
                aes(Dose, ymin = Lower, ymax = Upper, fill = Colour),
                alpha = .22, colour = NA) +
    geom_point(data = df_pts,
               aes(Dose, Viability, colour = Colour),
               size = 1.6, alpha = .65) +
    geom_line(data = df_line,
              aes(Dose, Viability, colour = Colour),
              linewidth = .9) +
    scale_x_log10() +
    coord_cartesian(ylim = c(0, 100)) +
    .plot_base() +
    labs(title = sprintf("%s (%s)", unique(df_line$CompoundName), cmpd),
         x = xlab, y = "Viability (%)",
         colour = legend_name, fill = legend_name) +
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2")

  if (!show_legend) p <- p + theme(legend.position = "none")
  if (is.null(plate) && facet_by_plate && dplyr::n_distinct(df_line$PlateID) > 1) {
    p <- p + facet_wrap(~ PlateID, scales = "free_x")
  }
  p
}

plot_plate_overview <- function(preds, viab, half = NULL, colour_by = c("Group","PlateID"),
                                show_legend = FALSE) {
  colour_by <- match.arg(colour_by)

  df_line <- preds
  df_pts  <- viab %>% dplyr::filter(Dose > 0)
  if (!is.null(half)) {
    df_line <- dplyr::filter(df_line, Half == half)
    df_pts  <- dplyr::filter(df_pts,  Half == half)
  }

  if (colour_by == "Group") {
    df_line <- dplyr::mutate(df_line, Colour = interaction(Patient, Group, sep=" / ", drop=TRUE))
    df_pts  <- dplyr::mutate(df_pts,  Colour = interaction(Patient, Group, sep=" / ", drop=TRUE))
    legend_name <- "Patient / Group"
  } else {
    df_line <- dplyr::mutate(df_line, Colour = PlateID)
    df_pts  <- dplyr::mutate(df_pts,  Colour = PlateID)
    legend_name <- "Plate"
  }

  ggplot() +
    geom_ribbon(data = df_line, aes(Dose, ymin = Lower, ymax = Upper, fill = Colour),
                alpha = .20, colour = NA) +
    geom_point(data = df_pts, aes(Dose, Viability, colour = Colour),
               size = 1.1, alpha = .55) +
    geom_line(data = df_line, aes(Dose, Viability, colour = Colour),
              linewidth = .7) +
    scale_x_log10() +
    coord_cartesian(ylim = c(0, 100)) +
    facet_grid(PlateID ~ CompoundID + CompoundName, scales = "free_x") +
    .plot_base() +
    labs(title = "Plate overview (per PlateID × Compound)",
         x = "Dose [per panel unit]", y = "Viability (%)",
         colour = legend_name, fill = legend_name) +
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
# Debug: Check specific anchor diagnostics (uncomment if needed)
# diagnose_anchors(all_viability) %>% dplyr::filter(PlateID == "ASC390cse_A", CompoundID == "A1")}

## ---- Methods blurb -----------------------------------------------------------
methods_blurb <- function(cfg){
  N  <- cfg$min_unique_doses
  bL <- cfg$slope_bounds[1]; bU <- cfg$slope_bounds[2]
  B  <- cfg$B_ci; span <- cfg$loess_span
  paste(
"Dose–response modelling and IC50 estimation.",
"",
"Signals were normalised per plate and compound using a two-point anchor. ",
"The upper anchor (U) was the median signal at 0-dose (optionally borrowed from a donor compound); ",
"the lower anchor (L) was the median signal at the highest tested dose (optionally borrowed). ",
"Viability was 100 × (S − L)/(U − L) truncated to [0, 100]. Robust fallbacks were used if donors were absent.",
"",
"Concentrations were standardised (µM when molarity was provided; otherwise original units such as % or mg/mL ",
"were preserved). For each compound, a four-parameter log-logistic model (LL.4) was fitted to positive doses ",
"(log-scaled x). The selected variant (free; c fixed to 0; d fixed to 100; or both) was chosen by AICc, ",
"preferring the simpler model when ΔAICc ≤ 2. Fits required at least ", N, " unique positive doses and used increased iteration limits.",
"",
"Uncertainty in the mean curve was quantified via percentile bootstrap from the asymptotic parameter covariance ",
"(", B, " draws; positive-definite repair applied as needed). Predictive 95% intervals incorporated heteroscedastic residual variance ",
"modelled by LOESS of |residuals| versus fitted values (span = ", sprintf("%.2f", span), "), with noise added in each bootstrap path.",
"",
"The IC50 (ED50) was extracted with 95% confidence intervals via the delta method when available, or bootstrap when not. ",
"An in-range flag indicates whether ED50 lies within the tested dose range.",
"",
"Quality control comprised replicate-count checks, convergence status, dose support (≥ ", N, "), slope bounds (|b| ∈ [",
sprintf("%.1f", bL), ", ", sprintf("%.1f", bU), "]), pseudo-R² and RMSE, and monotonicity (Spearman correlation between mean viability and log dose).",
sep = ""
  )
}

## ---- Pipeline ----------------------------------------------------------------
layout_long  <- build_layout(compound_path)
all_lines <- read_lines(raw_path)
rows <- which(str_detect(all_lines, regex("_fluo", ignore_case=TRUE)))
rows <- rows[rows + 8 <= length(all_lines)]
if (length(rows) == 0) stop("No \"_Fluo\" blocks found.")

# Viability with robust anchors (+ donor maps) and unit awareness
all_viability <- get_viability_data(rows, all_lines, layout_long, ref_map, min_map, max_map)

# Replicate validation
validate_replicates(all_viability)

# IC50 table with CIs and QC flags
all_ic50 <- compute_ic50(all_viability)

# QC summary (prints) + full annotated table
qc_df <- qc_summary(all_ic50)
# write_csv(qc_df, "IC50_QC_report_plate4x.csv")   # optional export
# write_csv(all_ic50, "IC50_results_plate4x.csv")  # optional export
diagnose_anchors(all_viability) %>% dplyr::filter(PlateID == "ASC390cse_A", CompoundID == "A1")

# Predictions (percentile bands)
all_preds_ci <- predict_curves(all_viability, band="ci", level=.95)
all_preds_pi <- predict_curves(all_viability, band="pi", level=.95)


# Methods text (console + optional file)
cat("\n================  METHODS BLURB  ==================\n")
cat(methods_blurb(cfg), "\n")
# writeLines(methods_blurb(cfg), "Methods_IC50_plate4x.txt")

# CI ribbons coloured by Group (closest to Exp-36 look)
print(plot_compound(all_preds_ci, all_viability, "A1", colour_by = "Group"))
print(plot_compound(all_preds_ci, all_viability, "A2", colour_by = "Group"))
print(plot_compound(all_preds_ci, all_viability, "A3", colour_by = "Group"))
print(plot_compound(all_preds_ci, all_viability, "A4", colour_by = "Group"))
print(plot_compound(all_preds_ci, all_viability, "B1", colour_by = "Group"))
print(plot_compound(all_preds_ci, all_viability, "B2", colour_by = "Group"))
print(plot_compound(all_preds_ci, all_viability, "B3", colour_by = "Group"))
print(plot_compound(all_preds_ci, all_viability, "B4", colour_by = "Group"))

# Overlay multiple plates for the same compound, coloured by PlateID
# print(plot_compound(all_preds_ci, all_viability, "A1", colour_by = "PlateID"))

# Facet each plate into its own panel (still coloured)
# print(plot_compound(all_preds_ci, all_viability, "A1", colour_by = "PlateID", facet_by_plate = TRUE))

# Plate overview, coloured by Group (or PlateID)
# print(plot_plate_overview(all_preds_ci, all_viability, colour_by = "Group"))
# print(plot_plate_overview(all_preds_ci, all_viability, colour_by = "PlateID"))
write_csv(all_ic50, "Resazurin_Exp30_IC50_results.csv")