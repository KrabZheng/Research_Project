##############################################################################
#   IC50 analysis — v16  (CI/PI ribbons with fast fallback)       2025-08-05
##############################################################################

## Paths
raw_path      <- "IC50/Resazurin_Exp_30/20250704_Resazurin exp 30_CSE_H2O2_KBrO3_Menadione_RSL3_RSL3CSE_Erastin_ErastinCSE.txt"
compound_path <- "IC50/Resazurin_Exp_30/Resazurin_Exp_30_Compound.txt"

# Optional - Allows using different min/max for defined compounds - "target" = "reference"
ref_map <- c()
min_map <- c("B2" = "B1", "B4" = "B3")  # optional: per-compound blank donor
max_map <- c("B1" = "B2", "B3" = "B4")  # optional: per-compound max donor

## Libraries
suppressPackageStartupMessages({
  library(MASS)
  library(readr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(drc)
  library(ggplot2)
  library(dplyr)
})

## Layout helpers
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

# For a plate “P676_A”, return c(A1 = "A1", …, A4 = "A4")
make_group_map <- function(plate_id) {
  suffix <- stringr::str_sub(plate_id, -1)           # "A" or "B"
  codes  <- if (suffix == "A") paste0("A", 1:4) else paste0("B", 1:4)
  setNames(codes, 1:4)
}

# Pick best LL.4 variant per group by AICc
.safe_fit_best <- function(df) {
  # use rows actually passed to drm
  df <- df[is.finite(df$Concentration) & is.finite(df$Viability) & df$Concentration > 0, ]
  doses <- unique(df$Concentration)
  if (length(doses) < 3) return(NULL)
  n <- nrow(df)  # <-- don't call nobs(); drm doesn't define it

  fit_variant <- function(fix_c = NA, fix_d = NA) {
    suppressWarnings(tryCatch(
                              drm(Viability ~ Concentration, data = df, fct = LL.4(fixed = c(NA, fix_c, fix_d, NA))),
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

  aicc <- function(f) {
    k <- length(coef(f))
    den <- n - k - 1
    AIC(f) + if (den > 0) 2*k*(k+1)/den else Inf
  }
  aics <- vapply(cand, aicc, numeric(1))
  ks   <- vapply(cand, function(f) length(coef(f)), numeric(1))

  best <- which.min(aics)
  ties <- which(aics - aics[best] <= 2)
  if (length(ties) > 1) best <- ties[ which.min(ks[ties]) ]

  best_fit <- cand[[best]]
  attr(best_fit, "variant") <- names(cand)[best]
  best_fit
}

# Use the best fit everywhere
safe_fit <- .safe_fit_best

## Extract viability
get_viability_data <- function(rows, all_lines, layout,
                               ref_map = NULL,   # legacy: both min & max borrow from here
                               min_map = NULL,   # optional: per-compound blank donor
                               max_map = NULL) { # optional: per-compound max donor
  raw <- purrr::map_dfr(rows, function(i){

    plate <- readr::read_delim(
      I(paste(all_lines[i:(i+8)], collapse = "\n")),
      "\t", trim_ws = TRUE, show_col_types = FALSE
    )

    fluo_start <- which(stringr::str_detect(names(plate), "_[Ff]luo$"))[1]
    sig <- plate[, (fluo_start + 1):(fluo_start + 12)] |>
      setNames(sprintf("%02d", 1:12))

    PlateID <- stringr::str_remove(names(plate)[fluo_start], "_[Ff]luo$")
    grp_map <- make_group_map(PlateID)

    sig |>
      dplyr::mutate(Row = LETTERS[1:8]) |>
      tidyr::pivot_longer(-Row, names_to = "Col", values_to = "Signal",
                          names_transform = list(Col = as.integer)) |>
      dplyr::mutate(
        PlateID     = PlateID,
        ColGroup    = ceiling(Col / 3),
        CompoundID  = grp_map[as.character(ColGroup)],
        Group       = stringr::str_remove(PlateID, "_?[AB]$"),
        Patient     = stringr::str_remove(Group,
                          stringr::regex("cse$", ignore_case = TRUE))
      ) |>
      dplyr::left_join(layout, by = c("CompoundID", "Row"))
  })

  # --- mapping priority: explicit min/max > ref_map > self ---
  choose_map <- function(x, map_specific, map_legacy) {
    if (!is.null(map_specific) && length(map_specific)) {
      out <- ifelse(x %in% names(map_specific),
                    as.character(map_specific[x]),
                    NA_character_)
    } else out <- NA_character_

    ifelse(!is.na(out), out,
           ifelse(!is.null(map_legacy) && length(map_legacy) && x %in% names(map_legacy),
                  as.character(map_legacy[x]),
                  x))
  }

  raw <- raw %>%
    mutate(
      BaseMinID = choose_map(CompoundID, min_map, ref_map),
      BaseMaxID = choose_map(CompoundID, max_map, ref_map)
    )

  # References computed once per plate & chosen base id
  blank_tbl <- raw %>%
    group_by(PlateID, BaseMinID) %>%
    summarise(Blank_ref = mean(Signal[Concentration == 0], na.rm = TRUE),
              .groups = "drop")

  max_tbl <- raw %>%
    group_by(PlateID, BaseMaxID) %>%
    summarise(Max_ref = mean(Signal[Concentration == max(Concentration, na.rm = TRUE)],
                             na.rm = TRUE),
              .groups = "drop")

  # Self fallbacks in case borrowed wells are missing
  self_tbl <- raw %>%
    group_by(PlateID, CompoundID) %>%
    summarise(Blank_self = mean(Signal[Concentration == 0], na.rm = TRUE),
              Max_self   = mean(Signal[Concentration == max(Concentration, na.rm = TRUE)],
                                na.rm = TRUE),
              .groups = "drop")

  normed <- raw %>%
    left_join(blank_tbl, by = c("PlateID", "BaseMinID")) %>%
    left_join(max_tbl,   by = c("PlateID", "BaseMaxID")) %>%
    left_join(self_tbl,  by = c("PlateID", "CompoundID")) %>%
    mutate(
      Blank = dplyr::coalesce(Blank_ref, Blank_self),
      Max   = dplyr::coalesce(Max_ref,   Max_self),
      Signal_bc   = Signal - Blank,
      denom       = Max - Blank,
      Signal_norm = ifelse(is.finite(denom) & denom != 0, Signal_bc / denom * 100, NA_real_),
      Viability   = pmin(pmax(100 - Signal_norm, 0), 100)
    ) %>%
    filter(!is.na(Concentration) & Concentration > 0) %>%
    select(-Blank_ref, -Max_ref, -Blank_self, -Max_self, -denom)

  normed
}

## IC50 table
compute_ic50 <- function(viab){
  viab %>%
    dplyr::group_by(Patient, Group, CompoundID, CompoundName, Unit) %>%
    dplyr::group_modify(~{
      f <- safe_fit(.x)
      if (is.null(f)) return(tibble(slope = NA,
                                    lower = NA,
                                    upper = NA,
                                    ic50  = NA))
      cf <- coef(f)
      # order in drc is b, [c], [d], e (plateaus may be absent if fixed)
      slope <- unname(cf[startsWith(names(cf), "b:")][1])
      lower <- if ("c:(Intercept)" %in% names(cf)) unname(cf["c:(Intercept)"]) else NA
      upper <- if ("d:(Intercept)" %in% names(cf)) unname(cf["d:(Intercept)"]) else NA
      ic50  <- unname(cf[startsWith(names(cf), "e:")][1])
      tibble(slope = slope,
             lower = lower,
             upper = upper,
             ic50 = ic50)
    }) %>%
    dplyr::ungroup() %>%
    dplyr::rename(IC50 = ic50)
}

## Mean + SE machinery
.ll4_mean <- function(x,b,c,d,e) c + (d - c) / (1 + (x / e)^b)

param_boot_se_curve <- function(fit, grid, B=3000, b_bounds=c(0.1,50), e_bounds=NULL){
  cf <- coef(fit); V <- try(vcov(fit), silent=TRUE); if (inherits(V,"try-error")) return(rep(NA,length(grid)))
  nm <- names(cf); draws <- MASS::mvrnorm(B, mu=cf, Sigma=V); ok <- is.finite(rowSums(draws))
  if ("b:(Intercept)" %in% nm) {
    b  <- draws[, nm == "b:(Intercept)"]
    ok <- ok & b > b_bounds[1] & b < b_bounds[2]
  }
  if ("e:(Intercept)" %in% nm) {
    e <- draws[, nm == "e:(Intercept)"]
    ok <- ok & e > 0
    if (!is.null(e_bounds)) ok <- ok & e > e_bounds[1] & e < e_bounds[2]
  }
  draws <- draws[ok,,drop=FALSE]; if (!nrow(draws)) return(rep(NA,length(grid)))
  getp <- function(letter,default) {
    i <- which(startsWith(nm, paste0(letter, ":")))
    if (length(i)) cf[i] else default
  }
  c0 <- getp("c",0); d0 <- getp("d",100); ib <- which(startsWith(nm,"b:")); ie <- which(startsWith(nm,"e:"))
  Y <- apply(draws,1,function(p){
    b<-p[ib]; e<-p[ie];
    c<-if("c:(Intercept)"%in%nm) p[nm=="c:(Intercept)"] else c0
    d<-if("d:(Intercept)"%in%nm) p[nm=="d:(Intercept)"] else d0
    .ll4_mean(grid,b,c,d,e)
  })
  if (is.null(dim(Y))) Y <- matrix(Y, ncol=1)
  apply(Y,1,sd,na.rm=TRUE)
}

predict_curves <- function(viab, n=200, band=c("ci","pi"), level=0.95,
                           hetero_pi=TRUE, param_B=3000, clamp_bands=FALSE){
  set.seed(1L)
  band <- match.arg(band)
  z <- qnorm(1 - (1 - level) / 2)

  viab %>%
    dplyr::group_by(Patient, Group, CompoundID, CompoundName, Unit) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      fit  = purrr::map(data, safe_fit),
      pred = purrr::map2(fit, data, ~{
        fit <- .x; df <- .y
        if (is.null(fit)) return(NULL)

        concs <- sort(unique(df$Concentration[df$Concentration>0 & is.finite(df$Concentration)]))
        if (length(concs) < 2) return(NULL)
        grid <- exp(seq(log(min(concs)), log(max(concs)), length.out=n))
        mu   <- as.numeric(predict(fit, newdata=data.frame(Concentration=grid)))

        # smooth heteroscedastic σ(x) for PI  (no warnings, no jagged edges)
        res <- residuals(fit); sigma <- sd(res, na.rm=TRUE); if (!is.finite(sigma)) sigma <- 0
        if (hetero_pi && length(res) >= 6) {
          tmp <- data.frame(fitv=fitted(fit), absr=abs(res))
          lo  <- stats::loess(absr ~ fitv, data=tmp, span=0.8, degree=1)
          sig_x <- predict(lo, newdata=data.frame(fitv=mu))
          if (anyNA(sig_x)) {                         # extrapolate smoothly at edges
            xs <- sort(unique(tmp$fitv)); ys <- predict(lo, newdata=data.frame(fitv=xs))
            sig_x <- approx(xs, ys, xout=mu, rule=2)$y
          }
          sig_x <- pmax(sig_x, 1e-6)
        } else sig_x <- rep(sigma, length(mu))

        # SE for mean curve: use parametric bootstrap always (smooth, no kinks)
        V <- try(as.matrix(vcov(fit)), silent = TRUE)
        rho_max <- if (!inherits(V, "try-error") && is.matrix(V) && ncol(V) > 1) {
          max(abs(cov2cor(V)[upper.tri(V)]), na.rm = TRUE)
        } else 0
        e_bounds <- range(concs) * c(0.25, 4)
        well_cond <- is.matrix(V) && all(is.finite(V)) && is.finite(kappa(V)) &&
          kappa(V) < 1e5 && rho_max < 0.95
        B_use <- if (well_cond) ceiling(param_B/2) else param_B
        se_mu <- param_boot_se_curve(fit, grid, B = B_use, e_bounds = e_bounds)

        lo <- mu - z * (if (band=="ci") se_mu else sqrt(se_mu^2 + sig_x^2))
        hi <- mu + z * (if (band=="ci") se_mu else sqrt(se_mu^2 + sig_x^2))
        lo <- pmax(pmin(lo,hi),0); hi <- pmin(pmax(lo,hi),100)

        if (clamp_bands) {
          lo <- pmax(pmin(lo, hi), 0)
          hi <- pmin(pmax(lo, hi), 100)
        }

        tibble::tibble(Concentration=grid, Viability=mu, Lower=lo, Upper=hi)
      })
    ) %>%
    dplyr::select(-data, -fit) %>%   # qualify to avoid MASS::select masking
    tidyr::unnest(pred)
}

## Plot helpers
plot_across <- function(preds, viab, cmpd){
  df_line <- dplyr::filter(preds, CompoundID==cmpd)
  df_pts  <- dplyr::filter(viab,  CompoundID==cmpd)
  ggplot() +
    geom_ribbon(data=df_line,
                aes(Concentration, ymin=Lower, ymax=Upper,
                    fill=interaction(Patient, Group)),
                alpha=.25, colour=NA, show.legend=FALSE) +
    geom_point(data=df_pts,
               aes(Concentration, Viability,
                   colour=interaction(Patient, Group)),
               size=1.5, alpha=.6) +
    geom_line(data=df_line,
              aes(Concentration, Viability,
                  colour=interaction(Patient, Group)),
              linewidth=.8) +
    scale_x_log10() +
    coord_cartesian(ylim = c(0, 100)) +   # ← clamp view only
    theme_minimal() +
    labs(title  = sprintf("%s across patients & conditions", unique(df_line$CompoundName)),
         x      = sprintf("Conc [%s]", unique(df_line$Unit)),
         y      = "Viability (%)",
         colour = "Patient / Group")
}

plot_within <- function(preds, viab, cmpd){
  df_line <- dplyr::filter(preds, CompoundID==cmpd)
  df_pts  <- dplyr::filter(viab,  CompoundID==cmpd)
  ggplot() +
    geom_ribbon(data=df_line,
                aes(Concentration, ymin=Lower, ymax=Upper, fill=Group),
                alpha=.25, colour=NA, show.legend=FALSE) +
    geom_point(data=df_pts,
               aes(Concentration, Viability, colour=Group),
               size=1.5, alpha=.6) +
    geom_line(data=df_line,
              aes(Concentration, Viability, colour=Group),
              linewidth=.8) +
    scale_x_log10() +
    coord_cartesian(ylim = c(0, 100)) +   # ← clamp view only
    facet_wrap(~ Patient) +
    theme_minimal() +
    labs(title = sprintf("%s by patient conditions", unique(df_line$CompoundName)),
         x     = sprintf("Conc [%s]", unique(df_line$Unit)),
         y     = "Viability (%)")
}

## Run pipeline
layout_long  <- build_layout(compound_path)
all_lines <- read_lines(raw_path)
rows <- which(str_detect(all_lines, regex("_fluo", ignore_case=TRUE)))
rows <- rows[rows + 8 <= length(all_lines)]
if (length(rows)==0) stop("No \"_Fluo\" blocks found.")
all_viability <- get_viability_data(rows, all_lines, layout_long, ref_map)
all_ic50 <- compute_ic50(all_viability)

## Bands (fast)
all_preds_ci <- predict_curves(all_viability, band="ci", level=.95, param_B=3000)
all_preds_pi <- predict_curves(all_viability, band="pi", level=.95, param_B=3000)

## Examples

print(plot_within(all_preds_pi, all_viability, "A1"))
print(plot_within(all_preds_pi, all_viability, "A2"))
print(plot_within(all_preds_pi, all_viability, "A3"))
print(plot_within(all_preds_pi, all_viability, "A4"))
print(plot_within(all_preds_pi, all_viability, "B1"))
print(plot_within(all_preds_pi, all_viability, "B2"))
print(plot_within(all_preds_pi, all_viability, "B3"))
print(plot_within(all_preds_pi, all_viability, "B4"))

print(plot_within(all_preds_ci, all_viability, "A1"))
print(plot_within(all_preds_ci, all_viability, "A2"))
print(plot_within(all_preds_ci, all_viability, "A3"))
print(plot_within(all_preds_ci, all_viability, "A4"))
print(plot_within(all_preds_ci, all_viability, "B1"))
print(plot_within(all_preds_ci, all_viability, "B2"))
print(plot_within(all_preds_ci, all_viability, "B3"))
print(plot_within(all_preds_ci, all_viability, "B4"))