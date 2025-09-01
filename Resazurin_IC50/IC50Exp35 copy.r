# Resazurin Exp35 — per-sample dose–response with drc IC50s + combined plots
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
raw_path    <- "Resazurin_IC50/Resazurin_Exp_36_Raw - Copy.txt"
cmp_path    <- "Resazurin_IC50/Resazurin_Exp_36_Compound Copy.txt"
sample_path <- "Resazurin_IC50/Resazurin_Exp_36_Sample.txt"
out_dir     <- "Resazurin_IC50/resazurin_outputs/Exp36"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "plots"), showWarnings = FALSE)
dir.create(file.path(out_dir, "csv"),   showWarnings = FALSE)

# Helper: parse RAW (Exp35: A1..A4/B1..B4 headers; each of 8 lines has 12 Fluo)
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
      parts   <- stringr::str_split(ln, "\\t")[[1]]  # keep empties
      lvl_raw <- suppressWarnings(as.integer(parts[1]))
      lvl     <- if (!is.na(lvl_raw) && lvl_raw %in% 1:8) lvl_raw else idx
      if (length(parts) < 13) stop("Raw data line has fewer than 13 tab fields:\n", ln)
      fl <- suppressWarnings(as.numeric(parts[2:13]))
      tibble(plate, pos, block, level = lvl, col = 1:12, meas = "Fluo", value = fl)
    })
  })
}

# Map column (1..12) -> sample_id from the Sample file --------------------
parse_sample_map <- function(path) {
  smp <- readr::read_tsv(path, col_types = readr::cols(.default = readr::col_character()))
  if ("ID" %in% names(smp)) {
    grp_cols <- setdiff(names(smp), "ID")
  } else {
    smp <- smp %>% dplyr::mutate(ID = dplyr::row_number())
    grp_cols <- setdiff(names(smp), "ID")
  }

  sample_row <- smp %>% dplyr::filter(tolower(ID) == "sample")
  if (nrow(sample_row) == 1) {
    sample_names <- sample_row %>%
      tidyr::pivot_longer(dplyr::all_of(grp_cols), names_to = "grp", values_to = "sample_id") %>%
      dplyr::select(grp, sample_id)
  } else {
    sample_names <- tibble::tibble(grp = grp_cols, sample_id = grp_cols)
  }

  flag_rows <- smp %>%
    dplyr::filter(grepl("^[0-9]+$", ID)) %>%
    dplyr::mutate(col = as.integer(ID)) %>%
    tidyr::pivot_longer(dplyr::all_of(grp_cols), names_to = "grp", values_to = "flag_chr") %>%
    dplyr::mutate(
      flag_chr = stringr::str_trim(flag_chr),
      flag_num = suppressWarnings(as.numeric(flag_chr)),
      flag = dplyr::case_when(
        is.na(flag_chr) ~ FALSE,
        !is.na(flag_num) ~ flag_num != 0,
        stringr::str_detect(stringr::str_to_lower(flag_chr), "^(true|t|x|yes|y)$") ~ TRUE,
        stringr::str_detect(stringr::str_to_lower(flag_chr), "^(false|f|no|n)$") ~ FALSE,
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

# Read + join layout -------------------------------------------------------
raw_long <- parse_raw_blocks(raw_path)
col_map  <- parse_sample_map(sample_path)

raw_long <- raw_long %>% dplyr::left_join(col_map, by = "col")
stopifnot("sample_id" %in% names(raw_long))
stopifnot(!any(is.na(raw_long$sample_id)))

# Compound sheet -> dose mapping ------------------------------------------
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

dat <- raw_long %>%
  dplyr::mutate(rep = ((col - 1L) %% 3L) + 1L) %>%
  dplyr::left_join(dils, by = c("pos","level","plate","block"))

if (any(!is.finite(dat$dose[dat$meas == "Fluo"]))) {
  bad <- dat %>% dplyr::filter(meas == "Fluo", !is.finite(dose)) %>% dplyr::count(pos, level)
  print(bad, n = 50)
  stop("Dose mapping failed for some (pos, level). Check Compound.txt positions/levels.")
}

# Fluorescence only --------------------------------------------------------
fluo <- dat %>%
  dplyr::filter(meas == "Fluo") %>%
  dplyr::select(sample_id, plate, pos, block, rep, level, dose, Unit, Compound, value)

# Robust top/bottom per curve ---------------------------------------------
top_bottom_map <- fluo %>%
  dplyr::group_by(sample_id, pos) %>%
  dplyr::summarise(
    top = {
      md <- suppressWarnings(min(dose, na.rm = TRUE))
      t1 <- mean(value[dose == md], na.rm = TRUE)
      if (!is.finite(t1)) suppressWarnings(max(value, na.rm = TRUE)) else t1
    },
    bottom = {
      b1 <- suppressWarnings(min(value[level >= 2], na.rm = TRUE))
      if (!is.finite(b1)) suppressWarnings(min(value, na.rm = TRUE)) else b1
    },
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    top    = dplyr::coalesce(ifelse(is.finite(top),    top,    NA_real_), bottom),
    bottom = dplyr::coalesce(ifelse(is.finite(bottom), bottom, NA_real_), top)
  )

# Normalised viability in [0,1] -------------------------------------------
fluo_n <- fluo %>%
  dplyr::left_join(top_bottom_map, by = c("sample_id","pos")) %>%
  dplyr::mutate(
    denom = pmax(top - bottom, .Machine$double.eps, na.rm = TRUE),
    viab  = (value - bottom) / denom,
    viab  = pmin(pmax(viab, 1e-6), 1 - 1e-6)
  )

# Fitting + prediction helpers --------------------------------------------
predict_ll4_ci <- function(fit, new_dose) {
  cf <- coef(fit); vc <- vcov(fit)
  b  <- unname(cf["b:(Intercept)"]); e <- unname(cf["e:(Intercept)"])
  if (is.na(b)) b <- unname(cf[grep("^b", names(cf))[1]])
  if (is.na(e)) e <- unname(cf[grep("^e", names(cf))[1]])
  x <- new_dose; t <- b * (log(x) - log(e)); et <- exp(t); den <- (1 + et); f <- 1/den
  K <- et/(den^2); ddb <- -K * (log(x) - log(e)); dde <- K * (b / e)
  G <- cbind(ddb, dde); vars <- rowSums((G %*% vc) * G); se <- sqrt(pmax(vars, 0))
  tibble(dose = x, fit = f, lwr = pmax(f - 1.96*se, 0), upr = pmin(f + 1.96*se, 1))
}

fit_one_curve <- function(df, fix_bounds = TRUE, B_boot = 600L) {
  dfit <- df %>% dplyr::filter(level >= 2, is.finite(dose), dose > 0, is.finite(viab))
  if (nrow(dfit) < 4L || dplyr::n_distinct(dfit$dose) < 3L) return(NULL)

  fct <- if (fix_bounds) LL.4(fixed = c(NA, 0, 1, NA)) else LL.4()
  fit <- try(drm(viab ~ dose, data = dfit, fct = fct, robust = "mean"), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)

  cf <- coef(fit)
  en <- grep("^e", names(cf), value = TRUE)[1]
  e_val <- as.numeric(cf[en]); se_val <- NA_real_; lower <- NA_real_; upper <- NA_real_

  vc <- try(vcov(fit), silent = TRUE)
  if (!inherits(vc, "try-error") && !is.null(dim(vc)) && is.finite(e_val) && en %in% rownames(vc)) {
    vv <- vc[en, en]; if (is.finite(vv) && vv >= 0) se_val <- sqrt(vv)
  }
  if (!is.finite(se_val)) {
    sm <- try(summary(fit), silent = TRUE)
    if (!inherits(sm, "try-error") && !is.null(sm$coefficients)) {
      rn <- rownames(sm$coefficients); idx <- grep("^e", rn)[1]
      if (!is.na(idx)) {
        se_candidate <- as.numeric(sm$coefficients[idx, "Std. Error"])
        if (is.finite(se_candidate)) se_val <- se_candidate
      }
    }
  }
  if (!is.finite(se_val)) {
    ed_fls <- try(ED(fit, respLev = 0.5, type = "absolute", interval = "fls"), silent = TRUE)
    if (!inherits(ed_fls, "try-error")) {
      df2 <- as.data.frame(ed_fls); cn <- tolower(colnames(df2))
      grab <- function(name) {
        j <- which(cn == tolower(name) | grepl(paste0("^", tolower(name)), cn))
        if (length(j)) as.numeric(df2[1, j[1]]) else NA_real_
      }
      est2 <- grab("estimate"); lower2 <- grab("lower"); upper2 <- grab("upper")
      if (!is.finite(e_val) && is.finite(est2)) e_val <- est2
      if (isTRUE(all(is.finite(c(lower2, upper2))))) {
        lower <- lower2; upper <- upper2; se_val <- (upper2 - lower2) / (2 * 1.96)
      }
    }
  }
  if (!is.finite(se_val)) {
    set.seed(1234); n <- nrow(dfit)
    boot_vals <- replicate(B_boot, {
      idx <- sample.int(n, n, replace = TRUE)
      dd <- dfit[idx, , drop = FALSE]
      ft <- try(drm(viab ~ dose, data = dd, fct = fct, robust = "mean"), silent = TRUE)
      if (inherits(ft, "try-error")) return(NA_real_)
      cfe <- try(coef(ft), silent = TRUE); if (inherits(cfe, "try-error")) return(NA_real_)
      en2 <- grep("^e", names(cfe), value = TRUE)[1]; as.numeric(cfe[en2])
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
    Lower = as.numeric(if (is.finite(lower)) lower else if (is.finite(se_val)) e_val - 1.96*se_val else NA_real_),
    Upper = as.numeric(if (is.finite(upper)) upper else if (is.finite(se_val)) e_val + 1.96*se_val else NA_real_)
  )

  rng  <- range(dfit$dose, na.rm = TRUE)
  grid <- tibble(dose = 10^seq(log10(rng[1]) - 0.05, log10(rng[2]) + 0.05, length.out = 200))
  pred <- try(predict_ll4_ci(fit, grid$dose), silent = TRUE)
  if (inherits(pred, "try-error")) {
    bpar <- as.numeric(cf[grep("^b", names(cf))[1]]); epar <- as.numeric(cf[grep("^e", names(cf))[1]])
    f    <- function(x) 1 / (1 + exp(bpar * (log(x) - log(epar))))
    pred <- tibble(dose = grid$dose, fit = f(grid$dose), lwr = NA_real_, upr = NA_real_)
  }

  p_trend <- tryCatch({
    sm <- summary(fit); cm <- if (!is.null(sm$coefficients)) sm$coefficients else coef(sm)
    rn <- rownames(cm); colp <- grep("Pr", colnames(cm), value = TRUE)[1]; idx <- grep("^b", rn)[1]
    as.numeric(cm[idx, colp])
  }, error = function(e) NA_real_)

  list(fit = fit, ic50 = ic50_df, pred = pred,
       points = dfit %>% dplyr::select(dose, viab, level, rep),
       p_trend = p_trend)
}

# IC50 ratio (delta method on log-IC50) -----------------------------------
ic50_compare <- function(e1, se1, e2, se2) {
  ok <- is.finite(e1) && is.finite(se1) && is.finite(e2) && is.finite(se2) && (e1 > 0) && (e2 > 0)
  if (!isTRUE(ok)) return(tibble(ratio=NA_real_, ratio_lwr=NA_real_, ratio_upr=NA_real_, p_value=NA_real_))
  l1 <- log(e1); l2 <- log(e2); s1l <- se1 / e1; s2l <- se2 / e2
  diff <- l1 - l2; se <- sqrt(s1l^2 + s2l^2); z <- diff / se
  tibble(ratio = exp(diff), ratio_lwr = exp(diff - 1.96*se), ratio_upr = exp(diff + 1.96*se),
         p_value = 2 * pnorm(abs(z), lower.tail = FALSE))
}

# Vectorised formatters (safe in mutate across vectors) -------------------
pretty_p <- function(p) {
  out <- character(length(p))
  bad <- !is.finite(p)
  out[bad] <- "p = n/a"
  tiny <- !bad & p < 1e-4
  out[tiny] <- paste0("p = ", format(p[tiny], digits = 2, scientific = TRUE))
  mid <- !bad & !tiny
  out[mid] <- paste0("p = ", signif(p[mid], 3))
  out
}
fmt_num <- function(x) {
  out <- character(length(x))
  bad <- !is.finite(x)
  out[bad] <- "n/a"
  big <- !bad & x >= 0.01
  out[big] <- formatC(x[big], format = "f", digits = 3)
  sml <- !bad & !big
  out[sml] <- format(x[sml], digits = 2, scientific = TRUE)
  out
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

# IC50 table ---------------------------------------------------------------
ic50_tbl <- res_list %>%
  dplyr::transmute(sample_id, Compound, Unit, pos, fit_ok, p_trend, ic50) %>%
  tidyr::unnest(ic50) %>%
  dplyr::mutate(across(c(ED, SE, Lower, Upper), ~ suppressWarnings(as.numeric(.)))) %>%
  dplyr::arrange(sample_id, Compound)

readr::write_csv(ic50_tbl, file.path(out_dir, "csv", "IC50_results.csv"))

# ===== Per-sample bundle plots: base + all co-treatments ==================
dir.create(file.path(out_dir, "plots_combined", "bundles"), recursive = TRUE, showWarnings = FALSE)

comp_base_of <- function(x) sub("_.+$", "", x)
is_baseline  <- function(x) !grepl("_", x)
base_blue <- "#0072B2"
cb_others <- c("#E69F00", "#009E73", "#CC79A7", "#D55E00", "#F0E442", "#000000", "#56B4E9")

fits <- res_list %>%
  dplyr::filter(fit_ok) %>%
  dplyr::group_by(sample_id, Compound) %>% dplyr::slice(1) %>% dplyr::ungroup() %>%
  dplyr::mutate(comp_base = comp_base_of(Compound), is_base = is_baseline(Compound))

bundles <- fits %>%
  dplyr::group_by(sample_id, comp_base) %>%
  dplyr::summarise(has_base = any(is_base), n_variants = sum(!is_base), .groups = "drop") %>%
  dplyr::filter(has_base, n_variants >= 1)

for (i in seq_len(nrow(bundles))) {
  sid <- bundles$sample_id[i]; cb <- bundles$comp_base[i]
  sub <- fits %>% dplyr::filter(sample_id == sid, comp_base == cb)

  preds <- purrr::map_dfr(seq_len(nrow(sub)), function(j) sub$pred[[j]] %>% dplyr::mutate(group = sub$Compound[[j]]))
  pts   <- purrr::map_dfr(seq_len(nrow(sub)), function(j) sub$points[[j]] %>% dplyr::mutate(group = sub$Compound[[j]]))
  vdat  <- tibble::tibble(dose = purrr::map_dbl(sub$ic50, ~ suppressWarnings(as.numeric(.x[1, "ED"]))),
                          group = sub$Compound) %>% dplyr::filter(is.finite(dose))

  groups <- unique(c(preds$group, pts$group))
  cols   <- setNames(rep_len(cb_others, length(groups)), groups)
  if (cb %in% names(cols)) cols[cb] <- base_blue

  unit_lbl <- sub$Unit[which(!is.na(sub$Unit))[1]]
  ttl <- glue::glue("{sid} — {cb} (base) ± pre-treatments")

  p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = preds, ggplot2::aes(dose, ymin = lwr, ymax = upr, fill = group), alpha = 0.18) +
    ggplot2::geom_line(  data = preds, ggplot2::aes(dose, fit, colour = group), linewidth = 0.9) +
    ggplot2::geom_point( data = pts,   ggplot2::aes(dose, viab, colour = group), size = 2, alpha = 0.85) +
    ggplot2::geom_vline(data = vdat,   ggplot2::aes(xintercept = dose, colour = group),
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

  fn <- file.path(out_dir, "plots_combined", "bundles",
                  glue::glue("bundle_{gsub('[^A-Za-z0-9]+','_', sid)}_{gsub('[^A-Za-z0-9]+','_', cb)}.png"))
  ggplot2::ggsave(fn, p, width = 7.2, height = 5.0, dpi = 300)
}

# ===== Between-sample plots: 397 vs 399; 397cse vs 399cse ==================
dir.create(file.path(out_dir, "plots_between_samples"), recursive = TRUE, showWarnings = FALSE)

res_fit <- res_list %>%
  dplyr::filter(fit_ok) %>%
  dplyr::group_by(sample_id, Compound) %>% dplyr::slice(1) %>% dplyr::ungroup()

make_pair_plots <- function(sid_a, sid_b, tag) {
  if (!(sid_a %in% res_fit$sample_id) || !(sid_b %in% res_fit$sample_id)) {
    message("Skipping ", tag, ": missing sample(s)."); return(invisible(NULL))
  }
  sub_a <- res_fit %>% dplyr::filter(sample_id == sid_a)
  sub_b <- res_fit %>% dplyr::filter(sample_id == sid_b)
  comps <- intersect(unique(sub_a$Compound), unique(sub_b$Compound))
  if (!length(comps)) { message("No common compounds for ", tag); return(invisible(NULL)) }

  # Fixed colours + labels: 397(±cse)=blue (current), 399(±cse)=red (former)
  pal_rule <- function(ids) {
    base <- c("ASC397"="#0072B2","ASC397cse"="#0072B2","ASC399"="#D55E00","ASC399cse"="#D55E00")
    base[ids]
  }
  status_lbl <- c(
    "ASC397"    = "ASC397 (current)",
    "ASC397cse" = "ASC397cse (current)",
    "ASC399"    = "ASC399 (former)",
    "ASC399cse" = "ASC399cse (former)"
  )
  pal <- pal_rule(c(sid_a, sid_b))

  out_sub <- file.path(out_dir, "plots_between_samples", tag)
  dir.create(out_sub, recursive = TRUE, showWarnings = FALSE)

  for (cmp_name in comps) {
    da <- sub_a %>% dplyr::filter(Compound == cmp_name) %>% dplyr::slice(1)
    db <- sub_b %>% dplyr::filter(Compound == cmp_name) %>% dplyr::slice(1)
    if (!nrow(da) || !nrow(db)) next

    preds <- dplyr::bind_rows(
      da$pred[[1]] %>% dplyr::mutate(group = sid_a),
      db$pred[[1]] %>% dplyr::mutate(group = sid_b)
    )
    pts <- dplyr::bind_rows(
      da$points[[1]] %>% dplyr::mutate(group = sid_a),
      db$points[[1]] %>% dplyr::mutate(group = sid_b)
    )
    vdat <- tibble::tibble(
      dose  = c(da$ic50[[1]][1, "ED", drop = TRUE], db$ic50[[1]][1, "ED", drop = TRUE]),
      group = c(sid_a, sid_b)
    ) %>% dplyr::filter(is.finite(dose))

    unit_lbl <- dplyr::coalesce(da$Unit[[1]], db$Unit[[1]])

    p <- ggplot2::ggplot() +
      ggplot2::geom_ribbon(data = preds, ggplot2::aes(dose, ymin = lwr, ymax = upr, fill = group),
                           alpha = 0.18, show.legend = FALSE) +
      ggplot2::geom_line(data = preds, ggplot2::aes(dose, fit, colour = group), linewidth = 0.9) +
      ggplot2::geom_point(data = pts, ggplot2::aes(dose, viab, colour = group), size = 2, alpha = 0.85) +
      ggplot2::geom_vline(data = vdat, ggplot2::aes(xintercept = dose, colour = group),
                          linetype = 2, linewidth = 0.7, alpha = 0.9) +
      ggplot2::scale_x_log10(labels = scales::label_number(),
                              expand = ggplot2::expansion(mult = c(0.02, 0.05))) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::scale_colour_manual(
        values = pal,
        breaks = c(sid_a, sid_b),
        labels = status_lbl[c(sid_a, sid_b)],
        name   = "Sample (smoking status)"
      ) +
      ggplot2::scale_fill_manual(values = scales::alpha(pal, 0.18), guide = "none") +
      ggplot2::labs(title = glue::glue("{cmp_name} — {sid_a} vs {sid_b}"),
                    x = glue::glue("Concentration ({unit_lbl})"),
                    y = "Normalised viability") +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                     panel.grid.minor.x = ggplot2::element_blank())

    fn <- file.path(out_sub, glue::glue(
      "between_{gsub('[^A-Za-z0-9]+','_', cmp_name)}_{gsub('[^A-Za-z0-9]+','_', sid_a)}_vs_{gsub('[^A-Za-z0-9]+','_', sid_b)}.png"
    ))
    ggplot2::ggsave(fn, p, width = 7.2, height = 5.0, dpi = 300)
  }
}

# Two requested comparisons; typically 8 plots if 4 compounds exist
make_pair_plots("ASC397",    "ASC399",    "ASC397_vs_ASC399")
make_pair_plots("ASC397cse", "ASC399cse", "ASC397cse_vs_ASC399cse")

# ===== IC50 stats: sample vs sample (per compound/co-treatment) ==========
between_sample_test <- function(idA, idB, out_stub = NULL) {
  d <- ic50_tbl %>%
    dplyr::mutate(sample_id = as.character(sample_id)) %>%
    dplyr::distinct(sample_id, Compound, Unit, .keep_all = TRUE) %>%
    dplyr::filter(sample_id %in% c(idA, idB)) %>%
    dplyr::mutate(ED = suppressWarnings(as.numeric(ED)),
                  SE = suppressWarnings(as.numeric(SE)))

  A <- d %>% dplyr::filter(sample_id == idA) %>% dplyr::select(Compound, Unit, ED_A = ED, SE_A = SE)
  B <- d %>% dplyr::filter(sample_id == idB) %>% dplyr::select(Compound, Unit, ED_B = ED, SE_B = SE)

  j <- dplyr::inner_join(A, B, by = c("Compound", "Unit"))

  res <- j %>%
    dplyr::rowwise() %>%
    dplyr::mutate(stats = ic50_compare(ED_A, SE_A, ED_B, SE_B)) %>%
    tidyr::unnest_wider(stats) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      sample_A = idA,
      sample_B = idB,
      ratio_desc = paste0(idA, "/", idB),
      summary = paste0(
        "ratio=", fmt_num(ratio),
        "  [", fmt_num(ratio_lwr), ", ", fmt_num(ratio_upr), "]  ",
        pretty_p(p_value)
      )
    ) %>%
    dplyr::select(sample_A, sample_B, Compound, Unit,
                  ED_A, SE_A, ED_B, SE_B,
                  ratio, ratio_lwr, ratio_upr, p_value, summary) %>%
    dplyr::arrange(Compound)

  if (!is.null(out_stub)) {
    readr::write_csv(res, file.path(out_dir, "csv", paste0("IC50_compare_", out_stub, ".csv")))
  }
  res
}

tab_397_vs_399        <- between_sample_test("ASC397",    "ASC399",    "ASC397_vs_ASC399")
tab_397cse_vs_399cse  <- between_sample_test("ASC397cse", "ASC399cse", "ASC397cse_vs_ASC399cse")

message("Done. IC50s: ", file.path(out_dir, "csv", "IC50_results.csv"))
message("Bundle plots in: ", file.path(out_dir, "plots_combined", "bundles"))
message("Between-sample plots in: ", file.path(out_dir, "plots_between_samples"))
message("Between-sample IC50 CSVs in: ", file.path(out_dir, "csv"))
