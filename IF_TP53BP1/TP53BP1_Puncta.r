# TP53BP1 — Puncta per nucleus-with-puncta only
# (total puncta) / (nuclei with puncta); robust to column-name variants
# Experiments discovered under: IF_TP53BP1/
# Outputs -> figures/, results/

suppressPackageStartupMessages({
  need <- c("tidyverse","readr","readxl","janitor","glue","emmeans","MASS","rlang","scales","svglite")
  to_install <- setdiff(need, rownames(installed.packages()))
  if (length(to_install)) install.packages(to_install, dep = TRUE)
  invisible(lapply(need, library, character.only = TRUE))
})

# ---------- Paths ----------
root <- "IF_TP53BP1"   # change to absolute path if needed
stopifnot(dir.exists(root))
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# ============================== Helpers ======================================
# ---- NEW: label recoder for CC/CX/XC/XX -> CTL-CTL etc. ----
recode_cc_labels <- function(x) {
  x <- as.character(x)
  dplyr::recode(x,
    "CC" = "CTL-CTL",
    "CX" = "CTL-CSE",
    "XC" = "CSE-CTL",
    "XX" = "CSE-CSE",
    .default = x
  )
}

.list_all <- function(root) {
  exts <- c("xlsx","xls","csv","tsv","txt")
  pats <- paste0("\\.(", paste(exts, collapse="|"), ")$")
  all <- list.files(root, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  all[grepl(pats, all, ignore.case = TRUE)]
}
find_one_any_strict <- function(stem_re, root = ".") {
  pat  <- paste0(stem_re, "\\.(xlsx|xls|csv|tsv|txt)$")
  hits <- list.files(root, pattern = pat, recursive = TRUE,
                     full.names = TRUE, ignore.case = TRUE)
  if (!length(hits)) stop(sprintf("STRICT: no match for /^%s\\.(xlsx|xls|csv|tsv|txt)$/ under: %s", stem_re, root))
  normalizePath(hits[1], mustWork = TRUE)
}
find_one_tokens <- function(tokens, root = ".") {
  all <- .list_all(root); bn <- basename(all)
  keep <- all[vapply(seq_along(all), function(i)
    all(sapply(tokens, function(t) grepl(t, bn[i], ignore.case = TRUE))), logical(1))]
  if (!length(keep)) stop(sprintf("TOKENS: no file with tokens: %s", paste(tokens, collapse = " + ")))
  depth    <- lengths(gregexpr("[/\\\\]", keep))
  ext_rank <- match(tolower(tools::file_ext(keep)), c("xlsx","xls","csv","tsv","txt"))
  keep <- keep[order(depth, ext_rank)]
  normalizePath(keep[1], mustWork = TRUE)
}
safe_find <- function(strict_re, tokens, root = ".") {
  tryCatch(
    find_one_any_strict(strict_re, root),
    error = function(e) {
      message(e$message)
      message("Falling back to tokens: [", paste(tokens, collapse = " + "), "]")
      find_one_tokens(tokens, root)
    }
  )
}

# Prefer non-senescence Exp70; if not found, optionally fall back to ANY Exp70
find_exp70_any <- function(root, allow_fallback_any = TRUE) {
  all <- .list_all(root)
  # Exclude saved outputs
  all <- all[!grepl("(^|[/\\\\])(figures|results)($|[/\\\\])", all, ignore.case = TRUE)]
  bn  <- basename(all)

  non_sen_mask <- grepl("(?i)exp\\s*[-_ ]?\\s*70", bn) & !grepl("(?i)sen", bn)
  cand1 <- all[non_sen_mask]
  if (length(cand1)) {
    cand1 <- cand1[order(lengths(gregexpr("[/\\\\]", cand1)),
                         match(tolower(tools::file_ext(cand1)), c("xlsx","xls","csv","tsv","txt")))]
    return(normalizePath(cand1[1], mustWork = TRUE))
  }

  if (allow_fallback_any) {
    cand2 <- all[grepl("(?i)exp\\s*[-_ ]?\\s*70", bn)]
    if (length(cand2)) {
      message("NOTE: Using a generic Exp70 file (may include 'Sen').")
      cand2 <- cand2[order(lengths(gregexpr("[/\\\\]", cand2)),
                           match(tolower(tools::file_ext(cand2)), c("xlsx","xls","csv","tsv","txt")))]
      return(normalizePath(cand2[1], mustWork = TRUE))
    }
  }
  stop("No Exp70 file found under: ", root)
}

parse_dose_num <- function(conc_char) {
  out <- suppressWarnings(as.numeric(sub("([0-9.]+).*", "\\1", as.character(conc_char))))
  ifelse(is.na(out), NA_real_, out)
}

# ---------- Universal reader ----------
read_cp <- function(path) {
  ext <- tolower(tools::file_ext(path))
  raw <- switch(ext,
    "xlsx" = readxl::read_excel(path, sheet = 1),
    "xls"  = readxl::read_excel(path, sheet = 1),
    "csv"  = readr::read_csv(path, show_col_types = FALSE),
              readr::read_tsv(path, show_col_types = FALSE)
  )
  raw
}

# ---------- Robust column normaliser ----------
# Standardised names after this function:
#   count_nuclei, count_nuclei_with_dots, count_puncta,
#   metadata_conc, metadata_sample
unify_cp_cols <- function(df) {
  df <- janitor::clean_names(df)
  nm <- names(df)
  pick_col <- function(...) {
    pats <- c(...)
    for (p in pats) {
      hit <- grep(p, nm, ignore.case = TRUE, value = TRUE)
      if (length(hit)) return(hit[1])
    }
    NULL
  }

  # total nuclei
  if (!"count_nuclei" %in% nm) {
    cand <- pick_col("^count_?nuclei($|_)", "nuclei(_count)?$", "^(total_)?cells?$", "cell(_count)?$")
    if (!is.null(cand)) df <- dplyr::rename(df, count_nuclei = !!rlang::sym(cand))
  }

  # nuclei WITH puncta/dots/foci
  if (!"count_nuclei_with_dots" %in% nm) {
    cand <- pick_col("nuclei.*(with_)?(dots|puncta|foci)",
                     "(puncta|foci).*nuclei",
                     "nuclei_?(pos|positive)",
                     "^count_?positive_?nuclei$")
    if (!is.null(cand)) df <- dplyr::rename(df, count_nuclei_with_dots = !!rlang::sym(cand))
  }

  # total puncta
  if (!"count_puncta" %in% nm) {
    cand <- pick_col("^count_?children_?puncta.*", "^count_?puncta($|_)",
                     "total_?puncta", "puncta($|_)")
    if (!is.null(cand)) df <- dplyr::rename(df, count_puncta = !!rlang::sym(cand))
  }

  # concentration / dose / condition
  if (!"metadata_conc" %in% nm) {
    cand <- pick_col("^conc($|_)", "concentration", "^dose($|_)",
                     "^treatment$", "metadata_?(concentration|dose)", "^condition$")
    if (!is.null(cand)) df <- dplyr::rename(df, metadata_conc = !!rlang::sym(cand))
  }

  # sample / line
  if (!"metadata_sample" %in% nm) {
    cand <- pick_col("^sample($|_)", "^well($|_)", "^group($|_)",
                     "metadata_?(well|group|sampleid)", "^line$")
    if (!is.null(cand)) df <- dplyr::rename(df, metadata_sample = !!rlang::sym(cand))
    else df$metadata_sample <- NA_character_
  }

  df
}

# ---------- EMMEANS helpers ----------
std_emm_cols <- function(df) {
  nm <- names(df)
  if ("response"    %in% nm) df <- dplyr::rename(df, estimate = response)
  else if ("prob"   %in% nm) df <- dplyr::rename(df, estimate = prob)
  else if ("rate"   %in% nm) df <- dplyr::rename(df, estimate = rate)
  else if ("ratio"  %in% nm) df <- dplyr::rename(df, estimate = ratio)
  else if ("emmean" %in% nm) df <- dplyr::rename(df, estimate = emmean)
  nm <- names(df); if ("asymp.LCL" %in% nm) df <- dplyr::rename(df, lcl = asymp.LCL)
  else if ("lower.CL" %in% nm) df <- dplyr::rename(df, lcl = lower.CL)
  else if ("LCL" %in% nm) df <- dplyr::rename(df, lcl = LCL)
  nm <- names(df); if ("asymp.UCL" %in% nm) df <- dplyr::rename(df, ucl = asymp.UCL)
  else if ("upper.CL" %in% nm) df <- dplyr::rename(df, ucl = upper.CL)
  else if ("UCL" %in% nm) df <- dplyr::rename(df, ucl = UCL)
  df
}
emm_pairs_vs_ctrl <- function(mod, specs, offset0 = FALSE) {
  if (offset0) {
    emmeans(mod, specs, type = "response", offset = 0) |>
      contrast("trt.vs.ctrl", ref = 1) |>
      summary(infer = TRUE, adjust = "BH") |> as_tibble()
  } else {
    emmeans(mod, specs, type = "response") |>
      contrast("trt.vs.ctrl", ref = 1) |>
      summary(infer = TRUE, adjust = "BH") |> as_tibble()
  }
}

# ---------- Little plotting bits ----------
theme_pub <- theme_minimal(base_size = 12) +
  theme(text = element_text(family = "Arial", size = 12),
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        legend.position = "bottom")

stars_map <- function(p) dplyr::case_when(is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")

add_brackets_to_plot <- function(p, br_df) {
  if (!nrow(br_df)) return(p)
  tick <- 0.012
  p +
    geom_segment(data = br_df, aes(x = x1, xend = x2, y = y, yend = y),
                 inherit.aes = FALSE, linewidth = 0.4) +
    geom_segment(data = br_df, aes(x = x1, xend = x1, y = y, yend = y - tick),
                 inherit.aes = FALSE, linewidth = 0.4) +
    geom_segment(data = br_df, aes(x = x2, xend = x2, y = y, yend = y - tick),
                 inherit.aes = FALSE, linewidth = 0.4) +
    geom_text(   data = br_df, aes(x = xm, y = y + 0.8*tick, label = star),
                 inherit.aes = FALSE, vjust = 0, size = 3.5)
}

# =========================== FILE DISCOVERY ==================================
file_9C_CSE <- safe_find("^9C[_-]?CSE$",   c("9C","CSE"),     root)
file_9X_CSE <- safe_find("^9X[_-]?CSE$",   c("9X","CSE"),     root)
file_9C_KBR <- safe_find("^9C[_-]?KBrO3$", c("9C","KBrO3"),   root)
file_9X_KBR <- safe_find("^9X[_-]?KBrO3$", c("9X","KBrO3"),   root)
file_Exp70  <- tryCatch(find_exp70_any(root, allow_fallback_any = TRUE),
                        error = function(e) { message(e$message); NA_character_ })
has_E70 <- is.character(file_Exp70) && !is.na(file_Exp70) && file.exists(file_Exp70)

cat("Using files:\n - ", file_9C_CSE,
    "\n - ", file_9X_CSE,
    "\n - ", file_9C_KBR,
    "\n - ", file_9X_KBR,
    if (has_E70) paste0("\n - ", file_Exp70) else "\n - (no Exp70 file found)",
    "\n", sep = "")

# =========================== LOAD & TIDY =====================================
as_dat <- function(path) {
  x <- read_cp(path) |> unify_cp_cols()
  nm <- names(x)

  x <- x |>
    mutate(file = basename(path)) |>
    mutate(
      experiment = dplyr::case_when(
        grepl("CSE", path, ignore.case=TRUE) ~ "CSE",
        grepl("KBrO3|KBrO", path, ignore.case=TRUE) ~ "KBrO3",
        grepl("SenExp", path, ignore.case=TRUE) ~ "SenExp",
        TRUE ~ "Exp"
      ),
      line = dplyr::case_when(
        grepl("9C", file, ignore.case=TRUE) ~ "ASC399cse",
        grepl("9X", file, ignore.case=TRUE) ~ "ASC399",
        TRUE ~ as.character(metadata_sample)
      ),
      conc = as.character(metadata_conc),
      dose_num = parse_dose_num(conc)
    )

  # (existing) puncta per nucleus-with-puncta
  if ("count_nuclei_with_dots" %in% names(x)) {
    x <- x |>
      mutate(puncta_per_withdots = dplyr::if_else(
        is.finite(count_nuclei_with_dots) & count_nuclei_with_dots > 0,
        count_puncta / count_nuclei_with_dots,
        NA_real_
      ))
  } else {
    x <- x |> mutate(puncta_per_withdots = NA_real_)
  }

  # ---- NEW: % puncta-positive nuclei (of total nuclei) ----
  if (all(c("count_nuclei_with_dots","count_nuclei") %in% names(x))) {
    x <- x |>
      mutate(pct_pos_nuclei = dplyr::if_else(
        is.finite(count_nuclei) & count_nuclei > 0,
        100 * count_nuclei_with_dots / count_nuclei,
        NA_real_
      ))
  } else {
    x <- x |> mutate(pct_pos_nuclei = NA_real_)
  }

  x
}

dat_CSE <- bind_rows(as_dat(file_9C_CSE), as_dat(file_9X_CSE)) %>%
  mutate(experiment = "CSE") %>% filter(!is.na(line))
dat_KBR <- bind_rows(as_dat(file_9C_KBR), as_dat(file_9X_KBR)) %>%
  mutate(experiment = "KBrO3") %>% filter(!is.na(line))
if (has_E70) {
  dat_E70 <- as_dat(file_Exp70) %>%
    mutate(experiment = "Exp70") %>% filter(!is.na(line))
}

# factor orders for doses where numeric present
levs_cse <- dat_CSE %>% distinct(conc, dose_num) %>% arrange(dose_num) %>% pull(conc)
if (all(is.na(dat_CSE$dose_num))) levs_cse <- sort(unique(dat_CSE$conc))
dat_CSE$conc <- factor(dat_CSE$conc, levels = levs_cse)

levs_kbr <- dat_KBR %>% distinct(conc, dose_num) %>% arrange(dose_num) %>% pull(conc)
if (all(is.na(dat_KBR$dose_num))) levs_kbr <- sort(unique(dat_KBR$conc))
dat_KBR$conc <- factor(dat_KBR$conc, levels = levs_kbr)

if (has_E70) {
  levs_e70 <- dat_E70 %>% distinct(conc, dose_num) %>% arrange(dose_num) %>% pull(conc)
  if (all(is.na(dat_E70$dose_num))) levs_e70 <- sort(unique(dat_E70$conc))
  dat_E70$conc <- factor(dat_E70$conc, levels = levs_e70)
}

# =========================== STATS CORE ======================================
# NB model for counts with offset(log(nuclei with puncta))
fit_nb_rate <- function(d, count_col = "count_puncta") {
  base_col <- if ("count_nuclei_with_dots" %in% names(d)) "count_nuclei_with_dots"
             else if ("count_nuclei" %in% names(d)) "count_nuclei"
             else stop("No nuclei count column found (neither 'count_nuclei_with_dots' nor 'count_nuclei').")

  d2 <- d %>%
    filter(!is.na(.data[[base_col]]), .data[[base_col]] > 0) %>%
    mutate(offset_log = log(.data[[base_col]])) %>%
    filter(is.finite(offset_log), !is.na(.data[[count_col]])) %>%
    droplevels()

  f_int <- as.formula(glue::glue("{count_col} ~ line * conc + offset(offset_log)"))
  f_add <- as.formula(glue::glue("{count_col} ~ line + conc + offset(offset_log)"))
  m <- try(MASS::glm.nb(f_int, data = d2), silent = TRUE)
  if (inherits(m, "try-error")) m <- try(MASS::glm.nb(f_add, data = d2), silent = TRUE)
  if (inherits(m, "try-error")) {
    message("NB failed; using quasi-Poisson.")
    m <- try(glm(f_int, data = d2, family = quasipoisson()), silent = TRUE)
    if (inherits(m, "try-error")) m <- glm(f_add, data = d2, family = quasipoisson())
  }
  m
}

# ================================ PLOTS ======================================
# One plot per experiment: puncta per nucleus-with-puncta (raw jitter + EMMEANS)
run_expt <- function(expt_name, dat_expt, pal_lines) {
  stopifnot(nrow(dat_expt) > 0)

  # recoded x labels (keeps original ordering)
  lv         <- levels(dat_expt$conc)
  lv_lab     <- recode_cc_labels(lv)
  dat_expt   <- dat_expt %>% mutate(conc_lab = factor(recode_cc_labels(conc), levels = lv_lab))

  m_rate   <- fit_nb_rate(dat_expt, "count_puncta")
  emm_rate <- emmeans(m_rate, ~ conc | line, type = "response", offset = 0) |>
              as_tibble() |> std_emm_cols() |>
              mutate(conc_lab = factor(recode_cc_labels(conc), levels = lv_lab))
  pairs_rt <- emm_pairs_vs_ctrl(m_rate, ~ conc | line, offset0 = TRUE)

  raw_pts <- dat_expt %>%
    filter(!is.na(puncta_per_withdots)) %>%
    dplyr::select(line, conc_lab, puncta_per_withdots)

  # Brackets with recoded targets
  pos_map <- setNames(seq_along(lv_lab), lv_lab)
  ci_rate <- emm_rate |> group_by(line, conc_lab) |> summarise(upr = max(ucl), .groups="drop")
  base_rate <- ci_rate |> group_by(line) |> summarise(y0 = max(upr, na.rm = TRUE) + 0.04)

  br_rate <- pairs_rt |>
    mutate(target_raw = sub("\\s*-\\s*.*$", "", contrast),
           target_lab = recode_cc_labels(target_raw),
           star = stars_map(p.value)) |>
    filter(star != "", target_lab %in% lv_lab[-1]) |>
    left_join(base_rate, by = "line") |>
    mutate(x1 = pos_map[lv_lab[1]], x2 = pos_map[target_lab], xm = (x1 + x2)/2,
           order = match(target_lab, lv_lab),
           idx = ave(order, line, FUN = function(z) rank(z, ties.method = "first")),
           y = y0 + (idx - 1) * 0.04) |>
    dplyr::select(line, x1, x2, xm, y, star)

  p <- ggplot() +
    geom_point(data = raw_pts,
               aes(conc_lab, puncta_per_withdots, colour = line),
               alpha = 0.7, size = 1.7,
               position = position_jitter(width = 0.08, height = 0, seed = 1)) +
    geom_errorbar(data = emm_rate, aes(conc_lab, ymin = lcl, ymax = ucl, colour = line),
                  width = 0.08, linewidth = 0.4,
                  position = position_dodge(width = 0.25)) +
    geom_point(data = emm_rate, aes(conc_lab, estimate, colour = line),
               size = 2.4, position = position_dodge(width = 0.25)) +
    scale_colour_manual(values = pal_lines, name = "Sample") +
    labs(x = "Condition", y = "Puncta per nucleus with puncta",
         title = "Puncta per nucleus with puncta") +
    facet_wrap(~ line, nrow = 1, scales = "free_y") + theme_pub

  p <- add_brackets_to_plot(p, br_rate)

  stub <- gsub("\\W+", "_", expt_name)
  ggsave(file.path("figures", paste0(stub, "_puncta_per_NucWithPuncta.png")), p, width = 7.0, height = 3.0, dpi = 600)
  ggsave(file.path("figures", paste0(stub, "_puncta_per_NucWithPuncta.svg")), p, width = 7.0, height = 3.0, device = svglite::svglite)

  # add readable labels to CSV output too
  write_csv(emm_rate %>% relocate(conc_lab, .after = conc),
            file.path("results", paste0(stub, "_puncta_per_NucWithPuncta_emmeans.csv")))
  write_csv(pairs_rt,  file.path("results", paste0(stub, "_puncta_per_NucWithPuncta_pairs_vs_ctrl.csv")))
}
# ---- NEW: Binomial model + plot for puncta-positive nuclei (%) ----
fit_binom_pos <- function(d) {
  stopifnot(all(c("count_nuclei_with_dots","count_nuclei") %in% names(d)))
  d2 <- d %>%
    filter(is.finite(count_nuclei), count_nuclei > 0,
           is.finite(count_nuclei_with_dots), count_nuclei_with_dots >= 0) %>%
    mutate(neg = pmax(count_nuclei - count_nuclei_with_dots, 0L)) %>%
    droplevels()
  # quasi-binomial for safety (overdispersion)
  glm(cbind(count_nuclei_with_dots, neg) ~ line * conc, data = d2, family = quasibinomial())
}

run_expt_pct <- function(expt_name, dat_expt, pal_lines) {
  stopifnot(nrow(dat_expt) > 0)

  lv         <- levels(dat_expt$conc)
  lv_lab     <- recode_cc_labels(lv)
  dat_expt   <- dat_expt %>% mutate(conc_lab = factor(recode_cc_labels(conc), levels = lv_lab))

  m_pos   <- fit_binom_pos(dat_expt)
  emm_pos <- emmeans(m_pos, ~ conc | line, type = "response") |>
             as_tibble() |> std_emm_cols() |>
             mutate(estimate = 100*estimate, lcl = 100*lcl, ucl = 100*ucl,
                    conc_lab = factor(recode_cc_labels(conc), levels = lv_lab))
  pairs_p <- emm_pairs_vs_ctrl(m_pos, ~ conc | line)  # BH-adjusted

  raw_pts <- dat_expt %>%
    filter(!is.na(pct_pos_nuclei)) %>%
    dplyr::select(line, conc_lab, pct_pos_nuclei)

  # Brackets at % scale
  pos_map <- setNames(seq_along(lv_lab), lv_lab)
  ci_pos  <- emm_pos |> group_by(line, conc_lab) |> summarise(upr = max(ucl), .groups="drop")
  base_y  <- ci_pos |> group_by(line) |> summarise(y0 = max(upr, na.rm = TRUE) + 1.5)

  br_pos <- pairs_p |>
    mutate(target_raw = sub("\\s*-\\s*.*$", "", contrast),
           target_lab = recode_cc_labels(target_raw),
           star = stars_map(p.value)) |>
    filter(star != "", target_lab %in% lv_lab[-1]) |>
    left_join(base_y, by = "line") |>
    mutate(x1 = pos_map[lv_lab[1]], x2 = pos_map[target_lab], xm = (x1 + x2)/2,
           order = match(target_lab, lv_lab),
           idx = ave(order, line, FUN = function(z) rank(z, ties.method = "first")),
           y = y0 + (idx - 1) * 1.5) |>
    dplyr::select(line, x1, x2, xm, y, star)

  p <- ggplot() +
    geom_point(data = raw_pts,
               aes(conc_lab, pct_pos_nuclei, colour = line),
               alpha = 0.7, size = 1.7,
               position = position_jitter(width = 0.08, height = 0, seed = 1)) +
    geom_errorbar(data = emm_pos, aes(conc_lab, ymin = lcl, ymax = ucl, colour = line),
                  width = 0.08, linewidth = 0.4,
                  position = position_dodge(width = 0.25)) +
    geom_point(data = emm_pos, aes(conc_lab, estimate, colour = line),
               size = 2.4, position = position_dodge(width = 0.25)) +
    scale_colour_manual(values = pal_lines, name = "Sample") +
    labs(x = "Condition", y = "Puncta-positive nuclei (%)",
         title = "Puncta-positive nuclei (%)") +
    facet_wrap(~ line, nrow = 1) + theme_pub

  p <- add_brackets_to_plot(p, br_pos)

  stub <- gsub("\\W+", "_", expt_name)
  ggsave(file.path("figures", paste0(stub, "_PctPosNuclei.png")), p, width = 7.0, height = 3.0, dpi = 600)
  ggsave(file.path("figures", paste0(stub, "_PctPosNuclei.svg")), p, width = 7.0, height = 3.0, device = svglite::svglite)

  write_csv(emm_pos %>% relocate(conc_lab, .after = conc),
            file.path("results", paste0(stub, "_PctPosNuclei_emmeans.csv")))
  write_csv(pairs_p, file.path("results", paste0(stub, "_PctPosNuclei_pairs_vs_ctrl.csv")))
}

# Palettes per experiment
pal_CSE <- setNames(c("#0072B2","#E69F00"), sort(unique(dat_CSE$line)))
pal_KBR <- setNames(c("#0072B2","#E69F00"), sort(unique(dat_KBR$line)))

# Run
# Existing
run_expt("CSE",  dat_CSE, pal_CSE)
run_expt("KBrO3",dat_KBR, pal_KBR)

# NEW: percentage-positive outputs
run_expt_pct("CSE",  dat_CSE, pal_CSE)
run_expt_pct("KBrO3",dat_KBR, pal_KBR)

# ---------- Exp70 (same metric) ----------
if (has_E70) {
  message("Running Exp70 section on: ", file_Exp70)

  # --- model data
  d_nb_e70 <- dat_E70 %>%
    filter(!is.na(count_nuclei_with_dots), count_nuclei_with_dots > 0) %>%
    mutate(offset_log = log(count_nuclei_with_dots)) %>%
    filter(is.finite(offset_log), !is.na(count_puncta)) %>% droplevels()

  f_int_e70 <- as.formula("count_puncta ~ line * conc + offset(offset_log)")
  f_add_e70 <- as.formula("count_puncta ~ line + conc + offset(offset_log)")
  m_e70 <- try(MASS::glm.nb(f_int_e70, data = d_nb_e70), silent = TRUE)
  if (inherits(m_e70, "try-error")) m_e70 <- try(MASS::glm.nb(f_add_e70, data = d_nb_e70), silent = TRUE)
  if (inherits(m_e70, "try-error")) {
    message("Exp70: NB failed; using quasi-Poisson.")
    m_e70 <- try(glm(f_int_e70, data = d_nb_e70, family = quasipoisson()), silent = TRUE)
    if (inherits(m_e70, "try-error")) m_e70 <- glm(f_add_e70, data = d_nb_e70, family = quasipoisson())
  }

  # --- recode CC/CX/XC/XX -> CTL-CTL / ...
  lv_e70      <- levels(dat_E70$conc)
  lv_e70_lab  <- recode_cc_labels(lv_e70)
  dat_E70     <- dat_E70 %>% mutate(conc_lab = factor(recode_cc_labels(conc), levels = lv_e70_lab))

  emm_e70  <- emmeans(m_e70, ~ conc | line, type = "response", offset = 0) %>%
              as_tibble() %>% std_emm_cols() %>%
              mutate(conc_lab = factor(recode_cc_labels(conc), levels = lv_e70_lab))
  pairs_e70<- emm_pairs_vs_ctrl(m_e70, ~ conc | line, offset0 = TRUE)

  pal_E70 <- setNames(scales::hue_pal()(length(sort(unique(dat_E70$line)))), sort(unique(dat_E70$line)))

  p_e70 <- ggplot() +
    geom_point(data = dat_E70 %>% filter(!is.na(puncta_per_withdots)),
               aes(conc_lab, puncta_per_withdots, colour = line),
               alpha = 0.7, size = 1.7,
               position = position_jitter(width = 0.08, height = 0, seed = 1)) +
    geom_errorbar(data = emm_e70, aes(conc_lab, ymin = lcl, ymax = ucl, colour = line),
                  width = 0.08, linewidth = 0.4,
                  position = position_dodge(width = 0.25)) +
    geom_point(data = emm_e70, aes(conc_lab, estimate, colour = line),
               size = 2.4, position = position_dodge(width = 0.25)) +
    scale_colour_manual(values = pal_E70, name = "Sample") +
    labs(x = "Condition", y = "Puncta per nucleus with puncta",
         title = "Puncta per nucleus with puncta") +
    facet_wrap(~ line, nrow = 1, scales = "free_y") + theme_pub

  # --- brackets with recoded labels
  ci_e70   <- emm_e70 %>% group_by(line, conc_lab) %>% summarise(upr = max(ucl), .groups="drop")
  base_y   <- ci_e70 %>% group_by(line) %>% summarise(y0 = max(upr, na.rm = TRUE) + 0.04)
  pos_map_e70 <- setNames(seq_along(lv_e70_lab), lv_e70_lab)

  br_e70 <- pairs_e70 %>%
    mutate(target_raw = sub("\\s*-\\s*.*$", "", contrast),
           target_lab = recode_cc_labels(target_raw),
           star = stars_map(p.value)) %>%
    filter(star != "", target_lab %in% lv_e70_lab[-1]) %>%
    left_join(base_y, by = "line") %>%
    mutate(x1 = pos_map_e70[lv_e70_lab[1]],
           x2 = pos_map_e70[target_lab],
           xm = (x1 + x2)/2,
           order = match(target_lab, lv_e70_lab),
           idx = ave(order, line, FUN = function(z) rank(z, ties.method = "first")),
           y  = y0 + (idx - 1) * 0.04) %>%
    dplyr::select(line, x1, x2, xm, y, star)

  p_e70 <- add_brackets_to_plot(p_e70, br_e70)

  ggsave("figures/Exp70_puncta_per_NucWithPuncta.png", p_e70, width = 7.0, height = 3.0, dpi = 600)
  ggsave("figures/Exp70_puncta_per_NucWithPuncta.svg", p_e70, width = 7.0, height = 3.0, device = svglite::svglite)

  # include the human-readable condition in the CSV
  write_csv(emm_e70 %>% relocate(conc_lab, .after = conc),
            "results/Exp70_puncta_per_NucWithPuncta_emmeans.csv")
  write_csv(pairs_e70, "results/Exp70_puncta_per_NucWithPuncta_pairs_vs_ctrl.csv")
}

# =============================== Console info ===============================
message("\nSaved outputs in: ", normalizePath("figures", winslash = "/"))
message("  - CSE_puncta_per_NucWithPuncta.(png|svg)")
message("  - KBrO3_puncta_per_NucWithPuncta.(png|svg)")
if (has_E70) message("  - Exp70_puncta_per_NucWithPuncta.(png|svg)")
message("And CSVs in results/: *_emmeans.csv, *_pairs_vs_ctrl.csv")
