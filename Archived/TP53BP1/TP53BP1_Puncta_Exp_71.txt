# TP53BP1_Puncta.R — tidy → QC → GLMs/NB → emmeans(BH) → figures (manual BH-FDR brackets) → CSVs

# ---------- Packages ----------
need <- c(
  "tidyverse","readr","readxl","janitor","glue","broom",
  "emmeans","multcomp","binom","scales","MASS","rlang","patchwork","svglite"
)
to_install <- setdiff(need, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, dep = TRUE)
invisible(lapply(need, library, character.only = TRUE))

# ============================== PALETTES ====================================
pal_sample <- c(ASC399 = "#0072B2", ASC399cse = "#E69F00")  # cool vs warm
pal_condition <- c(
  XX = "#D55E00",  XC = "#E69F00",
  CX = "#56B4E9",  CC = "#0072B2"
)
pal_pair_cool_warm <- c("#0072B2", "#E69F00")

# ---------- Locate data files (robust) ----------
root <- "F:/Repositories/ResearchGarbageBin/IF_TP53BP1"
stopifnot(dir.exists(root))

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
    all(sapply(tokens, function(t) grepl(t, bn[i], ignore.case = TRUE))),
    logical(1)
  )]
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

files <- c(
  safe_find("^9C[_-]?CSE$",   c("9C","CSE"),     root),
  safe_find("^9X[_-]?CSE$",   c("9X","CSE"),     root),
  safe_find("^9C[_-]?KBrO3$", c("9C","KBrO3"),   root),
  safe_find("^9X[_-]?KBrO3$", c("9X","KBrO3"),   root),
  safe_find("^SenExp7\\d$",   c("SenExp","70"),  root)
)
cat("Using files:\n", paste0(" - ", files, "\n"), sep = "")

# ---------- Helpers ----------
parse_context_from_filename <- function(fn) {
  tibble(
    file = fn,
    experiment = dplyr::case_when(
      grepl("CSE",   fn, ignore.case = TRUE) ~ "CSE",
      grepl("KBrO3|KBrO", fn, ignore.case = TRUE) ~ "KBrO3",
      TRUE ~ "SenExp"
    ),
    line = dplyr::case_when(
      grepl("9C", fn, ignore.case = TRUE) ~ "9C",
      grepl("9X", fn, ignore.case = TRUE) ~ "9X",
      TRUE ~ NA_character_
    )
  )
}

parse_dose_num <- function(conc_char) {
  out <- suppressWarnings(as.numeric(sub("([0-9.]+).*", "\\1", conc_char)))
  ifelse(is.na(out), NA_real_, out)
}

# level ordering without purrr::discard (avoid masking)
order_levels_by_numeric <- function(df, expt) {
  v <- df %>%
    dplyr::filter(experiment == expt) %>%
    dplyr::distinct(metadata_conc, dose_num) %>%
    dplyr::arrange(dose_num) %>%
    dplyr::pull(metadata_conc)
  v[!is.na(v)]
}
safe_levels <- function(data, expt) {
  lv <- order_levels_by_numeric(data, expt)
  if (length(lv)) lv else sort(unique(as.character(data$metadata_conc[data$experiment == expt])))
}

# Label helpers (fix: upscale fractional CSE to percent)
lab_cse <- function(x) {
  xx <- as.character(x)
  num <- suppressWarnings(as.numeric(xx))
  out <- ifelse(is.na(num), xx,
                ifelse(num <= 1, 100*num, num))
  paste0(format(out, trim = TRUE, nsmall = 1), "%")
}
lab_kbr <- function(x) {
  xx <- as.character(x)
  num <- suppressWarnings(as.numeric(xx))
  out <- ifelse(is.na(num), as.numeric(sub("([0-9.]+).*", "\\1", xx)), num)
  paste0(format(out, trim = TRUE, nsmall = 1), " mM")
}

# Column normaliser (tolerant)
unify_cp_cols <- function(df) {
  nm <- names(df)

  if (!"count_nuclei" %in% nm) {
    cand <- intersect(c("nuclei","nuclei_count","count_cells","cell_count","cells",
                        "count_nuclei_total","total_nuclei"), nm)
    if (length(cand)) df <- dplyr::rename(df, count_nuclei = !!rlang::sym(cand[1]))
  }

  if (!"count_nuclei_555pos" %in% nm) {
    cand <- intersect(c(
      "count_nuclei_555_pos","nuclei_555pos","nuclei_555_pos",
      "count_tp53bp1_pos","tp53bp1_pos_nuclei","count_nuclei_tp53bp1pos",
      "count_nuclei_tp53bp1_pos","count_555pos_nuclei",
      "count_nuclei_555positive","n_tp53bp1_pos","tp53bp1_positive_nuclei"
    ), nm)
    if (length(cand)) df <- dplyr::rename(df, count_nuclei_555pos = !!rlang::sym(cand[1]))
  }
  if (!"count_nuclei_with_dots" %in% nm) {
    cand <- intersect(c("count_nuclei_withdots","count_nuclei_with_puncta"), nm)
    if (length(cand)) df <- dplyr::rename(df, count_nuclei_with_dots = !!rlang::sym(cand[1]))
  }
  if (!"count_puncta" %in% nm) {
    cand <- intersect(c("count_children_puncta_count","puncta","count_puncta_total"), nm)
    if (length(cand)) df <- dplyr::rename(df, count_puncta = !!rlang::sym(cand[1]))
  }
  if (!"metadata_conc" %in% nm) {
    cand <- intersect(c("conc","concentration","dose","treatment",
                        "metadata_concentration","metadata_dose","metadata_cse","metadata_kbro3"), nm)
    if (length(cand)) df <- dplyr::rename(df, metadata_conc = !!rlang::sym(cand[1]))
  }
  if (!"metadata_sample" %in% nm) {
    cand <- intersect(c("sample","well","group","metadata_well","metadata_group","metadata_sampleid"), nm)
    if (length(cand)) df <- dplyr::rename(df, metadata_sample = !!rlang::sym(cand[1])) else
      df <- dplyr::mutate(df, metadata_sample = NA_character_)
  }
  if (!"imagenumber" %in% names(df)) df <- dplyr::mutate(df, imagenumber = dplyr::row_number())
  if (!"count_toobright_all" %in% names(df)) df <- dplyr::mutate(df, count_toobright_all = NA_real_)

  need <- c("count_nuclei","count_nuclei_555pos","count_nuclei_with_dots",
            "count_puncta","metadata_conc","metadata_sample","count_toobright_all")
  missing <- setdiff(need, names(df))
  if (length(missing)) stop("Missing expected columns after cleaning/normalising: ", paste(missing, collapse = ", "))
  df
}

wilson_df <- function(success, total) {
  res <- binom::binom.wilson(success, total)
  tibble(p = res$mean, lwr = res$lower, upr = res$upper)
}

std_emm_cols <- function(df) {
  nm <- names(df)
  if ("response"    %in% nm) df <- dplyr::rename(df, estimate = response)
  else if ("prob"   %in% nm) df <- dplyr::rename(df, estimate = prob)
  else if ("rate"   %in% nm) df <- dplyr::rename(df, estimate = rate)
  else if ("ratio"  %in% nm) df <- dplyr::rename(df, estimate = ratio)
  else if ("odds.ratio" %in% nm) df <- dplyr::rename(df, estimate = odds.ratio)
  else if ("emmean" %in% nm) df <- dplyr::rename(df, estimate = emmean)
  nm <- names(df); if ("asymp.LCL" %in% nm) df <- dplyr::rename(df, lcl = asymp.LCL)
  else if ("lower.CL" %in% nm) df <- dplyr::rename(df, lcl = lower.CL)
  else if ("LCL" %in% nm) df <- dplyr::rename(df, lcl = LCL)
  nm <- names(df); if ("asymp.UCL" %in% nm) df <- dplyr::rename(df, ucl = asymp.UCL)
  else if ("upper.CL" %in% nm) df <- dplyr::rename(df, ucl = upper.CL)
  else if ("UCL" %in% nm) df <- dplyr::rename(df, ucl = UCL)
  df
}

# ---------- Load & tidy ----------
read_cp <- function(path) {
  # protect readxl from seeing base::range
  rng <- NULL

  ctx <- parse_context_from_filename(basename(path))
  ext <- tolower(tools::file_ext(path))
  raw <- switch(ext,
    "xlsx" = readxl::read_excel(path, sheet = 1, range = rng),
    "xls"  = readxl::read_excel(path, sheet = 1, range = rng),
    "csv"  = readr::read_csv(path, show_col_types = FALSE),
    readr::read_tsv(path, show_col_types = FALSE)
  ) %>% janitor::clean_names()

  # Special case: SenExp sheet with 'sample' and 'treatment'
  if (identical(ctx$experiment, "SenExp")) {
    if ("sample" %in% names(raw))    raw <- dplyr::rename(raw, metadata_sample = sample)
    if ("treatment" %in% names(raw)) raw <- dplyr::rename(raw, metadata_conc   = treatment)
  }

  raw %>%
    unify_cp_cols() %>%
    dplyr::mutate(
      file = basename(path),
      experiment = ctx$experiment,
      line = dplyr::case_when(
        !is.na(ctx$line) ~ ctx$line,
        !is.na(metadata_sample) & metadata_sample %in% c("ASC399","ASC399cse") ~ metadata_sample,
        TRUE ~ NA_character_
      ),
      metadata_conc   = as.character(metadata_conc),
      metadata_sample = as.character(metadata_sample),
      dose_num = parse_dose_num(metadata_conc),
      p_555          = dplyr::if_else(count_nuclei > 0, count_nuclei_555pos / count_nuclei, NA_real_),
      p_withdots_all = dplyr::if_else(count_nuclei > 0, count_nuclei_with_dots / count_nuclei, NA_real_),
      p_withdots_555 = dplyr::if_else(count_nuclei_555pos > 0, count_nuclei_with_dots / count_nuclei_555pos, NA_real_),
      puncta_per_555 = dplyr::if_else(count_nuclei_555pos > 0, count_puncta / count_nuclei_555pos, NA_real_),
      too_bright_frac = dplyr::if_else(count_nuclei > 0 & !is.na(count_toobright_all),
                                       count_toobright_all / count_nuclei, NA_real_)
    )
}

dat <- purrr::map_dfr(files, read_cp)

# Recode lines for non-SenExp datasets
dat <- dat %>%
  dplyr::mutate(
    line = dplyr::case_when(
      experiment != "SenExp" & line == "9C" ~ "ASC399cse",
      experiment != "SenExp" & line == "9X" ~ "ASC399",
      TRUE ~ line
    )
  )

# ---------- Order factors ----------
lev_cse <- safe_levels(dat, "CSE")
lev_kbr <- safe_levels(dat, "KBrO3")

dat <- dat %>%
  dplyr::mutate(
    conc = as.character(metadata_conc),
    conc = dplyr::case_when(
      experiment == "CSE"   ~ factor(conc, levels = lev_cse),
      experiment == "KBrO3" ~ factor(conc, levels = lev_kbr),
      TRUE                  ~ factor(conc)
    ),
    line = factor(line, levels = c("ASC399","ASC399cse"))
  )

# ---------- QC ----------
qc_summary <- dat %>%
  dplyr::group_by(experiment, line, conc) %>%
  dplyr::summarise(
    images = dplyr::n(),
    total_nuclei    = sum(count_nuclei, na.rm = TRUE),
    total_555       = sum(count_nuclei_555pos, na.rm = TRUE),
    total_withdots  = sum(count_nuclei_with_dots, na.rm = TRUE),
    too_bright_pct  = 100 * mean(too_bright_frac, na.rm = TRUE),
    .groups = "drop"
  )
message("QC summary:\n"); print(qc_summary, n = Inf)

# ---------- Aggregated Wilson CIs ----------
agg_prop <- function(d, success_col, denom_col) {
  tmp <- d %>%
    dplyr::group_by(experiment, line, conc) %>%
    dplyr::summarise(
      success = sum(.data[[success_col]], na.rm = TRUE),
      denom   = sum(.data[[denom_col]],   na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(metric = success_col)
  dplyr::bind_cols(tmp, wilson_df(tmp$success, tmp$denom))
}
agg_555      <- agg_prop(dat, "count_nuclei_555pos", "count_nuclei")
agg_with_all <- agg_prop(dat, "count_nuclei_with_dots", "count_nuclei")

# ---------- Modelling ----------
check_overdisp <- function(m) { rdf <- df.residual(m); if (rdf <= 0) NA_real_ else deviance(m)/rdf }

fit_prop_glm <- function(d, success, total) {
  d2 <- d %>%
    dplyr::mutate(.success = .data[[success]], .failure = .data[[total]] - .data[[success]])
  m <- glm(cbind(.success, .failure) ~ line * conc, data = d2, family = binomial())
  disp <- check_overdisp(m)
  if (is.finite(disp) && disp > 1.5) {
    message(glue::glue("Overdispersion detected (phi ≈ {round(disp,2)}); refitting with quasibinomial."))
    m <- glm(cbind(.success, .failure) ~ line * conc, data = d2, family = quasibinomial())
  }
  m
}
fit_nb <- function(d, count_col) {
  d2 <- d %>%
    dplyr::filter(count_nuclei_555pos > 0) %>%
    dplyr::mutate(offset_log = log(count_nuclei_555pos)) %>%
    dplyr::filter(is.finite(offset_log)) %>%
    dplyr::filter(!is.na(.data[[count_col]])) %>%
    droplevels()
  f_int <- stats::as.formula(glue::glue("{count_col} ~ line * conc + offset(offset_log)"))
  f_add <- stats::as.formula(glue::glue("{count_col} ~ line + conc + offset(offset_log)"))
  m <- try(MASS::glm.nb(f_int, data = d2), silent = TRUE)
  if (inherits(m, "try-error")) m <- try(MASS::glm.nb(f_add, data = d2), silent = TRUE)
  if (inherits(m, "try-error")) {
    message("NB failed; falling back to quasi-Poisson.")
    m <- try(glm(f_int, data = d2, family = quasipoisson()), silent = TRUE)
    if (inherits(m, "try-error")) m <- glm(f_add, data = d2, family = quasipoisson())
  }
  m
}

run_models_for <- function(expt, data = dat) {
  d <- data %>% dplyr::filter(experiment == expt)
  list(
    p555      = fit_prop_glm(d, "count_nuclei_555pos",     "count_nuclei"),
    pdots_all = fit_prop_glm(d, "count_nuclei_with_dots",  "count_nuclei"),
    puncta_nb = fit_nb(d, "count_puncta")
  )
}
mods_cse <- run_models_for("CSE")
mods_kbr <- run_models_for("KBrO3")

# ---------- Pairwise tables ----------
emm_pairwise_tbl <- function(mod, specs, offset0 = FALSE) {
  if (offset0) {
    emmeans(mod, specs, type = "response", offset = 0) %>%
      contrast("trt.vs.ctrl", ref = 1) %>%
      summary(infer = TRUE, adjust = "BH") %>% as_tibble()
  } else {
    emmeans(mod, specs, type = "response") %>%
      contrast("trt.vs.ctrl", ref = 1) %>%
      summary(infer = TRUE, adjust = "BH") %>% as_tibble()
  }
}

contr_cse_pct  <- emm_pairwise_tbl(mods_cse$p555,       ~ conc | line)
contr_kbr_pct  <- emm_pairwise_tbl(mods_kbr$p555,       ~ conc | line)
contr_cse_rate <- emm_pairwise_tbl(mods_cse$puncta_nb,  ~ conc | line, offset0 = TRUE)
contr_kbr_rate <- emm_pairwise_tbl(mods_kbr$puncta_nb,  ~ conc | line, offset0 = TRUE)

# ---------- Theme ----------
theme_pub <- theme_minimal(base_size = 12) +
  theme(
    text = element_text(family = "Arial", size = 12),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )

# ---------- Manual bracket helpers (no ggsignif) ----------
stars_map <- function(p) dplyr::case_when(
  is.na(p) ~ "",
  p < 0.001 ~ "***",
  p < 0.01  ~ "**",
  p < 0.05  ~ "*",
  TRUE ~ ""
)

brackets_vs_ctrl_df <- function(contr_tbl, ci_df, levels_vec, facet_var = "line",
                                pad = 0.04) {
  if (!nrow(contr_tbl)) return(tibble())
  levels_vec <- levels_vec[!is.na(levels_vec)]
  ref <- levels_vec[1]
  contr_tbl <- contr_tbl %>%
    dplyr::mutate(
      target = sub("\\s*-\\s*.*$", "", contrast),
      star   = stars_map(p.value)
    ) %>%
    dplyr::filter(star != "", target %in% levels_vec[-1])

  if (!nrow(contr_tbl)) return(tibble())

  # y-base: max CI per facet
  base_y <- ci_df %>%
    dplyr::group_by(.data[[facet_var]]) %>%
    dplyr::summarise(y0 = max(upr, na.rm = TRUE) + pad, .groups = "drop")

  pos_map <- setNames(seq_along(levels_vec), levels_vec)

  contr_tbl %>%
    dplyr::left_join(base_y, by = setNames("line", facet_var)) %>%
    dplyr::mutate(
      order = match(target, levels_vec),
      idx   = ave(order, .data[[facet_var]], FUN = function(z) rank(z, ties.method = "first")),
      y     = y0 + (idx - 1) * pad,
      x1    = pos_map[ref],
      x2    = pos_map[target],
      xm    = (x1 + x2)/2
    ) %>%
    dplyr::select(!!facet_var, x1, x2, xm, y, star) %>%
    dplyr::rename(line = !!facet_var)
}

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

# ---------- CSE/KBrO3 Figures (with manual BH-FDR brackets) ----------
plot_p555 <- function(agg_df, raw_df, title_lab, pal, lab_fun) {
  ggplot() +
    geom_point(data = raw_df, aes(conc, p_555, colour = line),
               position = position_jitter(width = 0.08, height = 0),
               alpha = 0.7, size = 1.7) +
    geom_errorbar(data = agg_df, aes(conc, ymin = lwr, ymax = upr, colour = line),
                  width = 0.08, position = position_dodge(width = 0.25), linewidth = 0.4) +
    geom_point(data = agg_df, aes(conc, p, colour = line),
               size = 2.4, position = position_dodge(width = 0.25)) +
    geom_line(data = agg_df, aes(conc, p, colour = line, group = line),
              linewidth = 0.5, position = position_dodge(width = 0.25)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_x_discrete(labels = lab_fun, drop = TRUE) +
    scale_colour_manual(values = pal) +
    labs(x = "Condition", y = "% TP53BP1+ nuclei", title = title_lab) +
    facet_wrap(~ line, nrow = 1) + theme_pub
}

# % TP53BP1+ aggregated for both experiments
agg_p555 <- agg_555 %>% dplyr::filter(experiment %in% c("CSE","KBrO3"))

p1 <- plot_p555(
  agg_df  = agg_p555 %>% dplyr::filter(experiment == "CSE"),
  raw_df  = dat      %>% dplyr::filter(experiment == "CSE"),
  title_lab = "CSE: % TP53BP1+ nuclei",
  pal = pal_sample, lab_fun = lab_cse
)
p2 <- plot_p555(
  agg_df  = agg_p555 %>% dplyr::filter(experiment == "KBrO3"),
  raw_df  = dat      %>% dplyr::filter(experiment == "KBrO3"),
  title_lab = "KBrO₃: % TP53BP1+ nuclei",
  pal = pal_sample, lab_fun = lab_kbr
)

# Build bracket data (within line, vs control)
make_brackets_for_expt <- function(expt, contr_tbl, agg_df, raw_df, lab_fun) {
  lv <- levels(droplevels(raw_df$conc))
  ci <- agg_df %>% dplyr::select(line, conc, upr)
  brackets_vs_ctrl_df(contr_tbl, ci, lv, facet_var = "line", pad = 0.04)
}
br1 <- make_brackets_for_expt("CSE",  contr_cse_pct, agg_p555 %>% dplyr::filter(experiment=="CSE"),  dat %>% dplyr::filter(experiment=="CSE"),  lab_cse)
br2 <- make_brackets_for_expt("KBrO3",contr_kbr_pct, agg_p555 %>% dplyr::filter(experiment=="KBrO3"), dat %>% dplyr::filter(experiment=="KBrO3"), lab_kbr)

p1 <- add_brackets_to_plot(p1, br1)
p2 <- add_brackets_to_plot(p2, br2)

# Puncta per TP53BP1+ (NB means)
emm_puncta <- function(mod) as_tibble(emmeans(mod, ~ conc | line, type = "response", offset = 0)) %>% std_emm_cols()
emm_cse_puncta <- emm_puncta(mods_cse$puncta_nb)
emm_kbr_puncta <- emm_puncta(mods_kbr$puncta_nb)

plot_rate <- function(raw_df, emm_df, ttl, pal, lab_fun) {
  ggplot() +
    geom_point(data = raw_df %>% dplyr::filter(count_nuclei_555pos > 0),
               aes(conc, puncta_per_555, colour = line),
               alpha = 0.7, size = 1.7,
               position = position_jitter(width = 0.08, height = 0)) +
    geom_errorbar(data = emm_df,
                  aes(conc, ymin = lcl, ymax = ucl, colour = line),
                  width = 0.08, position = position_dodge(width = 0.25), linewidth = 0.4) +
    geom_point(data = emm_df, aes(conc, estimate, colour = line),
               size = 2.4, position = position_dodge(width = 0.25)) +
    scale_colour_manual(values = pal) +
    scale_x_discrete(labels = lab_fun, drop = TRUE) +
    labs(x = "Condition", y = "Puncta per TP53BP1+ nucleus", title = ttl) +
    facet_wrap(~ line, nrow = 1) + theme_pub
}

p3 <- plot_rate(dat %>% dplyr::filter(experiment=="CSE"),   emm_cse_puncta, "CSE: puncta per TP53BP1+ nucleus", pal_sample, lab_cse)
p4 <- plot_rate(dat %>% dplyr::filter(experiment=="KBrO3"), emm_kbr_puncta, "KBrO₃: puncta per TP53BP1+ nucleus", pal_sample, lab_kbr)

# Brackets for rate plots
ci_cse_rate <- emm_cse_puncta %>% dplyr::group_by(line, conc) %>% dplyr::summarise(upr = max(ucl), .groups="drop")
ci_kbr_rate <- emm_kbr_puncta %>% dplyr::group_by(line, conc) %>% dplyr::summarise(upr = max(ucl), .groups="drop")
br3 <- brackets_vs_ctrl_df(contr_cse_rate, ci_cse_rate, levels(droplevels((dat %>% dplyr::filter(experiment=="CSE"))$conc)))
br4 <- brackets_vs_ctrl_df(contr_kbr_rate, ci_kbr_rate, levels(droplevels((dat %>% dplyr::filter(experiment=="KBrO3"))$conc)))
p3 <- add_brackets_to_plot(p3, br3)
p4 <- add_brackets_to_plot(p4, br4)

# ---------- Across-sample (ASC399 vs ASC399cse) — plot only significant ----------
# For each experiment and metric, test sample (line) difference at each conc
across_sig_plot <- function(expt, mod_prop, mod_rate, title_prefix, lab_fun) {
  # prop metric
  comp_prop <- emmeans(mod_prop, ~ line | conc, type = "response") %>%
    contrast("pairwise") %>% summary(infer = TRUE, adjust = "BH") %>% as_tibble()
  comp_rate <- emmeans(mod_rate, ~ line | conc, type = "response", offset = 0) %>%
    contrast("pairwise") %>% summary(infer = TRUE, adjust = "BH") %>% as_tibble()

  keep_prop <- comp_prop %>% dplyr::filter(p.value < 0.05) %>% dplyr::pull(conc) %>% unique()
  keep_rate <- comp_rate %>% dplyr::filter(p.value < 0.05) %>% dplyr::pull(conc) %>% unique()

  make_panel <- function(metric, keep_levels, ylab) {
    if (!length(keep_levels)) return(NULL)
    d <- dat %>% dplyr::filter(experiment == expt, conc %in% keep_levels)
    agg <- if (metric == "pct") {
      agg_555 %>% dplyr::filter(experiment == expt, conc %in% keep_levels) %>% dplyr::rename(est = p, lwr = lwr, upr = upr)
    } else {
      if (expt == "CSE") emm_cse_puncta else emm_kbr_puncta
    }
    pos <- position_dodge(width = 0.35)

    p <- ggplot() +
      geom_point(data = d %>% dplyr::distinct(line, conc),
                 aes(x = line, y = 0, colour = line), alpha = 0) + # for legend
      {
        if (metric == "pct")
          geom_point(data = agg, aes(conc, est, colour = line), position = pos, size = 2.4)
        else
          geom_point(data = agg, aes(conc, estimate, colour = line), position = pos, size = 2.4)
      } +
      {
        if (metric == "pct")
          geom_errorbar(data = agg, aes(conc, ymin = lwr, ymax = upr, colour = line),
                        width = 0.08, position = pos, linewidth = 0.4)
        else
          geom_errorbar(data = agg, aes(conc, ymin = lcl, ymax = ucl, colour = line),
                        width = 0.08, position = pos, linewidth = 0.4)
      } +
      scale_colour_manual(values = pal_sample) +
      scale_x_discrete(labels = lab_fun, drop = TRUE) +
      labs(x = "Condition", y = ylab,
           title = paste0(title_prefix, if (metric=="pct") ": % TP53BP1+ (sig only)" else ": puncta/TP53BP1+ (sig only)")) +
      facet_wrap(~ line, nrow = 1) + theme_pub
    p
  }

  list(
    pct  = make_panel("pct",  keep_prop, "% TP53BP1+ nuclei"),
    rate = make_panel("rate", keep_rate, "Puncta per TP53BP1+ nucleus")
  )
}

ac_cse  <- across_sig_plot("CSE",  mods_cse$p555, mods_cse$puncta_nb,  "CSE",  lab_cse)
ac_kbr  <- across_sig_plot("KBrO3",mods_kbr$p555, mods_kbr$puncta_nb,  "KBrO₃", lab_kbr)

# ---------- Senescence (CC/CX/XC/XX by sample; brackets vs CC) ----------
sen <- dat %>% dplyr::filter(experiment == "SenExp")
if (nrow(sen)) {
  # ensure sample names come from file
  sen2 <- sen %>%
    dplyr::mutate(
      code = factor(as.character(conc), levels = c("CC","CX","XC","XX")),
      sample = factor(dplyr::coalesce(metadata_sample, line), levels = c("ASC399","ASC399cse"))
    ) %>% dplyr::filter(!is.na(code), !is.na(sample))

  # Model
  mod_sen <- glm(
    cbind(count_nuclei_555pos, count_nuclei - count_nuclei_555pos) ~ sample * code,
    data = sen2, family = binomial()
  )
  if (is.finite(check_overdisp(mod_sen)) && check_overdisp(mod_sen) > 1.5) {
    mod_sen <- glm(
      cbind(count_nuclei_555pos, count_nuclei - count_nuclei_555pos) ~ sample * code,
      data = sen2, family = quasibinomial()
    )
  }

  emm_sen <- emmeans(mod_sen, ~ code | sample, type = "response")
  sen_pairs <- contrast(emm_sen, "trt.vs.ctrl", ref = 1) %>%
    summary(infer = TRUE, adjust = "BH") %>% as_tibble()

  agg_sen <- sen2 %>%
    dplyr::group_by(sample, code) %>%
    dplyr::summarise(success = sum(count_nuclei_555pos), denom = sum(count_nuclei), .groups = "drop") %>%
    dplyr::bind_cols(wilson_df(.$success, .$denom))

  lv_sen <- c("CC","CX","XC","XX")
  ci_sen <- agg_sen %>% dplyr::select(line = sample, conc = code, upr)
  br_sen <- sen_pairs %>% dplyr::mutate(line = sample)
  br_sen <- brackets_vs_ctrl_df(br_sen, ci_sen, lv_sen, facet_var = "line", pad = 0.05)

  p5 <- ggplot() +
    geom_point(data = sen2, aes(code, p_555, colour = code),
               position = position_jitter(width = 0.08, height = 0), alpha = 0.7, size = 1.7) +
    geom_errorbar(data = agg_sen, aes(code, ymin = lwr, ymax = upr, colour = code),
                  width = 0.10, linewidth = 0.4) +
    geom_point(data = agg_sen, aes(code, p, colour = code), size = 2.4) +
    scale_colour_manual(values = pal_condition, breaks = lv_sen) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.05),
                       expand = expansion(mult = c(0.02, 0.08))) +
    labs(x = "Condition code", y = "% TP53BP1+ nuclei",
         title = "Senescence: % TP53BP1+ nuclei (BH-FDR vs CC)") +
    facet_wrap(~ sample, nrow = 1) +
    theme_pub + theme(legend.position = "bottom")

  p5 <- add_brackets_to_plot(p5, br_sen)

  dir.create("figures", showWarnings = FALSE)
  ggsave("figures/Sen_TP53BP1pos.png", p5, width = 7.8, height = 3.2, dpi = 600)
  ggsave("figures/Sen_TP53BP1pos.svg", p5, width = 7.8, height = 3.2, device = svglite::svglite)
}

# ---------- Export overview figures ----------
dir.create("figures", showWarnings = FALSE)
ggsave("figures/CSE_TP53BP1pos.png",       p1, width = 7.0, height = 3.0, dpi = 600)
ggsave("figures/KBrO3_TP53BP1pos.png",     p2, width = 7.0, height = 3.0, dpi = 600)
ggsave("figures/CSE_puncta_TP53BP1pos.png",  p3, width = 7.0, height = 3.0, dpi = 600)
ggsave("figures/KBrO3_puncta_TP53BP1pos.png", p4, width = 7.0, height = 3.0, dpi = 600)
ggsave("figures/CSE_TP53BP1pos.svg",       p1, width = 7.0, height = 3.0, device = svglite::svglite)
ggsave("figures/KBrO3_TP53BP1pos.svg",     p2, width = 7.0, height = 3.0, device = svglite::svglite)
ggsave("figures/CSE_puncta_TP53BP1pos.svg",  p3, width = 7.0, height = 3.0, device = svglite::svglite)
ggsave("figures/KBrO3_puncta_TP53BP1pos.svg", p4, width = 7.0, height = 3.0, device = svglite::svglite)

# ---------- Export across-sample (sig only) if present ----------
if (!is.null(ac_cse$pct))  ggsave("figures/CSE_sig_across_samples_pct.png",  ac_cse$pct,  width = 7.0, height = 3.0, dpi = 600)
if (!is.null(ac_cse$rate)) ggsave("figures/CSE_sig_across_samples_rate.png", ac_cse$rate, width = 7.0, height = 3.0, dpi = 600)
if (!is.null(ac_kbr$pct))  ggsave("figures/KBrO3_sig_across_samples_pct.png",  ac_kbr$pct,  width = 7.0, height = 3.0, dpi = 600)
if (!is.null(ac_kbr$rate)) ggsave("figures/KBrO3_sig_across_samples_rate.png", ac_kbr$rate, width = 7.0, height = 3.0, dpi = 600)

# ---------- RESULTS CSV ----------
dir.create("results", showWarnings = FALSE)

add_dose_label <- function(df) {
  df %>% dplyr::mutate(
    dose_label = dplyr::case_when(
      experiment == "CSE"   ~ lab_cse(conc),
      experiment == "KBrO3" ~ lab_kbr(conc),
      TRUE ~ as.character(conc)
    )
  )
}

collect_prop_results <- function(mod, experiment, metric) {
  means <- emmeans(mod, ~ conc | line, type = "response") %>% as_tibble() %>% std_emm_cols() %>%
    dplyr::mutate(result_type = "mean", experiment = experiment, metric = metric)
  pairs <- emmeans(mod, ~ conc | line, type = "response") %>%
    contrast("trt.vs.ctrl", ref = 1) %>%
    summary(infer = TRUE, adjust = "BH") %>% as_tibble() %>% std_emm_cols() %>%
    dplyr::mutate(result_type = "contrast_vs_ctrl", experiment = experiment, metric = metric)
  dplyr::bind_rows(means, pairs)
}
collect_rate_results <- function(mod, experiment, metric) {
  means <- emmeans(mod, ~ conc | line, type = "response", offset = 0) %>% as_tibble() %>% std_emm_cols() %>%
    dplyr::mutate(result_type = "mean", experiment = experiment, metric = metric)
  pairs <- emmeans(mod, ~ conc | line, type = "response", offset = 0) %>%
    contrast("trt.vs.ctrl", ref = 1) %>%
    summary(infer = TRUE, adjust = "BH") %>% as_tibble() %>% std_emm_cols() %>%
    dplyr::mutate(result_type = "contrast_vs_ctrl", experiment = experiment, metric = metric)
  dplyr::bind_rows(means, pairs)
}

results_all <- dplyr::bind_rows(
  collect_prop_results(mods_cse$p555,      "CSE",   "pct_555_pos"),
  collect_prop_results(mods_cse$pdots_all, "CSE",   "pct_with_dots_all"),
  collect_rate_results(mods_cse$puncta_nb, "CSE",   "puncta_per_555"),
  collect_prop_results(mods_kbr$p555,      "KBrO3", "pct_555_pos"),
  collect_prop_results(mods_kbr$pdots_all, "KBrO3", "pct_with_dots_all"),
  collect_rate_results(mods_kbr$puncta_nb, "KBrO3", "puncta_per_555")
) %>%
  dplyr::mutate(conc = as.character(conc)) %>%
  add_dose_label() %>%
  dplyr::relocate(experiment, line, metric, result_type, conc, dose_label)

readr::write_csv(results_all, "results/summary_results.csv")

# ---------- JOURNAL-READY WIDE TABLES ----------
metric_pretty_map <- c(
  pct_555_pos       = "% TP53BP1+ nuclei",
  pct_with_dots_all = "% with dots (all)",
  puncta_per_555    = "Puncta per TP53BP1+ nucleus"
)
fmt_ci_pct <- function(est, lcl, ucl) paste0(scales::percent(est, accuracy = 0.1), " (",
                                             scales::percent(lcl, accuracy = 0.1), "–",
                                             scales::percent(ucl, accuracy = 0.1), ")")
fmt_ci_num <- function(est, lcl, ucl) sprintf("%.2f (%.2f–%.2f)", est, lcl, ucl)

means_all <- results_all %>%
  dplyr::filter(result_type == "mean") %>%
  dplyr::mutate(
    metric_pretty = dplyr::recode(metric, !!!metric_pretty_map),
    display = dplyr::case_when(
      metric == "puncta_per_555" ~ fmt_ci_num(estimate, lcl, ucl),
      TRUE                       ~ fmt_ci_pct(estimate, lcl, ucl)
    )
  )

dose_cols_cse <- lab_cse(lev_cse)
dose_cols_kbr <- lab_kbr(lev_kbr)

write_means_wide <- function(df, expt, dose_cols, out_path) {
  if (!any(df$experiment == expt)) return(invisible(NULL))
  wide <- df %>%
    dplyr::filter(experiment == expt) %>%
    dplyr::mutate(dose_label = factor(dose_label, levels = dose_cols)) %>%
    dplyr::select(metric_pretty, line, dose_label, display) %>%
    dplyr::arrange(metric_pretty, line, dose_label) %>%
    tidyr::pivot_wider(names_from = dose_label, values_from = display) %>%
    dplyr::arrange(metric_pretty, line)
  readr::write_csv(wide, out_path)
}

write_means_wide(means_all, "CSE",   dose_cols_cse, file.path("results","CSE_means_wide.csv"))
write_means_wide(means_all, "KBrO3", dose_cols_kbr, file.path("results","KBrO3_means_wide.csv"))

message("\nDone. Overview plots → ./figures/ (with BH-FDR brackets), across-sample significant plots → ./figures/*sig_across_*.png, tables → ./results/")
