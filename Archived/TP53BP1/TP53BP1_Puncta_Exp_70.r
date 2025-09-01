# SenExp70_only.R — CC/CX/XC/XX within SenExp70 (TP53BP1+)
# Outputs:
#  - figures/SenExp70_TP53BP1pos.png/.svg (per-sample, BH-FDR vs CC)
#  - results/SenExp70_pairs_vs_CC.csv, SenExp70_pairs_all.csv
#  - results/SenExp70_means_wide.csv

# ---------- Packages ----------
need <- c("tidyverse","readr","readxl","janitor","glue","emmeans","binom","scales","MASS","rlang","svglite")
to_install <- setdiff(need, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, dep = TRUE)
invisible(lapply(need, library, character.only = TRUE))

# ---------- Paths ----------
root <- "F:/Repositories/ResearchGarbageBin/IF_TP53BP1"
stopifnot(dir.exists(root))

# ---------- Helpers ----------
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

# ---------- Load & tidy ----------
parse_context_from_filename <- function(fn) tibble(file = basename(fn), experiment = "SenExp")
parse_dose_num <- function(conc_char) {
  out <- suppressWarnings(as.numeric(sub("([0-9.]+).*", "\\1", conc_char)))
  ifelse(is.na(out), NA_real_, out)
}
unify_cp_cols <- function(df) {
  nm <- names(df)
  if (!"count_nuclei" %in% nm) {
    cand <- intersect(c("nuclei","nuclei_count","count_cells","cell_count","cells",
                        "count_nuclei_total","total_nuclei"), nm)
    if (length(cand)) df <- dplyr::rename(df, count_nuclei = !!rlang::sym(cand[1]))
  }
  if (!"count_nuclei_555pos" %in% nm) {
    cand <- intersect(c("count_nuclei_555_pos","nuclei_555pos","nuclei_555_pos",
                        "count_tp53bp1_pos","tp53bp1_pos_nuclei","count_nuclei_tp53bp1pos",
                        "count_nuclei_tp53bp1_pos","count_555pos_nuclei",
                        "count_nuclei_555positive","n_tp53bp1_pos","tp53bp1_positive_nuclei"), nm)
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
                        "metadata_concentration","metadata_dose"), nm)
    if (length(cand)) df <- dplyr::rename(df, metadata_conc = !!rlang::sym(cand[1]))
  }
  if (!"metadata_sample" %in% nm) {
    cand <- intersect(c("sample","well","group","metadata_well","metadata_group","metadata_sampleid"), nm)
    if (length(cand)) df <- dplyr::rename(df, metadata_sample = !!rlang::sym(cand[1])) else
      df <- dplyr::mutate(df, metadata_sample = NA_character_)
  }
  if (!"count_toobright_all" %in% nm) df <- dplyr::mutate(df, count_toobright_all = NA_real_)
  need <- c("count_nuclei","count_nuclei_555pos","count_nuclei_with_dots","count_puncta",
            "metadata_conc","metadata_sample","count_toobright_all")
  missing <- setdiff(need, names(df))
  if (length(missing)) stop("Missing expected columns after cleaning/normalising: ", paste(missing, collapse = ", "))
  df
}
wilson_df <- function(success, total) {
  res <- binom::binom.wilson(success, total)
  tibble(p = res$mean, lwr = res$lower, upr = res$upper)
}
check_overdisp <- function(m) { rdf <- df.residual(m); if (rdf <= 0) NA_real_ else deviance(m)/rdf }

# robust reader (avoid base::range collision)
read_cp <- function(path) {
  rng <- NULL
  ctx <- parse_context_from_filename(basename(path))
  ext <- tolower(tools::file_ext(path))
  raw <- switch(ext,
                "xlsx" = readxl::read_excel(path, sheet = 1, range = rng),
                "xls"  = readxl::read_excel(path, sheet = 1, range = rng),
                "csv"  = readr::read_csv(path, show_col_types = FALSE),
                readr::read_tsv(path, show_col_types = FALSE)) %>%
    janitor::clean_names()

  # SenExp sheet often uses 'sample' and 'treatment'
  if ("sample" %in% names(raw))    raw <- dplyr::rename(raw, metadata_sample = sample)
  if ("treatment" %in% names(raw)) raw <- dplyr::rename(raw, metadata_conc   = treatment)

  raw %>%
    unify_cp_cols() %>%
    dplyr::mutate(
      file = basename(path),
      experiment = ctx$experiment,
      line = dplyr::coalesce(metadata_sample, NA_character_),
      metadata_conc   = as.character(metadata_conc),
      metadata_sample = as.character(metadata_sample),
      dose_num = parse_dose_num(metadata_conc),
      p_555          = dplyr::if_else(count_nuclei > 0, count_nuclei_555pos / count_nuclei, NA_real_),
      puncta_per_555 = dplyr::if_else(count_nuclei_555pos > 0, count_puncta / count_nuclei_555pos, NA_real_),
      too_bright_frac = dplyr::if_else(count_nuclei > 0 & !is.na(count_toobright_all),
                                       count_toobright_all / count_nuclei, NA_real_)
    )
}

# ---------- Get SenExp70 file only ----------
sen_file <- safe_find("^SenExp70$", c("SenExp","70"), root)
cat("Using file:\n - ", sen_file, "\n", sep = "")

dat <- read_cp(sen_file) %>%
  dplyr::mutate(
    code = factor(as.character(metadata_conc), levels = c("CC","CX","XC","XX")),
    sample = factor(dplyr::coalesce(metadata_sample, line))
  ) %>%
  dplyr::filter(!is.na(code), !is.na(sample))

stopifnot(nrow(dat) > 0)

# ---------- QC ----------
qc <- dat %>%
  dplyr::group_by(sample, code) %>%
  dplyr::summarise(
    n_images   = dplyr::n(),
    nuclei     = sum(count_nuclei, na.rm = TRUE),
    tp53bp1pos = sum(count_nuclei_555pos, na.rm = TRUE),
    .groups = "drop"
  )
message("QC summary:\n"); print(qc, n = Inf)

# ---------- Aggregated means (Wilson) ----------
agg <- dat %>%
  dplyr::group_by(sample, code) %>%
  dplyr::summarise(success = sum(count_nuclei_555pos), denom = sum(count_nuclei), .groups = "drop") %>%
  dplyr::bind_cols(wilson_df(.$success, .$denom))

# ---------- Modelling (per-sample contrasts vs CC; BH) ----------
mod <- glm(
  cbind(count_nuclei_555pos, count_nuclei - count_nuclei_555pos) ~ sample * code,
  data = dat, family = binomial()
)
phi <- check_overdisp(mod)
if (is.finite(phi) && phi > 1.5) {
  message(glue::glue("Overdispersion detected (phi ≈ {round(phi,2)}); refitting with quasibinomial."))
  mod <- glm(
    cbind(count_nuclei_555pos, count_nuclei - count_nuclei_555pos) ~ sample * code,
    data = dat, family = quasibinomial()
  )
}

emm_code_by_sample <- emmeans(mod, ~ code | sample, type = "response")
pairs_vs_CC <- contrast(emm_code_by_sample, "trt.vs.ctrl", ref = 1) %>%
  summary(infer = TRUE, adjust = "BH") %>% as_tibble()

pairs_all <- contrast(emm_code_by_sample, "pairwise") %>%
  summary(infer = TRUE, adjust = "BH") %>% as_tibble()

# ---------- Plot (per-sample with BH-FDR brackets vs CC) ----------
pal_condition <- c(CC="#0072B2", CX="#56B4E9", XC="#E69F00", XX="#D55E00")
theme_pub <- theme_minimal(base_size = 12) +
  theme(
    text = element_text(family = "Arial", size = 12),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )
stars_map <- function(p) dplyr::case_when(
  is.na(p) ~ "",
  p < 0.001 ~ "***",
  p < 0.01  ~ "**",
  p < 0.05  ~ "*",
  TRUE ~ ""
)
brackets_vs_ctrl_df <- function(contr_tbl, ci_df, levels_vec, facet_var = "sample", pad = 0.05) {
  if (!nrow(contr_tbl)) return(tibble())
  contr_tbl <- contr_tbl %>%
    dplyr::mutate(target = sub("\\s*-\\s*.*$", "", contrast),
                  star = stars_map(p.value)) %>%
    dplyr::filter(star != "", target %in% levels_vec[-1])
  if (!nrow(contr_tbl)) return(tibble())

  base_y <- ci_df %>%
    dplyr::group_by(.data[[facet_var]]) %>%
    dplyr::summarise(y0 = max(upr, na.rm = TRUE) + pad, .groups = "drop")

  pos_map <- setNames(seq_along(levels_vec), levels_vec)
  contr_tbl %>%
    dplyr::left_join(base_y, by = setNames("sample", facet_var)) %>%
    dplyr::mutate(
      order = match(target, levels_vec),
      idx   = ave(order, .data[[facet_var]], FUN = function(z) rank(z, ties.method = "first")),
      y     = y0 + (idx - 1) * pad,
      x1    = pos_map[levels_vec[1]],
      x2    = pos_map[target],
      xm    = (x1 + x2)/2
    ) %>%
    dplyr::select(!!facet_var, x1, x2, xm, y, star) %>%
    dplyr::rename(sample = !!facet_var)
}

lv_codes <- c("CC","CX","XC","XX")
ci_df <- agg %>% dplyr::select(sample, code, upr)
br_df <- brackets_vs_ctrl_df(pairs_vs_CC %>% dplyr::mutate(sample = as.character(sample)),
                             ci_df %>% dplyr::rename(code = code), lv_codes, facet_var = "sample")

p <- ggplot() +
  geom_point(data = dat, aes(code, p_555, colour = code),
             position = position_jitter(width = 0.08, height = 0), alpha = 0.7, size = 1.7) +
  geom_errorbar(data = agg, aes(code, ymin = lwr, ymax = upr, colour = code),
                width = 0.10, linewidth = 0.4) +
  geom_point(data = agg, aes(code, p, colour = code), size = 2.4) +
  scale_colour_manual(values = pal_condition, breaks = lv_codes) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1.05),
                     expand = expansion(mult = c(0.02, 0.08))) +
  labs(x = "Condition code", y = "% TP53BP1+ nuclei",
       title = "SenExp70: % TP53BP1+ nuclei (BH-FDR vs CC)") +
  facet_wrap(~ sample, nrow = 1) +
  theme_pub

# add brackets
if (nrow(br_df)) {
  tick <- 0.012
  p <- p +
    geom_segment(data = br_df, aes(x = x1, xend = x2, y = y, yend = y),
                 inherit.aes = FALSE, linewidth = 0.4) +
    geom_segment(data = br_df, aes(x = x1, xend = x1, y = y, yend = y - tick),
                 inherit.aes = FALSE, linewidth = 0.4) +
    geom_segment(data = br_df, aes(x = x2, xend = x2, y = y, yend = y - tick),
                 inherit.aes = FALSE, linewidth = 0.4) +
    geom_text(data = br_df, aes(x = xm, y = y + 0.8*tick, label = star),
              inherit.aes = FALSE, vjust = 0, size = 3.5)
}

# ---------- Save ----------
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
ggsave("figures/SenExp70_TP53BP1pos.png", p, width = 7.8, height = 3.2, dpi = 600)
ggsave("figures/SenExp70_TP53BP1pos.svg", p, width = 7.8, height = 3.2, device = svglite::svglite)

readr::write_csv(pairs_vs_CC, "results/SenExp70_pairs_vs_CC.csv")
readr::write_csv(pairs_all,   "results/SenExp70_pairs_all.csv")

# Pretty means (journal-ready wide)
fmt_ci_pct <- function(est, lcl, ucl) paste0(scales::percent(est, accuracy = 0.1), " (",
                                             scales::percent(lcl, accuracy = 0.1), "–",
                                             scales::percent(ucl, accuracy = 0.1), ")")
means_disp <- agg %>%
  dplyr::mutate(display = fmt_ci_pct(p, lwr, upr)) %>%
  dplyr::select(sample, code, display) %>%
  tidyr::pivot_wider(names_from = code, values_from = display) %>%
  dplyr::arrange(sample)

readr::write_csv(means_disp, "results/SenExp70_means_wide.csv")

message("\nDone. SenExp70-only contrasts and figure saved in ./results and ./figures.")