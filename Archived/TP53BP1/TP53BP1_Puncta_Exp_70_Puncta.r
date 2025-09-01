# ======================= Exp70 — puncta per TP53BP1+ nucleus =================
# Outputs:
#   figures/Exp70_puncta_per_TP53BP1pos.png / .svg
#   results/Exp70_puncta_per_TP53BP1pos_emmeans.csv
#   results/Exp70_puncta_per_TP53BP1pos_pairs_vs_ctrl.csv

suppressPackageStartupMessages({
  need <- c("tidyverse","readr","readxl","janitor","glue","emmeans","MASS","svglite","rlang","scales")
  to_install <- setdiff(need, rownames(installed.packages()))
  if (length(to_install)) install.packages(to_install, dep = TRUE)
  invisible(lapply(need, library, character.only = TRUE))
})

# ---------- Helpers ----------
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
parse_dose_num <- function(conc_char) {
  out <- suppressWarnings(as.numeric(sub("([0-9.]+).*", "\\1", as.character(conc_char))))
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
                        "count_nuclei_555positive","n_tp53bp1_pos","tp53bp1_positive_nuclei"), nm)
    if (length(cand)) df <- dplyr::rename(df, count_nuclei_555pos = !!rlang::sym(cand[1]))
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
  need <- c("count_nuclei","count_puncta","count_nuclei_555pos","metadata_conc","metadata_sample")
  missing <- setdiff(need, names(df))
  if (length(missing)) stop("Missing expected columns: ", paste(missing, collapse = ", "))
  df
}
theme_pub <- theme_minimal(base_size = 12) +
  theme(
    text = element_text(family = "Arial", size = 12),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )
stars_map <- function(p) dplyr::case_when(
  is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ ""
)

# ---------- Paths ----------
root <- "F:/Repositories/ResearchGarbageBin/IF_TP53BP1"
stopifnot(dir.exists(root))
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

.list_all <- function(root) {
  exts <- c("xlsx","xls","csv","tsv","txt")
  pats <- paste0("\\.(", paste(exts, collapse="|"), ")$")
  all <- list.files(root, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  all[grepl(pats, all, ignore.case = TRUE)]
}
find_exp70 <- function(root) {
  all <- .list_all(root)
  bn  <- basename(all)
  keep <- all[grepl("exp[_-]?70", bn, ignore.case = TRUE)]
  if (!length(keep)) stop("No file found for Exp70 in: ", root,
                          "\nExpected filename to contain 'Exp70' or 'Exp_70'.")
  normalizePath(keep[1], mustWork = TRUE)
}

# ---------- Load Exp70 ----------
exp70_path <- find_exp70(root)
cat("Using Exp70 file:\n - ", exp70_path, "\n", sep = "")

rng <- NULL
ext <- tolower(tools::file_ext(exp70_path))
raw <- switch(ext,
  "xlsx" = readxl::read_excel(exp70_path, sheet = 1, range = rng),
  "xls"  = readxl::read_excel(exp70_path, sheet = 1, range = rng),
  "csv"  = readr::read_csv(exp70_path, show_col_types = FALSE),
  readr::read_tsv(exp70_path, show_col_types = FALSE)
) %>% janitor::clean_names()

dat70 <- raw %>%
  unify_cp_cols() %>%
  dplyr::mutate(
    conc   = as.character(metadata_conc),
    line   = as.character(metadata_sample),
    dose_n = parse_dose_num(conc),
    puncta_per_555 = dplyr::if_else(count_nuclei_555pos > 0, count_puncta / count_nuclei_555pos, NA_real_)
  ) %>%
  dplyr::filter(!is.na(conc), !is.na(line))

# factor order: numeric if possible, else alphabetical
lev_conc <- dat70 %>% dplyr::distinct(conc, dose_n) %>% dplyr::arrange(dose_n) %>% dplyr::pull(conc)
if (all(is.na(dat70$dose_n))) lev_conc <- sort(unique(dat70$conc))
dat70 <- dat70 %>% dplyr::mutate(conc = factor(conc, levels = lev_conc))

# palette for lines (dynamic; falls back if names differ)
u_lines <- sort(unique(dat70$line))
pal_sample <- if (setequal(u_lines, c("ASC399","ASC399cse"))) {
  c(ASC399 = "#0072B2", ASC399cse = "#E69F00")
} else {
  setNames(scales::hue_pal()(length(u_lines)), u_lines)
}

# ---------- Model: NB with offset(log(TP53BP1+ nuclei)) ----------
d_nb <- dat70 %>%
  dplyr::filter(count_nuclei_555pos > 0) %>%
  dplyr::mutate(offset_log = log(count_nuclei_555pos)) %>%
  dplyr::filter(is.finite(offset_log), !is.na(count_puncta))
f_int <- as.formula("count_puncta ~ line * conc + offset(offset_log)")
f_add <- as.formula("count_puncta ~ line + conc + offset(offset_log)")
m_nb <- try(MASS::glm.nb(f_int, data = d_nb), silent = TRUE)
if (inherits(m_nb, "try-error")) m_nb <- try(MASS::glm.nb(f_add, data = d_nb), silent = TRUE)
if (inherits(m_nb, "try-error")) {
  message("NB failed; using quasi-Poisson.")
  m_nb <- try(glm(f_int,  data = d_nb, family = quasipoisson()), silent = TRUE)
  if (inherits(m_nb, "try-error")) m_nb <- glm(f_add, data = d_nb, family = quasipoisson())
}

# EMMEANS: predicted puncta per TP53BP1+ nucleus (robust column names)
emm70 <- emmeans(m_nb, ~ conc | line, type = "response", offset = 0) %>%
  as_tibble() %>%
  std_emm_cols()

# Pairwise vs control (first conc level per line), BH-FDR
pairs70 <- emmeans(m_nb, ~ conc | line, type = "response", offset = 0) %>%
  contrast("trt.vs.ctrl", ref = 1) %>%
  summary(infer = TRUE, adjust = "BH") %>%
  as_tibble()

readr::write_csv(emm70,  "results/Exp70_puncta_per_TP53BP1pos_emmeans.csv")
readr::write_csv(pairs70,"results/Exp70_puncta_per_TP53BP1pos_pairs_vs_ctrl.csv")

# ---------- Brackets ----------
ci_top <- emm70 %>%
  dplyr::group_by(line, conc) %>%
  dplyr::summarise(upr = max(ucl), .groups = "drop")
levels_vec <- levels(dat70$conc)
pos_map <- setNames(seq_along(levels_vec), levels_vec)
br_df <- pairs70 %>%
  dplyr::mutate(target = sub("\\s*-\\s*.*$", "", contrast),
                star = stars_map(p.value)) %>%
  dplyr::filter(star != "", target %in% levels_vec[-1]) %>%
  dplyr::left_join(ci_top %>% dplyr::group_by(line) %>% dplyr::summarise(y0 = max(upr, na.rm = TRUE) + 0.04),
                   by = "line") %>%
  dplyr::mutate(
    x1 = pos_map[levels_vec[1]],
    x2 = pos_map[target],
    xm = (x1 + x2)/2,
    order = match(target, levels_vec),
    idx   = ave(order, line, FUN = function(z) rank(z, ties.method = "first")),
    y     = y0 + (idx - 1) * 0.04
  ) %>%
  dplyr::select(line, x1, x2, xm, y, star)

# ---------- Plot ----------
p_exp70 <- ggplot() +
  geom_point(data = dat70 %>% dplyr::filter(!is.na(puncta_per_555)),
             aes(conc, puncta_per_555, colour = line),
             alpha = 0.7, size = 1.7,
             position = position_jitter(width = 0.08, height = 0)) +
  geom_errorbar(data = emm70,
                aes(conc, ymin = lcl, ymax = ucl, colour = line),
                width = 0.08, position = position_dodge(width = 0.25), linewidth = 0.4) +
  geom_point(data = emm70,
             aes(conc, estimate, colour = line),
             size = 2.4, position = position_dodge(width = 0.25)) +
  scale_colour_manual(values = pal_sample) +
  labs(x = "Condition", y = "Puncta per TP53BP1+ nucleus",
       title = "Exp70: puncta per TP53BP1+ nucleus (NB means \u00B1 95% CI; BH-FDR vs control)") +
  facet_wrap(~ line, nrow = 1, scales = "free_y") +
  theme_pub

if (nrow(br_df)) {
  tick <- 0.012
  p_exp70 <- p_exp70 +
    geom_segment(data = br_df, aes(x = x1, xend = x2, y = y, yend = y),
                 inherit.aes = FALSE, linewidth = 0.4) +
    geom_segment(data = br_df, aes(x = x1, xend = x1, y = y, yend = y - tick),
                 inherit.aes = FALSE, linewidth = 0.4) +
    geom_segment(data = br_df, aes(x = x2, xend = x2, y = y, yend = y - tick),
                 inherit.aes = FALSE, linewidth = 0.4) +
    geom_text(   data = br_df, aes(x = xm, y = y + 0.8*tick, label = star),
                 inherit.aes = FALSE, vjust = 0, size = 3.5)
}

ggsave("figures/Exp70_puncta_per_TP53BP1pos.png", p_exp70, width = 7.0, height = 3.0, dpi = 600)
ggsave("figures/Exp70_puncta_per_TP53BP1pos.svg", p_exp70, width = 7.0, height = 3.0, device = svglite::svglite)

message("Saved: figures/Exp70_puncta_per_TP53BP1pos.(png|svg)")
