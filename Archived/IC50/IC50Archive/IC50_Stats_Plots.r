###############################################################################
# Resazurin Exp35 — IC50 publication plots + robust stats          2025-08-26
# Input: Resazurin_Exp35_IC50_results.csv  (columns: Patient, Group, CompoundName, Unit/UnitStd, IC50)
# Outputs (in Exp35_IC50_Plots/):
#   - Stats_Modifiers_vsBase.csv
#   - Plot_B_Modifier_Effects.pdf
#   - Plot_C_PairedSpaghetti_Modifiers.pdf
#   - Plot_D_Forest_Modifiers_vsBase.pdf
#   - (optional) Stats_CSE_vs_Base_byCompound.csv and Plot_A_CSE_Effect_byCompound.pdf if CSE data exist
###############################################################################

## ---- Paths ------------------------------------------------------------------
in_file <- "IC50/Resazurin_Exp_35/Resazurin_Exp35_IC50_results.csv"  # change if needed
out_dir <- "Exp35_IC50_Plots"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- Libraries ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(purrr); library(scales)
})

set.seed(1L)

## ---- Settings ----------------------------------------------------------------
treated_regex <- "(?i)(?:^|[ _\\-\\(])cse[\\)]?$"   # define CSE in Group labels
min_pairs <- 2                                      # Exp35: 2 donors → paired tests OK
B_boot    <- 4000                                   # BCa bootstrap replicates

## ---- Helpers -----------------------------------------------------------------
numish <- function(x) readr::parse_number(as.character(x))
pretty_p <- function(p) ifelse(is.na(p), "n/a",
                               ifelse(p < 1e-4, "p < 1e-4", format.pval(p, digits = 2, eps = 1e-4)))
stars <- function(p) ifelse(is.na(p),"", ifelse(p<.001,"***", ifelse(p<.01,"**", ifelse(p<.05,"*",""))))

theme_pub <- function(){
  theme_classic(base_size = 11) %+replace% theme(
    strip.text = element_text(face="bold"),
    legend.position = "right",
    panel.grid.major.y = element_line(colour="grey90", linewidth=.3),
    panel.grid.minor.y = element_blank()
  )
}
pal_groups <- function() list(
  scale_colour_brewer(palette = "Dark2"),
  scale_fill_brewer(palette = "Dark2")
)

# robust donor id (removes trailing “cse” if it creeps into IDs)
clean_donor <- function(x){
  x <- as.character(x); x <- trimws(x)
  stringr::str_replace(x, "(?i)[ _-]*cse$", "")
}

# Paired Wilcoxon p-value
paired_wilcox <- function(x, y){
  if (length(x) != length(y) || length(x) < 2 || any(is.na(x)) || any(is.na(y))) return(NA_real_)
  out <- try(stats::wilcox.test(x, y, paired = TRUE, exact = FALSE, correct = FALSE), silent = TRUE)
  if (inherits(out, "try-error")) NA_real_ else unname(out$p.value)
}
# BCa CI for the median difference (x - y) on log10 scale
boot_bca <- function(x, y, R = B_boot) {
  eps <- 1e-12
  # Validate R
  if (!is.numeric(R) || length(R) != 1 || is.na(R) || R < 1) return(c(NA_real_, NA_real_))
  R <- as.integer(R)
  # Validate x, y
  if (!is.numeric(x) || !is.numeric(y)) return(c(NA_real_, NA_real_))
  # Remove NA pairs
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) != length(y) || length(x) < 2) return(c(NA_real_, NA_real_))
  n <- length(x)
  d <- x - y
  if (!all(is.finite(d)) || n < 2) return(c(NA_real_, NA_real_))
  # Degenerate: all identical
  if (all(abs(d - d[1]) < eps)) return(rep(stats::median(d, na.rm = TRUE), 2))
  stat <- function(ix) stats::median(d[ix], na.rm = TRUE)
  t0 <- stat(seq_len(n))
  # Jackknife
  jack <- vapply(seq_len(n), function(i) stat(setdiff(seq_len(n), i)), numeric(1))
  jack_var <- sum((mean(jack) - jack)^2)
  if (jack_var < eps) return(rep(t0, 2))
  # Bootstrap
  t_boot <- replicate(R, stat(sample.int(n, n, TRUE)))
    UnitStd      = dplyr::coalesce(
      if ("UnitStd" %in% names(raw)) as.character(UnitStd) else NA_character_,
      if ("Unit"    %in% names(raw)) as.character(Unit)    else NA_character_
    ),  prop <- mean(t_boot < t0)
  prop <- min(max(prop, eps), 1 - eps)
  z0 <- qnorm(prop)
  # aacc: protect denominator
  aacc_num <- sum((mean(jack) - jack)^3)
  aacc_den <- 6 * (jack_var)^(3/2) + eps
  aacc <- aacc_num / aacc_den
  # CI percentiles
  alpha <- c(.025, .975)
  zalpha <- qnorm(alpha)
  pct <- pnorm(z0 + (z0 + zalpha)/(1 - aacc * (z0 + zalpha)))
  pct <- pmin(pmax(pct, eps), 1 - eps) # strictly within (0,1)
  # Quantile
  stats::quantile(t_boot, probs = pct, na.rm = TRUE, names = FALSE)
}

## ---- Load & standardise ------------------------------------------------------
raw <- suppressMessages(readr::read_csv(in_file, show_col_types = FALSE))
stopifnot(all(c("Patient","Group","CompoundName","IC50") %in% names(raw)))

df <- raw %>%
  transmute(
    Patient      = as.character(Patient),
    Group        = as.character(Group),
    CompoundName = as.character(CompoundName),
    UnitStd      = dplyr::coalesce(as.character(UnitStd), as.character(Unit)),
    IC50         = numish(IC50)
  ) %>%
  filter(is.finite(IC50), IC50 > 0) %>%
  mutate(
      diff <- median(pair$base - pair$cse, na.rm = TRUE)
      tibble(n_pairs = np,
             p = p,
             diff_log10 = diff,
             diff_L = ci[1],
             diff_U = ci[2],
             ratio = 10^(-diff),
             ratio_L = 10^(-ci[2]),
             ratio_U = 10^(-ci[1]))    token        = toupper(str_extract(CompoundName, "(?<=_)\\w+$"))        # e.g., "NAC","SFN","ERASTIN" or NA
  )

# Aggregate to one value per donor × arm × compound (median log10(IC50))
df_aggr <- df %>%
  mutate(donor_id = clean_donor(coalesce(Patient, Group))) %>%
  group_by(base_compound, CompoundName, GroupType, UnitStd, donor_id) %>%
  summarise(logIC50 = median(logIC50, na.rm = TRUE), .groups = "drop")

message("Aggregated counts (1 row per donor/arm/compound):")
print(df_aggr %>% count(CompoundName, GroupType) %>% arrange(CompoundName, GroupType))

## ---- (A) CSE vs Base per compound (paired) ----------------------------------
stats_cse <- df_aggr %>%
  distinct(CompoundName, UnitStd) %>%
  group_by(CompoundName, UnitStd) %>%
  group_modify(~{
    sub <- df_aggr %>% filter(CompoundName == .y$CompoundName)
    base_tbl <- sub %>% filter(GroupType == "Base") %>% select(donor_id, base = logIC50)
    cse_tbl  <- sub %>% filter(GroupType == "CSE")  %>% select(donor_id, cse  = logIC50)
    pair <- inner_join(base_tbl, cse_tbl, by = "donor_id")
    np <- nrow(pair)
    if (np >= min_pairs) {
      p  <- paired_wilcox(pair$base, pair$cse)
      ci <- boot_bca(pair$base, pair$cse)
      diff <- median(pair$base - pair$cse)
      tibble(n_pairs=np, p=p, diff_log10=diff, diff_L=ci[1], diff_U=ci[2],
             ratio=10^(-diff), ratio_L=10^(-ci[2]), ratio_U=10^(-ci[1]))
    } else {
      tibble(n_pairs=np, p=NA_real_, diff_log10=NA_real_, diff_L=NA_real_, diff_U=NA_real_,
             ratio=NA_real_, ratio_L=NA_real_, ratio_U=NA_real_)
    }
  }) %>% ungroup()

if (nrow(stats_cse) && any(is.finite(stats_cse$p))) {
  stats_cse <- stats_cse %>% mutate(p_adj_BH = p.adjust(p, method = "BH"),
                                    sig = stars(p_adj_BH))
  readr::write_csv(stats_cse, file.path(out_dir, "Stats_CSE_vs_Base_byCompound.csv"))
}

## ---- (B) Modifier vs base within each arm (paired) --------------------------
base_list <- sort(unique(df_aggr$base_compound))

mk_mod_one <- function(bc){
  mods <- df_aggr %>%
    filter(base_compound == bc) %>%
    transmute(tok = toupper(str_extract(CompoundName, "(?<=_)\\w+$"))) %>%
    filter(!is.na(tok)) %>% distinct() %>% arrange(tok) %>% pull(tok)

  map_dfr(mods, function(tok){
    map_dfr(c("Base","CSE"), function(gtype){
      base_tbl <- df_aggr %>% filter(CompoundName == bc,            GroupType == gtype) %>%
        select(donor_id, base = logIC50, UnitStd) %>% distinct()
      mod_tbl  <- df_aggr %>% filter(CompoundName == paste0(bc,"_",tok), GroupType == gtype) %>%
        select(donor_id, mod = logIC50, UnitStd) %>% distinct()
      pair <- inner_join(base_tbl, mod_tbl, by = "donor_id")
      np <- nrow(pair)
      if (np >= min_pairs) {
        p  <- paired_wilcox(pair$base, pair$mod)
        ci <- boot_bca(pair$base, pair$mod)
        diff <- median(pair$base - pair$mod)
        tibble(
          base_compound = bc, contrast_token = tok, GroupType = gtype,
          UnitStd = paste(unique(c(base_tbl$UnitStd, mod_tbl$UnitStd)), collapse="/"),
          n_pairs = np, p = p, diff_log10 = diff, diff_L = ci[1], diff_U = ci[2],
          ratio = 10^(-diff), ratio_L = 10^(-ci[2]), ratio_U = 10^(-ci[1])
        )
      } else {
        tibble(
          base_compound = bc, contrast_token = tok, GroupType = gtype,
          UnitStd = paste(unique(c(base_tbl$UnitStd, mod_tbl$UnitStd)), collapse="/"),
          n_pairs = np, p = NA_real_, diff_log10 = NA_real_, diff_L = NA_real_, diff_U = NA_real_,
          ratio = NA_real_, ratio_L = NA_real_, ratio_U = NA_real_
        )
      }
    })
  })
}

stats_mod <- map_dfr(base_list, mk_mod_one) %>%
  group_by(base_compound, GroupType) %>%                          # BH within each arm
  mutate(p_adj_BH = ifelse(is.finite(p), p.adjust(p, "BH"), NA_real_)) %>%
  ungroup() %>%
ggsave(file.path(out_dir, "Plot_A_CSE_Effect_byCompound.pdf"),
       p_cse, width = 11, height = 3.8, 
       device = if (capabilities("cairo")) grDevices::cairo_pdf else "pdf")
readr::write_csv(stats_mod, file.path(out_dir, "Stats_Modifiers_vsBase.csv"))

## ---- Plot data (aggregated, cleaner visuals) --------------------------------
df_aggr <- df_aggr %>%
  mutate(
    token   = toupper(stringr::str_extract(CompoundName, "(?<=_)\\w+$")),
    Variant = dplyr::if_else(is.na(token), base_compound, paste0(base_compound, " + ", token))
  ) %>%
  mutate(Variant = factor(Variant, levels = unique(Variant)))

## ---- Plot A: CSE effect per compound (only if CSE present) ------------------
if (any(df_aggr$GroupType == "CSE")) {
  lab_cse <- stats_cse %>%
    transmute(CompoundName,
              label = ifelse(is.finite(p_adj_BH),
                             paste0("FDR: ", pretty_p(p_adj_BH), " ", sig),
                             "FDR: n/a"))
  p_cse <- ggplot(df_aggr, aes(GroupType, logIC50, colour = GroupType, fill = GroupType)) +
    geom_violin(width = 0.9, alpha = 0.15, linewidth = 0.25, colour = NA, trim = FALSE) +
    geom_boxplot(width = 0.55, alpha = 0.25, outlier.shape = NA, linewidth = 0.4) +
    geom_point(size = 2, alpha = 0.85, position = position_jitter(width = 0.06)) +
    geom_line(aes(group = donor_id), alpha = 0.5, colour = "grey40", linewidth = 0.5) +
    facet_wrap(~ CompoundName, scales = "free_y", nrow = 1) +
    pal_groups() +
    scale_y_continuous("log10(IC50)", breaks = pretty_breaks()) +
    labs(x = NULL, title = "CSE effect by compound (paired by donor)") +
    theme_pub() +
    geom_text(data = lab_cse, inherit.aes = FALSE,
              aes(x = -Inf, y = Inf, label = label),
              hjust = -0.03, vjust = 1.1, size = 3)
  ggsave(file.path(out_dir, "Plot_A_CSE_Effect_byCompound.pdf"),
         p_cse, width = 11, height = 3.8, device = cairo_pdf)
}

## ---- Plot B: Modifiers vs base (box/violin + paired points) -----------------
single_base <- length(base_list) == 1L

if (single_base) {
  # precise labels per arm: n pairs, FDR, Δ, CI, ratio
  lab_mod <- stats_mod %>%
    transmute(
      GroupType,
      token = tolower(contrast_token),
      label = ifelse(
        is.finite(p_adj_BH),
        sprintf("n=%d pairs\nFDR: %s  (Δ=%.3f, %.3f–%.3f; ratio=%.2fx)",
                n_pairs, pretty_p(p_adj_BH), diff_log10, diff_L, diff_U, ratio),
        sprintf("n=%d pairs\nFDR: n/a", n_pairs)
      )
    )

  p_mod <- ggplot(df_aggr, aes(Variant, logIC50, colour = Variant, fill = Variant)) +
    geom_violin(width = 0.9, alpha = 0.15, linewidth = 0.25, colour = NA, trim = FALSE) +
    geom_boxplot(width = 0.55, alpha = 0.25, outlier.shape = NA, linewidth = 0.4) +
    geom_point(size = 2, alpha = 0.85, position = position_jitter(width = 0.06)) +
    geom_line(aes(group = donor_id), alpha = 0.5, colour = "grey40", linewidth = 0.5) +
    facet_grid(GroupType ~ token, labeller = labeller(token = function(x) paste0("+", toupper(x))),
               scales = "free_y") +
    pal_groups() +
    scale_y_continuous("log10(IC50)", breaks = pretty_breaks()) +
    labs(x = NULL, title = paste0("Modifier effect vs ", base_list, " (paired by donor within arm)")) +
    theme_pub() +
    geom_text(data = lab_mod, inherit.aes = FALSE,
              aes(x = -Inf, y = Inf, label = label),
              hjust = -0.03, vjust = 1.1, size = 3)

} else {

  lab_mod <- stats_mod %>%
    transmute(
      GroupType, base_compound,
      label = ifelse(
        is.finite(p_adj_BH),
        sprintf("n=%d pairs\nFDR: %s  (Δ=%.3f, %.3f–%.3f; ratio=%.2fx)",
                n_pairs, pretty_p(p_adj_BH), diff_log10, diff_L, diff_U, ratio),
        sprintf("n=%d pairs\nFDR: n/a", n_pairs)
      )
    )

  p_mod <- ggplot(df_aggr, aes(Variant, logIC50, colour = Variant, fill = Variant)) +
    geom_violin(width = 0.9, alpha = 0.15, linewidth = 0.25, colour = NA, trim = FALSE) +
    geom_boxplot(width = 0.55, alpha = 0.25, outlier.shape = NA, linewidth = 0.4) +
    geom_point(size = 2, alpha = 0.85, position = position_jitter(width = 0.06)) +
    geom_line(aes(group = donor_id), alpha = 0.5, colour = "grey40", linewidth = 0.5) +
    facet_grid(GroupType ~ base_compound, scales = "free_y") +
    pal_groups() +
    scale_y_continuous("log10(IC50)", breaks = pretty_breaks()) +
    labs(x = NULL, title = "Modifier effects vs base compound (paired by donor within arm)") +
    theme_pub() +
    geom_text(data = lab_mod, inherit.aes = FALSE,
              aes(x = -Inf, y = Inf, label = label),
              hjust = -0.03, vjust = 1.1, size = 3)
}

ggsave(file.path(out_dir, "Plot_B_Modifier_Effects.pdf"),
       p_mod, width = 11, height = 6.5, device = cairo_pdf)

## ---- Plot C: Paired spaghetti (base vs each modifier, within arm) -----------
# Build explicit base-vs-variant pairs per donor/token/arm → two points & one line each
tokens <- df_aggr %>% filter(!is.na(token)) %>% distinct(base_compound, token) %>% arrange(base_compound, token)

make_spaghetti_df <- function(bc, tok, gtype){
  base_tbl <- df_aggr %>% filter(CompoundName == bc, GroupType == gtype) %>%
    select(base_compound, GroupType, donor_id, base = logIC50)
  mod_tbl  <- df_aggr %>% filter(CompoundName == paste0(bc, "_", tok), GroupType == gtype) %>%
    select(base_compound, GroupType, donor_id, mod  = logIC50)
  pair <- inner_join(base_tbl, mod_tbl, by = c("base_compound","GroupType","donor_id"))
  if (!nrow(pair)) return(NULL)
  long <- bind_rows(
    pair %>% transmute(base_compound, GroupType, token = tok, donor_id,
                       Variant = base_compound, logIC50 = base),
    pair %>% transmute(base_compound, GroupType, token = tok, donor_id,
                       Variant = paste0(base_compound, " + ", tok), logIC50 = mod)
  )
  long
}

spag_df <- pmap_dfr(list(tokens$base_compound, tokens$token), function(bc, tok){
  bind_rows(make_spaghetti_df(bc, tok, "Base"),
            make_spaghetti_df(bc, tok, "CSE"))
})

if (nrow(spag_df)) {
  spag_df <- spag_df %>% mutate(Variant = factor(Variant, levels = unique(Variant)))
  p_spag <- ggplot(spag_df, aes(Variant, logIC50, group = donor_id)) +
    geom_line(alpha = 0.45, colour = "grey50", linewidth = 0.6) +
    geom_point(aes(colour = Variant, fill = Variant), size = 2) +
    facet_grid(GroupType ~ token, labeller = labeller(token = function(x) paste0("+", toupper(x))),
               scales = "free_y") +
    pal_groups() +
    labs(x = NULL, y = "log10(IC50)", title = "Paired spaghetti: modifiers vs base (within arm)") +
    theme_pub()
  ggsave(file.path(out_dir, "Plot_C_PairedSpaghetti_Modifiers.pdf"),
         p_spag, width = 11, height = 6.5, device = cairo_pdf)
}

## ---- Plot D: Forest of ratios (modifiers) -----------------------------------
sf <- stats_mod %>%
  filter(is.finite(ratio) | is.finite(ratio_L) | is.finite(ratio_U)) %>%
  mutate(y_lab = paste0(base_compound, " + ", contrast_token, " (", GroupType, ")",
                        ifelse(is.na(UnitStd) | UnitStd=="", "", paste0(" [", UnitStd, "]")))) %>%
  mutate(y_lab = factor(y_lab, levels = rev(unique(y_lab))))

if (nrow(sf)) {
  p_forest <- ggplot(sf, aes(x = log10(ratio), y = y_lab)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    geom_errorbarh(aes(xmin = log10(ratio_L), xmax = log10(ratio_U)),
                   height = 0.25, linewidth = 0.6, colour = "black", alpha = 0.85) +
    geom_point(aes(fill = p_adj_BH < 0.05), shape = 21, size = 3) +
    scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "black"), guide = "none") +
    scale_x_continuous("log10 ratio (Variant / Base)",
                       breaks = pretty_breaks(),
                       labels = function(z) paste0(z, "  (", number(10^z, accuracy = 0.01), "×)")) +
    labs(y = NULL, title = "Modifier effect on IC50 (paired within arm)") +
    theme_pub()
  ggsave(file.path(out_dir, "Plot_D_Forest_Modifiers_vsBase.pdf"),
         p_forest, width = 7.2, height = 5.2, device = cairo_pdf)
}

## ---- Summary -----------------------------------------------------------------
message("\n=== EXP35 SUMMARY ===")
message("Rows (raw): ", nrow(raw), "; rows (after QC): ", nrow(df))
message("Donors: ", paste(sort(unique(df_aggr$donor_id)), collapse = ", "))
message("Compounds: ", paste(sort(unique(df$CompoundName)), collapse = ", "))
if (exists("stats_cse")) {
  message("CSE contrasts tested: ", sum(is.finite(stats_cse$p), na.rm=TRUE))
}

# Ensure output directory exists
if (!dir.exists(out_dir)) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
}

# Robust file output with error handling
safe_write_csv <- function(obj, path, ...) {
  tryCatch({
    readr::write_csv(obj, path, ...)
    message("[Success] Wrote CSV: ", path)
  }, error = function(e) {
    warning("[Failure] Could not write CSV to ", path, ": ", conditionMessage(e))
  })
}

safe_ggsave <- function(filename, plot, ...) {
  tryCatch({
    ggsave(filename, plot, ...)
    message("[Success] Saved PDF: ", filename)
  }, error = function(e) {
    warning("[Failure] Could not save PDF to ", filename, ": ", conditionMessage(e))
  })
}

# Example usage for post-summary outputs (add as needed)
if (exists("stats_cse")) {
  safe_write_csv(stats_cse, file.path(out_dir, "Stats_CSE_vs_Base_byCompound.csv"))
}
safe_write_csv(stats_mod, file.path(out_dir, "Stats_Modifiers_vsBase.csv"))

if (exists("p_cse")) {
  safe_ggsave(file.path(out_dir, "Plot_A_CSE_Effect_byCompound.pdf"), p_cse, width = 11, height = 3.8, device = cairo_pdf)
}
if (exists("p_mod")) {
  safe_ggsave(file.path(out_dir, "Plot_B_Modifier_Effects.pdf"), p_mod, width = 11, height = 6.5, device = cairo_pdf)
}
if (exists("p_spag")) {
  safe_ggsave(file.path(out_dir, "Plot_C_PairedSpaghetti_Modifiers.pdf"), p_spag, width = 11, height = 6.5, device = cairo_pdf)
}
if (exists("p_forest")) {
  safe_ggsave(file.path(out_dir, "Plot_D_Forest_Modifiers_vsBase.pdf"), p_forest, width = 7.2, height = 5.2, device = cairo_pdf)
}

message("Modifier contrasts tested: ", sum(is.finite(stats_mod$p), na.rm=TRUE))
message("Saved figures to: ", normalizePath(out_dir))
