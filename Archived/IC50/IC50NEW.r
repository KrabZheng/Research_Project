# ============================
# Resazurin IC50 pipeline — robust & loop-based
# ============================

suppressPackageStartupMessages({
  library(readr);  library(dplyr);  library(tidyr);  library(stringr)
  library(purrr);  library(ggplot2); library(drc);   library(broom)
  library(glue);   library(scales);  library(tibble); library(rlang)
})
options(readr.show_col_types = FALSE)

# ----------------------------
# 0) User config
# ----------------------------
dir_in             <- "Resazurin_IC50"                 # folder with Resazurin_Exp_* files
out_dir            <- file.path(dir_in, "ic50_outputs")
channel_mode       <- "Fluo"                           # "Fluo" or "AbsMinusBG"
exclude_compounds  <- c("Rotenone")                    # case-insensitive
min_reps_for_fit   <- 5                                # min positive conc points for drc
save_svg           <- TRUE

# ----------------------------
# 1) Utilities
# ----------------------------
`%ni%` <- Negate(`%in%`)
safe_num   <- function(x) suppressWarnings(as.numeric(x))
ensure_dir <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)

.detect_format <- function(tok1) {
  if (str_detect(tok1, "^[A-H][0-9]+_Fluo$")) return("A1_format")            # e.g. A1_Fluo
  if (str_detect(tok1, "^[A-Za-z0-9]+_[A-H]_Fluo$")) return("sample_format") # e.g. P676_A_Fluo
  "unknown"
}

# ----------------------------
# 2) Companion maps (robust to whitespace / case)
# ----------------------------
read_compound_map <- function(path_compound) {
  raw <- readr::read_tsv(path_compound, col_names = FALSE, progress = FALSE, trim_ws = TRUE)
  stopifnot(nrow(raw) >= 12)

  get_row_idx <- function(label) which(tolower(trimws(raw$X1)) == tolower(label))[1]

  # Header with slot names (A1, A2, ...)
  hdr   <- raw %>% dplyr::slice(1) %>% as.character()
  slots <- hdr[-1]; slots <- slots[!is.na(slots) & nzchar(slots)]
  slots <- toupper(trimws(slots))

  # Compound names row (case-insensitive)
  i_comp <- get_row_idx("Compound")
  if (is.na(i_comp)) stop("Compound row not found in: ", basename(path_compound))
  comp_nm <- raw %>% dplyr::slice(i_comp) %>% as.character()
  comp_nm <- comp_nm[-1][seq_along(slots)]

  # Unit row (Unit/Units)
  i_unit <- get_row_idx("Unit"); if (is.na(i_unit)) i_unit <- get_row_idx("Units")
  unit_nm <- if (!is.na(i_unit)) {
    u <- raw %>% dplyr::slice(i_unit) %>% as.character()
    u[-1][seq_along(slots)]
  } else rep(NA_character_, length(slots))

  comp_df <- tibble(
    slot       = slots,
    compound   = unname(comp_nm),
    set        = stringr::str_sub(slots, 1, 1),
    slot_index = suppressWarnings(as.integer(stringr::str_sub(slots, 2, 3)))
  )
  unit_df  <- tibble(slot = slots, unit = unname(unit_nm))

  # Dose table = rows with numeric in X1 (dose index)
  dose_tbl <- raw %>%
    dplyr::mutate(idx = suppressWarnings(as.integer(X1))) %>%
    dplyr::filter(!is.na(idx)) %>%
    dplyr::select(-X1)

  n_dose <- min(8L, nrow(dose_tbl))
  dose_vals <- dose_tbl %>% dplyr::slice(1:n_dose) %>% dplyr::mutate(dplyr::across(everything(), as.character))

  doses_long <- as_tibble(dose_vals)
  names(doses_long) <- paste0("col", seq_along(doses_long))
  doses_long <- dplyr::bind_cols(dose_index = 1:n_dose, doses_long) %>%
    tidyr::pivot_longer(-dose_index, names_to = "col", values_to = "conc_raw",
                        values_transform = list(conc_raw = as.character),
                        values_ptypes    = list(conc_raw = character())) %>%
    dplyr::mutate(col_idx = suppressWarnings(as.integer(stringr::str_remove(.data$col, "col")))) %>%
    dplyr::filter(.data$col_idx <= length(slots)) %>%
    dplyr::mutate(slot = slots[.data$col_idx], conc = safe_num(.data$conc_raw)) %>%
    dplyr::transmute(slot = .data$slot, dose_index = .data$dose_index, conc = .data$conc)

  comp_map <- comp_df %>%
    dplyr::transmute(slot = .data$slot, compound = .data$compound, set = .data$set, slot_index = .data$slot_index) %>%
    dplyr::left_join(unit_df,  by = "slot") %>%
    dplyr::left_join(doses_long, by = "slot")

  list(slots = comp_df, doses = comp_map)
}

# ----------------------------
# 3) Parse raw files (two formats) — from raw text
# ----------------------------
parse_raw_file <- function(path_raw, channel_mode = "Fluo") {
  lines <- readr::read_lines(path_raw)
  lines <- lines[nchar(trimws(lines)) > 0]
  head_tokens <- stringr::str_split(lines[1], "\\t+")[[1]]
  fmt <- .detect_format(head_tokens[1])
  if (fmt == "unknown") stop(glue::glue("Unknown raw format in {basename(path_raw)}"))

  ncols_block <- 13L  # 1 label + 12 values for each channel

  read_block <- function(h) {
    mats <- list(fluo = NULL, abs = NULL, bg = NULL)
    for (r in (h + 1):(h + 8)) {
      row_toks <- stringr::str_split(lines[r], "\\t+")[[1]]
      fluo <- safe_num(row_toks[2:(1+12)])
      abs  <- safe_num(row_toks[(2+ncols_block):(1+2*ncols_block)])
      bg   <- safe_num(row_toks[(2+2*ncols_block):(1+3*ncols_block)])
      mats$fluo <- rbind(mats$fluo, fluo)
      mats$abs  <- rbind(mats$abs,  abs)
      mats$bg   <- rbind(mats$bg,   bg)
    }
    mats
  }

  out <- list()

  if (fmt == "A1_format") {
    hdr_idx <- which(stringr::str_detect(lines, "^[A-H][0-9]+_Fluo(\\t|$)"))
    for (h in hdr_idx) {
      lab <- stringr::str_match(lines[h], "^([A-H])(\\d+)_Fluo")[,2:3]
      set_letter  <- lab[1]
      slot_index  <- suppressWarnings(as.integer(lab[2]))
      mats <- read_block(h)
      sig <- if (identical(channel_mode,"AbsMinusBG")) pmax(0, mats$abs - mats$bg) else mats$fluo

      sig_df <- tibble::as_tibble(sig, .name_repair = "minimal")
      colnames(sig_df) <- paste0("C", seq_len(ncol(sig_df)))

      df <- sig_df %>%
        dplyr::mutate(dose_index = dplyr::row_number()) %>%
        tidyr::pivot_longer(dplyr::starts_with("C"), names_to = "col", values_to = "signal") %>%
        dplyr::mutate(
          col        = suppressWarnings(as.integer(stringr::str_remove(.data$col, "^C"))),
          set        = set_letter,
          slot_index = slot_index,
          plate_id   = tools::file_path_sans_ext(basename(path_raw)),
          sample     = NA_character_,
          fmt        = "A1"
        )
      out[[length(out)+1]] <- df
    }
  } else { # sample_format
    hdr_idx <- which(stringr::str_detect(lines, "^[A-Za-z0-9]+_[A-H]_Fluo(\\t|$)"))
    for (h in hdr_idx) {
      lab <- stringr::str_match(lines[h], "^([A-Za-z0-9]+)_([A-H])_Fluo")[,2:3]
      sample_id  <- lab[1]
      set_letter <- lab[2]
      mats <- read_block(h)
      sig <- if (identical(channel_mode,"AbsMinusBG")) pmax(0, mats$abs - mats$bg) else mats$fluo

      sig_df <- tibble::as_tibble(sig, .name_repair = "minimal")
      colnames(sig_df) <- paste0("C", seq_len(ncol(sig_df)))

      df <- sig_df %>%
        dplyr::mutate(dose_index = dplyr::row_number()) %>%
        tidyr::pivot_longer(dplyr::starts_with("C"), names_to = "col", values_to = "signal") %>%
        dplyr::mutate(
          col        = suppressWarnings(as.integer(stringr::str_remove(.data$col, "^C"))),
          set        = set_letter,
          slot_index = NA_integer_,  # compute later from compound map
          plate_id   = tools::file_path_sans_ext(basename(path_raw)),
          sample     = sample_id,
          fmt        = "sample"
        )
      out[[length(out)+1]] <- df
    }
  }

  dplyr::bind_rows(out)
}

# ----------------------------
# 4) Sample map reader (robust)
# ----------------------------
read_sample_map <- function(path_sample) {
  raw <- readr::read_tsv(path_sample, col_names = FALSE, progress = FALSE, trim_ws = TRUE)
  stopifnot(nrow(raw) >= 6)

  header <- raw %>% dplyr::slice(1) %>% unlist(use.names = FALSE)
  stopifnot(header[1] == "ID")
  n_groups <- length(header) - 1L
  idx <- seq_len(n_groups)

  # Group names from a row labeled "Compound" (case-insensitive).
  row_label <- which(tolower(trimws(raw$X1)) == "compound")
  sample_names <- if (length(row_label)) {
    nm <- raw %>% dplyr::slice(row_label[1]) %>% unlist(use.names = FALSE)
    nm[-1][seq_along(idx)]
  } else {
    paste0("Group", idx)
  }

  # Rows 1..12 encode the column-to-group layout (one '1' per column)
  map_rows <- raw %>%
    dplyr::mutate(ID = suppressWarnings(as.integer(X1))) %>%
    dplyr::filter(!is.na(ID), ID >= 1, ID <= 12) %>%
    dplyr::select(-X1)

  mat <- as.matrix(map_rows[, seq_len(n_groups), drop = FALSE])
  col_to_group <- apply(mat, 1, function(v) { g <- which(v == 1); if (length(g) == 1) g else NA_integer_ })

  df <- tibble::tibble(col = 1:12, group = col_to_group) %>%
    dplyr::mutate(sample = ifelse(!is.na(group), sample_names[group], NA_character_)) %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(rep = dplyr::row_number()) %>%
    dplyr::ungroup()

  list(map = df)
}

# ----------------------------
# 5) Attach maps to a raw file (no bare names)
# ----------------------------
read_experiment_from_raw <- function(raw_path) {
  base      <- sub("_Raw\\d*\\.txt$", "", raw_path)
  comp_path <- paste0(base, "_Compound.txt")
  samp_path <- paste0(base, "_Sample.txt")
  stem_name <- basename(base)

  stopifnot(file.exists(raw_path), file.exists(comp_path))

  comp   <- read_compound_map(comp_path)
  df_raw <- parse_raw_file(raw_path, channel_mode = channel_mode)

  # How many slots per set in the map?
  comp_slots <- comp$slots %>%
    dplyr::transmute(slot = .data$slot,
                     compound = .data$compound,
                     set = .data$set,
                     slot_index = .data$slot_index)

  slots_per_set <- comp_slots %>%
    dplyr::group_by(.data$set) %>%
    dplyr::summarise(n_slots = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(n_slots = pmax(1L, .data$n_slots))

  # Compute slot_index for sample_format rows based on map's n_slots
  df_raw <- df_raw %>%
    dplyr::mutate(set = stringr::str_sub(.data$set, 1, 1)) %>%
    dplyr::left_join(slots_per_set, by = "set") %>%
    dplyr::mutate(
      n_slots   = dplyr::coalesce(.data$n_slots, 1L),
      rep_block = pmax(1L, as.integer(12L / .data$n_slots))
    ) %>%
    dplyr::mutate(
      slot_index = dplyr::if_else(
        is.na(.data$slot_index),
        pmin(pmax(1L, ((.data$col - 1L) %/% .data$rep_block) + 1L), .data$n_slots),
        .data$slot_index
      ),
      slot = toupper(paste0(.data$set, .data$slot_index))
    )

  comp_doses <- comp$doses  %>% dplyr::transmute(slot = .data$slot, dose_index = .data$dose_index, conc = .data$conc, unit = .data$unit)
  comp_names <- comp_slots  %>% dplyr::transmute(slot = .data$slot, compound = .data$compound)

  df <- df_raw %>%
    dplyr::left_join(comp_names, by = "slot") %>%
    dplyr::left_join(comp_doses, by = c("slot", "dose_index"))

  # Guarantee 'compound' column; drop unmapped rows with a warning
  if (!"compound" %in% names(df)) df$compound <- NA_character_
  if (all(is.na(df$compound))) {
    message(stem_name)
    message("Raw slots: ", paste(sort(unique(df_raw$slot)), collapse = ", "))
    message("Map slots: ", paste(sort(unique(comp_slots$slot)), collapse = ", "))
    stop("Compound mapping produced only NA for this plate")
  }
  na_before <- sum(is.na(df$compound))
  if (na_before > 0) {
    warning(glue::glue("{stem_name}: dropping {na_before} rows with unmapped compound slots"))
    df <- df %>% dplyr::filter(!is.na(.data$compound))
  }

  # Attach sample map if present
  if (file.exists(samp_path)) {
    sm <- read_sample_map(samp_path)$map %>% dplyr::rename(sample_map = sample, rep_map = rep)
    df <- df %>% dplyr::left_join(sm, by = "col")

    # Ensure targets exist → coalesce
    if (!"sample"     %in% names(df)) df$sample     <- NA_character_
    if (!"rep"        %in% names(df)) df$rep        <- NA_integer_
    if (!"sample_map" %in% names(df)) df$sample_map <- NA_character_
    if (!"rep_map"    %in% names(df)) df$rep_map    <- NA_integer_

    df <- df %>%
      dplyr::mutate(
        sample = dplyr::coalesce(.data$sample, .data$sample_map),
        rep    = dplyr::coalesce(.data$rep,    .data$rep_map)
      ) %>%
      dplyr::select(-sample_map, -rep_map)
  } else {
    df <- df %>% dplyr::mutate(rep = ((.data$col - 1L) %% 3L) + 1L)
  }

  df %>%
    dplyr::mutate(
      experiment = stem_name,
      compound   = as.factor(.data$compound),
      unit       = dplyr::if_else(is.na(.data$unit) | .data$unit == "", NA_character_, .data$unit),
      sample     = as.factor(.data$sample)
    )
}

# ----------------------------
# 6) Discover raw files and read all
# ----------------------------
raw_files <- list.files(dir_in, pattern = "^Resazurin_Exp_\\d+_Raw\\d*\\.txt$", full.names = TRUE)
if (!length(raw_files)) stop("No 'Resazurin_Exp_*_Raw*.txt' files found in dir_in")

message(glue::glue("Found {length(raw_files)} raw file(s):\n- {paste(basename(raw_files), collapse = '\n- ')}"))

ensure_dir(out_dir); ensure_dir(file.path(out_dir, "plots"))

all_df <- purrr::map_dfr(raw_files, read_experiment_from_raw)

# Hard sanity check: 'compound' must exist here
if (!"compound" %in% names(all_df)) stop("Internal: column 'compound' missing in all_df.")

all_df <- all_df %>%
  dplyr::filter(is.finite(.data$signal)) %>%
  dplyr::filter(tolower(as.character(.data$compound)) %ni% tolower(exclude_compounds))

# ----------------------------
# 7) % Viability per plate triplicate
# ----------------------------
all_df <- all_df %>%
  dplyr::mutate(sample_chr  = as.character(.data$sample),
                is_cse      = stringr::str_detect(tolower(.data$sample_chr), "cse$"),
                sample_base = toupper(stringr::str_remove(tolower(.data$sample_chr), "_?cse$")))

viab <- all_df %>%
  dplyr::group_by(dplyr::across(c(experiment, sample, compound, slot))) %>%
  dplyr::mutate(ctrl = median(.data$signal[.data$dose_index == 1], na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(viability = 100 * .data$signal / .data$ctrl) %>%
  dplyr::filter(is.finite(.data$viability)) %>%
  dplyr::mutate(is_cse      = stringr::str_detect(tolower(as.character(.data$sample)), "cse$"),
                sample_base = toupper(stringr::str_remove(tolower(as.character(.data$sample)), "_?cse$")))

# ----------------------------
# 8) Borrowing-as-scale (RSL3/Erastin)
# ----------------------------
anchor_tbl <- dplyr::bind_rows(
  viab %>% dplyr::filter(!.data$is_cse,  as.character(.data$compound) == "RSL3") %>%
    dplyr::group_by(.data$sample_base) %>%
    dplyr::summarise(base = "RSL3", upper = max(.data$viability, na.rm = TRUE), .groups = "drop"),
  viab %>% dplyr::filter( .data$is_cse,  as.character(.data$compound) == "RSL3_CSE") %>%
    dplyr::group_by(.data$sample_base) %>%
    dplyr::summarise(base = "RSL3", lower = min(.data$viability, na.rm = TRUE), .groups = "drop"),
  viab %>% dplyr::filter(!.data$is_cse,  as.character(.data$compound) == "Erastin") %>%
    dplyr::group_by(.data$sample_base) %>%
    dplyr::summarise(base = "Erastin", upper = max(.data$viability, na.rm = TRUE), .groups = "drop"),
  viab %>% dplyr::filter( .data$is_cse,  as.character(.data$compound) == "Erastin_CSE") %>%
    dplyr::group_by(.data$sample_base) %>%
    dplyr::summarise(base = "Erastin", lower = min(.data$viability, na.rm = TRUE), .groups = "drop")
) %>%
  dplyr::group_by(.data$sample_base, .data$base) %>%
  dplyr::summarise(upper = max(.data$upper, na.rm = TRUE),
                   lower = min(.data$lower, na.rm = TRUE), .groups = "drop")

global_anchor <- anchor_tbl %>%
  dplyr::group_by(.data$base) %>%
  dplyr::summarise(global_upper = max(.data$upper, na.rm = TRUE),
                   global_lower = min(.data$lower, na.rm = TRUE), .groups = "drop")

get_anchors <- function(sb, base) {
  row <- anchor_tbl %>% dplyr::filter(.data$sample_base == sb, .data$base == base)
  if (nrow(row) == 0) {
    g <- global_anchor %>% dplyr::filter(.data$base == base)
    return(c(upper = g$global_upper[1], lower = g$global_lower[1]))
  }
  c(upper = row$upper[1], lower = row$lower[1])
}

viab_scaled <- viab %>%
  dplyr::mutate(base_flag = dplyr::case_when(
    as.character(.data$compound) %in% c("RSL3", "Erastin") ~ as.character(.data$compound),
    TRUE ~ NA_character_
  )) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    scaled_viab = {
      if (!is.na(.data$base_flag) && !.data$is_cse) {
        anc <- get_anchors(.data$sample_base, .data$base_flag)
        upp <- anc["upper"]; low <- anc["lower"]
        rng <- ifelse(is.finite(upp - low) && (upp - low) > 1e-8, upp - low, NA_real_)
        if (is.na(rng)) NA_real_ else pmax(0, pmin(100, (.data$viability - low) / rng * 100))
      } else {
        .data$viability
      }
    }
  ) %>%
  dplyr::ungroup()

# ----------------------------
# 9) LL.4 fits per sample × compound (pooled) — safe loops
# ----------------------------
compound_scales <- viab_scaled %>%
  dplyr::group_by(.data$compound) %>%
  dplyr::summarise(unit = dplyr::first(stats::na.omit(.data$unit)), .groups = "drop")

# safe model builder
safe_drm <- function(df) {
  tryCatch(
    suppressWarnings(drm(y ~ conc, data = df, fct = LL.4())),
    error   = function(e) NULL,
    warning = function(w) NULL
  )
}

# all sample×compound combos with any finite conc
combos <- viab_scaled %>%
  dplyr::filter(is.finite(.data$conc)) %>%
  dplyr::distinct(sample, compound) %>%
  dplyr::arrange(.data$compound, .data$sample)

# build fits
fit_rows <- vector("list", nrow(combos))
for (i in seq_len(nrow(combos))) {
  smp <- combos$sample[i]
  cmp <- combos$compound[i]
  d   <- viab_scaled %>% dplyr::filter(.data$sample == smp, .data$compound == cmp)

  df  <- d %>%
    dplyr::transmute(conc = .data$conc, y = .data$scaled_viab) %>%
    dplyr::distinct() %>%
    dplyr::filter(is.finite(conc), conc > 0, is.finite(y))

  npos <- nrow(df)
  m <- NULL
  if (npos >= min_reps_for_fit && length(unique(df$conc)) >= 3) {
    m <- safe_drm(df)  # may be NULL if convergence fails
  }

  fit_rows[[i]] <- tibble(sample = smp, compound = cmp, n = npos, model = list(m))
}
fits <- dplyr::bind_rows(fit_rows)

# robust prediction helper
predict_with_ci <- function(m, newd) {
  if (is.null(m)) return(NULL)
  out <- tryCatch({
    pr  <- predict(m, newdata = newd, se.fit = TRUE)
    fit <- if (is.list(pr) && !is.null(pr$fit)) as.numeric(pr$fit) else as.numeric(pr)
    se  <- if (is.list(pr) && !is.null(pr$se.fit)) as.numeric(pr$se.fit) else rep(NA_real_, length(fit))
    lwr <- fit - 1.96 * se; upr <- fit + 1.96 * se
    dplyr::bind_cols(newd, tibble(fit = fit, lwr = lwr, upr = upr))
  }, error = function(e) NULL, warning = function(w) NULL)
  if (is.null(out) || !all(c("conc","fit") %in% names(out))) return(NULL)
  out
}

# IC50 table
ic50_tbl <- fits %>%
  dplyr::mutate(ED = purrr::map(.data$model, function(m) {
    if (is.null(m)) return(tibble(estimate = NA_real_, lower = NA_real_, upper = NA_real_))
    suppressWarnings({
      ed <- tryCatch(ED(m, 50, interval = "delta", display = FALSE), error = function(e) NA)
      if (all(is.na(ed))) tibble(estimate = NA_real_, lower = NA_real_, upper = NA_real_)
      else tibble(estimate = as.numeric(ed[1]), lower = as.numeric(ed[2]), upper = as.numeric(ed[3]))
    })
  })) %>%
  tidyr::unnest(ED) %>%
  dplyr::left_join(compound_scales, by = "compound") %>%
  dplyr::arrange(.data$compound, .data$sample)

ensure_dir(out_dir)
readr::write_csv(ic50_tbl, file.path(out_dir, "IC50_results.csv"))

# ----------------------------
# 10) Plotting (SVG) with prediction ribbons, conc > 0 only
# ----------------------------
plot_one <- function(smp, cmp) {
  d <- viab_scaled %>% dplyr::filter(.data$sample == smp, .data$compound == cmp)
  d_pos <- d %>% dplyr::filter(is.finite(.data$conc), .data$conc > 0, is.finite(.data$scaled_viab))
  if (!nrow(d_pos)) return(NULL)

  unit_lab <- compound_scales %>% dplyr::filter(.data$compound == cmp) %>% dplyr::pull(.data$unit)
  unit_lab <- ifelse(length(unit_lab) && !is.na(unit_lab), paste0(" (", unit_lab, ")"), "")

  conc_pos <- sort(unique(d_pos$conc))
  if (length(conc_pos) < 2) return(NULL)

  newd <- tibble(conc = seq(min(conc_pos), max(conc_pos), length.out = 200))
  m <- fits %>% dplyr::filter(.data$sample == smp, .data$compound == cmp) %>% dplyr::pull(.data$model)
  m <- if (length(m)) m[[1]] else NULL
  pred <- predict_with_ci(m, newd)

  p <- ggplot(d_pos, aes(x = conc, y = scaled_viab)) +
    geom_point(alpha = 0.65, size = 1.6) +
    stat_summary(fun = median, geom = "line", linewidth = 0.5, alpha = 0.6) +
    scale_y_continuous(labels = scales::label_number(accuracy = 1), limits = c(0, 120)) +
    scale_x_log10(labels = scales::label_number(), breaks = scales::pretty_breaks(6)) +
    labs(title = glue::glue("{cmp} — {as.character(smp)}"),
         x = glue::glue("Concentration{unit_lab}"),
         y = "% viability (scaled)") +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank())

  if (!is.null(pred)) {
    p <- p +
      geom_ribbon(data = pred, aes(y = fit, ymin = lwr, ymax = upr), inherit.aes = FALSE, alpha = 0.18) +
      geom_line(data = pred, aes(y = fit), linewidth = 0.7, alpha = 0.9)
  }
  p
}

save_plot <- function(p, path_base) {
  if (is.null(p)) return(invisible(FALSE))
  if (save_svg) ggsave(paste0(path_base, ".svg"), p, width = 5.5, height = 4.2)
  invisible(TRUE)
}

plot_dir <- file.path(out_dir, "plots"); ensure_dir(plot_dir)

n_ok <- 0L
for (i in seq_len(nrow(combos))) {
  smp <- combos$sample[i]
  cmp <- combos$compound[i]
  gp  <- plot_one(smp, cmp)
  if (inherits(gp, "ggplot")) {
    fn_base <- file.path(plot_dir, sprintf("%s__%s", as.character(cmp), as.character(smp)))
    try(save_plot(gp, fn_base), silent = TRUE)
    n_ok <- n_ok + 1L
  } else {
    message("Skipped (no plottable data): ", as.character(cmp), " — ", as.character(smp))
  }
}
message("Saved ", n_ok, " plot(s) to ", normalizePath(plot_dir, winslash = "/"))

# ----------------------------
# 11) Diagnostics & ranked suggestions
# ----------------------------
fit_diag <- {
  # compute RMSE for models that converged
  rmse_vec <- purrr::map_dbl(fits$model, function(m) {
    if (is.null(m)) return(NA_real_)
    r <- tryCatch(residuals(m), error = function(e) NA)
    if (all(is.na(r))) return(NA_real_)
    sqrt(mean(r^2, na.rm = TRUE))
  })

  base <- fits %>%
    dplyr::mutate(RMSE = rmse_vec,
                  AIC  = purrr::map_dbl(.data$model, ~ tryCatch(AIC(.x), error = function(e) NA_real_)))

  # add simple data diagnostics per group
  diag_dat <- viab_scaled %>%
    dplyr::group_by(.data$sample, .data$compound) %>%
    dplyr::summarise(
      n_points = sum(is.finite(.data$scaled_viab) & is.finite(.data$conc) & .data$conc > 0),
      dose_span = ifelse(sum(.data$conc > 0, na.rm = TRUE) > 1,
                         diff(range(log10(.data$conc[.data$conc > 0]), na.rm = TRUE)), 0),
      monotonic = suppressWarnings(abs(cor(log10(.data$conc[.data$conc > 0]),
                                           .data$scaled_viab[.data$conc > 0],
                                           method = "spearman", use = "complete.obs"))),
      .groups = "drop"
    )

  base %>%
    dplyr::left_join(diag_dat, by = c("sample", "compound")) %>%
    dplyr::mutate(score = 0.5*dose_span +
                         0.4*coalesce(monotonic, 0) +
                         0.2*coalesce(n_points, 0)/10 -
                         0.2*coalesce(RMSE, 0)/20)
}

readr::write_csv(fit_diag, file.path(out_dir, "fit_diagnostics.csv"))

suggested <- fit_diag %>%
  dplyr::filter(n_points >= 6, is.finite(score)) %>%
  dplyr::arrange(dplyr::desc(score)) %>%
  dplyr::select(sample, compound, n_points, dose_span, monotonic, RMSE, AIC, score)

readr::write_csv(suggested, file.path(out_dir, "suggested_combinations.csv"))

message("Outputs directory: ", normalizePath(out_dir, winslash = "/"))
