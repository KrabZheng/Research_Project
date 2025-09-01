## palettes.R  — generic, reusable palette logic ------------------------------

# A colour-blind–safe qualitative base (Okabe–Ito order)
pal_base <- c(
  blue      = "#0072B2",
  orange    = "#E69F00",
  green     = "#009E73",
  vermilion = "#D55E00",
  purple    = "#CC79A7",
  yellow    = "#F0E442",
  sky       = "#56B4E9",
  black     = "#000000"
)

# Get the first n colours (wraps if n > length)
pal_n <- function(n, base = pal_base) {
  unname(base[seq_len(n %% length(base) + (n %% length(base) == 0) * length(base))])
}

# Build a named vector of colours for a set of labels.
# - labels: character or factor values you will map (e.g., levels(df$group))
# - mapping: optional named vector to pin specific labels to specific hex codes
#            e.g., c("control"="#0072B2","treated"="#E69F00")
# - ignore_case: match mapping keys case-insensitively if TRUE
pal_for <- function(labels, mapping = NULL, base = pal_base, ignore_case = TRUE) {
  labs <- as.character(labels)
  uniq <- unique(labs)
  out  <- setNames(rep(NA_character_, length(uniq)), uniq)

  if (!is.null(mapping)) {
    m <- mapping
    if (ignore_case) names(m) <- tolower(names(m))
    key <- if (ignore_case) tolower(uniq) else uniq
    hit_idx <- match(key, names(m))
    has_hit <- !is.na(hit_idx)
    out[has_hit] <- unname(m[hit_idx[has_hit]])
  }

  # Fill the rest from base, avoiding duplicates; repeat base if needed
  pool <- unname(base)
  pool <- pool[!pool %in% out]
  need <- which(is.na(out))
  if (length(need)) {
    fill <- rep_len(pool, length(need))
    out[need] <- fill
  }
  out
}

# ggplot2 helpers: drop-in scales that just need the labels/limits you’re using
scale_colour_by <- function(labels, mapping = NULL, base = pal_base, ignore_case = TRUE, ...) {
  vals <- pal_for(labels, mapping = mapping, base = base, ignore_case = ignore_case)
  ggplot2::scale_colour_manual(values = vals, limits = names(vals), ...)
}
scale_fill_by <- function(labels, mapping = NULL, base = pal_base, ignore_case = TRUE, ...) {
  vals <- pal_for(labels, mapping = mapping, base = base, ignore_case = ignore_case)
  ggplot2::scale_fill_manual(values = vals, limits = names(vals), ...)
}

# (Optional) Make these the session defaults for discrete scales
# NOTE: enable only if you want ALL discrete scales to use pal_base by default.
# scale_colour_discrete <- function(...) ggplot2::scale_colour_manual(values = unname(pal_base), ...)
# scale_fill_discrete   <- function(...) ggplot2::scale_fill_manual(values   = unname(pal_base), ...)