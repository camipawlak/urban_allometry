## =============================================================================
## percent_difference_lines.R
##
## Three percent-difference spaghetti panels (inventory vs. USFS predicted
## height) in a row, one per prediction-interval level, sharing a common
## y-axis:
##
##     [ PI lower bound ] [ predicted height ] [ PI upper bound ]
##
##     pct(dbh) = 100 * ( fit_Inventory - fit_USFS ) / fit_USFS
##
## at the matching interval level in each panel (inventory lower vs USFS
## lower, etc.). This is the percent-difference sibling of
## 15_spaghetti_log_ratio_patchwork.R -- same alignment machinery, same panel
## layout, linear percent scale instead of the log ratio.
##
## Above 0 = the inventory model predicts taller trees than the USFS pooled
## equation at that dbh. -50 = inventory predicts half the USFS height.
##
## Inventory curves are interpolated onto each species' USFS dbh grid with
## approx() over the overlapping span, as in 12-16. The three panels share one
## y-range so the spread is directly comparable; only the left panel carries
## the y-axis title and tick labels.
##
## Run with working directory = usfs_comparison/
##
## INPUT   data/_processed/line_df.rds        (written by combine_data.R)
## OUTPUT  _figures/pct_diff_patchwork.png
## =============================================================================

library(tidyverse)
library(patchwork)
source("scripts/r/fig-theme.R")

## ---- settings ---------------------------------------------------------------
## Species drawn bold + labeled. Everything else is a thin grey line.
label_species <- c("Eucalyptus globulus",
                   "Acacia melanoxylon",
                   "Ulmus parvifolia")

## Percent difference is unbounded as fit_usfs -> 0, so a few species can blow
## up at the smallest dbh and squash every other line. Y_CLIP sets the visible
## window with coord_cartesian(), which clips the DRAWING without dropping the
## rows (unlike scale_y_continuous(limits =), which would delete them and break
## the lines). Set to NULL to auto-range on the data instead.
Y_CLIP <- c(-60, 85)

OUT_DIR <- "_figures"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

line_df <- readRDS("data/_processed/line_df.rds")

## ---- interpolate inventory onto the USFS grid, keep all three levels ---------
usfs_l <- line_df |> filter(eq_source == "USFS")
inv_l  <- line_df |> filter(eq_source == "Inventory")

aligned <- usfs_l |>
  nest(usfs = -ScientificName) |>
  inner_join(inv_l |> nest(inv = -ScientificName), by = "ScientificName") |>
  mutate(d = map2(usfs, inv, function(u, v) {
    lo_x <- max(min(u$dbh_inch), min(v$dbh_inch))
    hi_x <- min(max(u$dbh_inch), max(v$dbh_inch))
    u |>
      filter(dbh_inch >= lo_x, dbh_inch <= hi_x) |>
      transmute(dbh_inch,
                fit_usfs = fit,  lo_usfs = lo,  hi_usfs = hi,
                fit_inv  = approx(v$dbh_inch, v$fit, xout = dbh_inch)$y,
                lo_inv   = approx(v$dbh_inch, v$lo,  xout = dbh_inch)$y,
                hi_inv   = approx(v$dbh_inch, v$hi,  xout = dbh_inch)$y)
  })) |>
  select(ScientificName, d) |>
  unnest(d) |>
  mutate(is_labeled = ScientificName %in% label_species)

message("", n_distinct(aligned$ScientificName), " species curves")

missing_labels <- setdiff(label_species, unique(aligned$ScientificName))
if (length(missing_labels))
  warning("label_species not present in line_df: ",
          paste(missing_labels, collapse = ", "))

## ---- the three percent-difference tables, so we can share a y-range ---------
levels_spec <- tibble::tribble(
  ~key,  ~inv,      ~usfs,      ~label,
  "lo",  "lo_inv",  "lo_usfs",  "95% PI lower limit",
  "fit", "fit_inv", "fit_usfs", "Predicted height",
  "hi",  "hi_inv",  "hi_usfs",  "95% PI upper limit")

pct_tbl <- function(inv_col, usfs_col) {
  aligned |>
    filter(.data[[usfs_col]] > 0) |>          # denominator must be positive
    mutate(pct_diff = 100 * (.data[[inv_col]] - .data[[usfs_col]]) /
             .data[[usfs_col]])
}
dfs <- map2(levels_spec$inv, levels_spec$usfs, pct_tbl) |>
  set_names(levels_spec$key)

## common y-window across all three panels
y_rng <- if (is.null(Y_CLIP)) {
  r <- range(map_dbl(dfs, \(d) range(d$pct_diff)))
  r + c(-1, 1) * 0.04 * diff(r)
} else Y_CLIP

## common x-window too, so the panels line up on dbh
x_rng <- range(aligned$dbh_inch)

## ---- panel builder ----------------------------------------------------------
## show_y : left panel keeps the y title + tick labels; the others drop them
build_panel <- function(df, level_label, show_y) {
  
  lab_df <- df |>
    filter(is_labeled) |>
    slice_max(dbh_inch, n = 1, by = ScientificName)
  
  p <- ggplot(df, aes(dbh_inch, pct_diff, group = ScientificName)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
    geom_line(data = \(d) filter(d, !is_labeled),
              alpha = 0.4, linewidth = 0.4, color = "grey25") +
    geom_line(data = \(d) filter(d, is_labeled),
              alpha = 1, linewidth = 1.1, color = "grey10") +
    ggrepel::geom_text_repel(
      data = lab_df, aes(label = ScientificName),
      size = 4, fontface = "italic", max.overlaps = Inf,
       min.segment.length = 0,
      direction = "y", hjust = 0, nudge_x = 1, seed = 1) +
    coord_cartesian(xlim = x_rng, ylim = y_rng, clip = "off") +
    labs(x = "DBH (in)",
         y = "% difference (inventory \u2212 USFS) / USFS",
         title = level_label) +
    theme_minimal() +
    theme(plot.title = element_text(size = 15, face = "plain"))
  
  if (!show_y) {
    p <- p + theme(axis.title.y = element_blank(),
                   axis.text.y  = element_blank(),
                   axis.ticks.y = element_blank())
  }
  p
}

panels <- pmap(list(dfs, levels_spec$label, show_y = c(TRUE, FALSE, FALSE)),
               build_panel)

## ---- assemble: 3 in a row, shared title ------------------------------------
combined <- panels[[1]] + panels[[2]] + panels[[3]] +
  plot_layout(nrow = 1) +
  plot_annotation(
    title    = "Percent difference of inventory vs. USFS predicted heights",
    theme    = theme(plot.title   = element_text(face = "bold", size = 18),
                     axis.text = element_text(size = 20)))

# three panels across the full text width
save_fig(combined, file.path(OUT_DIR, "pct_diff_patchwork.png"), height = 2.6)

message("percent_difference_lines: wrote ",
        file.path(OUT_DIR, "pct_diff_patchwork.png"))