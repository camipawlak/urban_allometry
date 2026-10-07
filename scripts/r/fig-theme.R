## Shared figure styling.
##
## THE PROBLEM THIS SOLVES
##
## Every figure was drawn at base_size 11, but on canvases ranging from 3x3 to
## 16x6 inches. Placed at a common column width, each is rescaled by a different
## factor, so a 3-inch figure enlarges ~2x while a 16-inch one shrinks to ~0.4x --
## its labels then render roughly five times smaller. Nothing about the themes
## caused that; the canvas sizes did.
##
## The fix is to draw each figure at the width it will actually be placed at.
## Text is then consistent across figures without touching base_size.
##
## Heights are per-figure, since aspect ratio is driven by content: a four-panel
## row is genuinely wider than a square scatter.

FIG_WIDTH_FULL <- 7.0    # full text width
FIG_WIDTH_HALF <- 3.5    # single column
FIG_DPI        <- 300
FIG_BASE_SIZE  <- 11

## Palette. Two different blues were in use for what reads as the same thing --
## #1f5c99 for the population curve and #1f78b4 for inventory points. Unified on
## the ColorBrewer Paired blue/red, which the comparison figures already used.
fig_colors <- c(
  inventory  = "#1f78b4",
  usfs       = "#e31a1c",
  population = "#1f78b4",
  species    = "gray70"
)

theme_allometry <- function(base_size = FIG_BASE_SIZE) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.2)
    )
}

## Save at a standard width so the placed text size is consistent.
save_fig <- function(plot, path, width = FIG_WIDTH_FULL, height, ...) {
  ggplot2::ggsave(path, plot, width = width, height = height,
                  dpi = FIG_DPI, ...)
  message("saved ", path, sprintf("  (%.1f x %.1f in)", width, height))
  invisible(path)
}
