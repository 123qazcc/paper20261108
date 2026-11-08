# Measure the actual rendered grid viewports, including faceted ggplots.
# All measurements use the same PDF device, font, and physical size as export.
source("private/render_tools/nature-figure/scripts/panel_alignment.R")
require_rendered_panel_alignment <- function(plot, stem, width_in, height_in,
                                            family = "Helvetica") {
  grob <- if (inherits(plot, "patchwork")) patchwork::patchworkGrob(plot) else ggplot2::ggplotGrob(plot)
  cells <- grob$layout[grepl("^panel(?:-[0-9]+)*$", grob$layout$name, perl = TRUE), , drop = FALSE]
  cells <- cells[order(cells$t, cells$l), , drop = FALSE]
  stopifnot(nrow(cells) > 0L)
  probe <- tempfile(fileext = ".pdf")
  grDevices::pdf(probe, width = width_in, height = height_in, family = family, useDingbats = FALSE)
  is_open <- TRUE
  on.exit({if (is_open) grDevices::dev.off(); unlink(probe)}, add = TRUE)
  grid::grid.newpage(); grid::grid.draw(grob); grid::grid.force()
  panels <- lapply(seq_len(nrow(cells)), function(i) {
    z <- cells[i, ]
    vp <- paste0(z$name, ".", z$t, "-", z$l, "-", z$b, "-", z$r)
    grid::seekViewport(vp)
    lower <- grid::deviceLoc(grid::unit(0, "npc"), grid::unit(0, "npc"), valueOnly = TRUE)
    upper <- grid::deviceLoc(grid::unit(1, "npc"), grid::unit(1, "npc"), valueOnly = TRUE)
    grid::upViewport(0)
    bounds <- c(lower$x, lower$y, upper$x, upper$y) * 72
    stopifnot(all(is.finite(bounds)), bounds[3] > bounds[1], bounds[4] > bounds[2])
    list(id = letters[i], bbox_pt = bounds, grid_id = "rendered-grid",
         row_start = as.integer(z$t - 1L), row_stop = as.integer(z$b),
         col_start = as.integer(z$l - 1L), col_stop = as.integer(z$r))
  })
  grDevices::dev.off(); is_open <- FALSE
  manifest <- list(schema_version = 1L, backend = "r-grid-rendered-viewports",
                   figure = list(width_pt = width_in * 72, height_pt = height_in * 72),
                   panels = panels, exemptions = list())
  jsonlite::write_json(manifest, paste0(stem, ".alignment-layout.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
  python <- "python"
  auditor <- "private/render_tools/nature-figure/scripts/audit_panel_alignment.py"
  status <- system2(python, c(shQuote(auditor), shQuote(paste0(stem, ".alignment-layout.json")),
                             "--json-out", shQuote(paste0(stem, ".alignment.json")),
                             "--overlay-svg", shQuote(paste0(stem, ".alignment.svg")), "--strict"))
  if (status != 0L) stop("Rendered panel alignment failed: ", stem)
  invisible(manifest)
}
