# Presentation-only shared infrastructure for additional English figures.
if (.Platform$OS.type == "windows") Sys.setlocale("LC_CTYPE", "English_United States.utf8")
options(charls.visual.helpers_only = TRUE)
source("R/17i_charls_dynamic_figures_EN.R", encoding = "UTF-8")
export_single <- function(plot, stem, height_mm) {
  path <- file.path(figdir, stem); w <- 183 / 25.4; h <- height_mm / 25.4
  write_json(list(applicable = FALSE, auditable = TRUE, verdict = "NOT_APPLICABLE_SINGLE_PANEL",
    reason = "One plot area; no comparable panel edges", panel_count = 1L),
    paste0(path, ".alignment.json"), pretty = TRUE, auto_unbox = TRUE)
  cairo_pdf(paste0(path, ".pdf"), width = w, height = h, family = font); print(plot); dev.off()
  svglite(paste0(path, ".svg"), width = w, height = h); print(plot); dev.off()
  render_raster(plot, paste0(path, ".png"), w, h, 600)
  render_raster(plot, paste0(path, ".tiff"), w, h, 600, tiff = TRUE)
  render_raster(plot, paste0(path, ".preview.png"), w, h, 300)
  files <- paste0(path, c(".pdf", ".svg", ".png", ".tiff"))
  write_json(list(backend = "R", width_mm = 183, height_mm = height_mm, dpi = 600,
    font = font, files = basename(files), sha256 = as.list(vapply(files, sha, ""))),
    paste0(path, ".export.json"), pretty = TRUE, auto_unbox = TRUE)
}
assert_scientific_unchanged <- function() {
  stopifnot(identical(vapply(input_files, sha, ""), input_hashes),
    identical(unname(vapply(names(frozen), sha, "")), unname(frozen)))
}
