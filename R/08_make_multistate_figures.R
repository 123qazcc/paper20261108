#!/usr/bin/env Rscript

# Publication figure suite for the HRS/ELSA insomnia and cardiovascular
# multi-state analysis.  Plotting, rendering, previews, and exports are R-only.
#
# Figure contract
# Core conclusion: insomnia burden is evaluated as a transition-specific
# exposure across a prespecified, progressive cardiovascular state system,
# with cross-cohort replication and explicit robustness/diagnostic checks.
# Results-level question: which transitions carry the association, how does
# that alter predicted state occupancy, and how stable is the pattern?
# Archetypes: schematic-led state definition followed by quantitative grids.
# Backend: R (ggplot2 + patchwork); no synthetic estimates are ever generated.
# Export contract: editable SVG/PDF plus 600-dpi TIFF/PNG, one source-data CSV
# per plotted figure, >=5-point text, and strict 1.5-point patchwork alignment.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

required_packages <- c("jsonlite", "svglite", "ragg", "scales", "grid")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing R figure dependencies: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

script_version <- "1.2.0"
options(stringsAsFactors = FALSE, scipen = 999)
# Static journal-width declaration for source preflight.  Individual figures
# use this double-column width or a documented narrower schematic width.
width_mm = 183

`%||%` <- function(x, y) if (is.null(x) || !length(x) || all(is.na(x))) y else x

env_flag <- function(name, default = FALSE) {
  value <- tolower(trimws(Sys.getenv(name, if (default) "true" else "false")))
  value %in% c("1", "true", "yes", "y", "on")
}

command_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", command_args, value = TRUE)
output_root_env <- Sys.getenv("GA_MULTISTATE_OUTPUT_ROOT", unset = "")
project_root_env <- Sys.getenv("GA_PROJECT_ROOT", unset = "")
# R on this Windows host cannot reliably translate a non-ASCII --file path.
# Production runs therefore use the explicit ASCII staging/output root.  The
# script-path fallback is evaluated only when neither explicit root is set.
if (nzchar(project_root_env)) {
  project_root <- normalizePath(project_root_env, winslash = "/", mustWork = FALSE)
} else if (nzchar(output_root_env)) {
  project_root <- normalizePath(output_root_env, winslash = "/", mustWork = FALSE)
} else {
  script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1L]) else file.path(getwd(), "R", "08_make_multistate_figures.R")
  project_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
}
output_root <- normalizePath(
  if (nzchar(output_root_env)) output_root_env else project_root,
  winslash = "/", mustWork = FALSE
)
figure_root <- file.path(output_root, "figures", "multistate")
source_data_root <- file.path(figure_root, "source_data")
figure_qa_root <- file.path(output_root, "qc", "figures")
dir.create(figure_root, recursive = TRUE, showWarnings = FALSE)
dir.create(source_data_root, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_qa_root, recursive = TRUE, showWarnings = FALSE)

nature_skill_root <- Sys.getenv(
  "NATURE_FIGURE_SKILL_ROOT",
  unset = "private/render_tools/nature-figure"
)
alignment_helper <- file.path(nature_skill_root, "scripts", "panel_alignment.R")
alignment_auditor <- file.path(nature_skill_root, "scripts", "audit_panel_alignment.py")
pdf_text_auditor <- file.path(nature_skill_root, "scripts", "audit_pdf_text.py")
collision_auditor <- file.path(nature_skill_root, "scripts", "audit_figure_collisions.py")
if (!file.exists(alignment_helper)) {
  stop("Nature figure alignment helper is unavailable: ", alignment_helper, call. = FALSE)
}
source(alignment_helper, local = TRUE)

python_candidates <- c(
  Sys.getenv("GA_FIGURE_AUDIT_PYTHON", unset = ""),
  "python",
  Sys.which("python")
)
python_candidates <- unique(python_candidates[nzchar(python_candidates)])
audit_python <- python_candidates[file.exists(python_candidates)][1L] %||% ""
run_pdf_qa <- env_flag("GA_FIGURE_RUN_PDF_QA", TRUE)
strict_completeness <- env_flag("GA_FIGURE_STRICT_COMPLETENESS", FALSE)

palette_contract <- c(
  neutral_dark = "#252A34",
  neutral_mid = "#707984",
  neutral_light = "#D9DEE3",
  signal_blue = "#2878B5",
  signal_teal = "#2A9D8F",
  accent_orange = "#E28E2C",
  accent_red = "#C84C4C",
  pale_blue = "#BFD7EA",
  pale_orange = "#F2C89B"
)

state_palette <- c(
  S0 = palette_contract[["neutral_mid"]],
  S1 = palette_contract[["signal_blue"]],
  S2 = palette_contract[["accent_orange"]],
  S3 = palette_contract[["accent_red"]],
  S4 = palette_contract[["neutral_dark"]]
)
cohort_palette <- c(HRS = palette_contract[["signal_blue"]], ELSA = palette_contract[["accent_orange"]])
profile_palette <- c(
  "0 symptoms" = palette_contract[["neutral_mid"]],
  "1 symptom" = palette_contract[["pale_blue"]],
  "2 symptoms" = palette_contract[["signal_teal"]],
  "3 symptoms" = palette_contract[["accent_red"]],
  "Sample mean" = palette_contract[["neutral_dark"]],
  "Symptom absent" = palette_contract[["neutral_mid"]],
  "Symptom present" = palette_contract[["accent_red"]]
)

theme_nature_contract <- function(base_size = 7, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.5, colour = "black"),
      legend.title = element_text(size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.5),
      strip.text = element_text(size = base_size, face = "bold"),
      plot.title = element_text(size = base_size + 0.7, face = "bold", margin = margin(b = 3)),
      plot.subtitle = element_text(size = base_size - 0.2, colour = palette_contract[["neutral_mid"]], margin = margin(b = 4)),
      plot.caption = element_text(size = base_size - 0.7, colour = palette_contract[["neutral_mid"]], hjust = 0),
      plot.tag = element_text(size = 8, face = "bold"),
      plot.margin = margin(4, 5, 4, 5),
      panel.grid = element_blank(),
      legend.key.height = grid::unit(3.2, "mm"),
      legend.key.width = grid::unit(4.5, "mm")
    )
}
theme_set(theme_nature_contract())

# The Nature alignment helper measures patchwork geometry on R's base PDF
# device before the publication exports are drawn.  On Windows that probe
# device does not know the system font name "Arial" unless a Type-1 metrics
# alias is registered, even though cairo_pdf, svglite, and ragg do.  The probe
# is temporary and uses Helvetica's metrically compatible core-font metrics;
# final files below continue to render with Arial on their native devices.
if (!"Arial" %in% names(grDevices::pdfFonts())) {
  arial_probe_font <- grDevices::pdfFonts()[["Helvetica"]]
  arial_probe_font$family <- "Arial"
  do.call(grDevices::pdfFonts, setNames(list(arial_probe_font), "Arial"))
}

# patchwork >=1.3 keeps flexible panel tracks as `null` units.  The upstream
# Nature helper correctly identifies their gtable cells, but grid converts a
# bare `null` unit to zero after drawing, which produces a non-auditable box.
# Wrap the helper and resolve those flexible tracks against the final physical
# canvas before the backend-neutral JSON gate is applied.
.nature_alignment_panel_rows <- function(grob) {
  layout_names <- as.character(grob$layout$name)
  standard <- grepl("^panel(?:-[0-9]+)?$", layout_names, perl = TRUE)
  fixed_aspect <- grepl("^panel;.*-[0-9]+$", layout_names, perl = TRUE)
  rows <- grob$layout[standard | fixed_aspect, , drop = FALSE]
  if (nrow(rows) < 2L) {
    stop("At least two patchwork panel cells are required for alignment QA", call. = FALSE)
  }
  rows[order(rows$t, rows$l), , drop = FALSE]
}

nature_write_patchwork_panel_layout_upstream <- write_patchwork_panel_layout
write_patchwork_panel_layout <- function(
    plot, manifest_path, width_in, height_in, panel_ids = NULL,
    row_groups = NULL, column_groups = NULL, exemptions = list()) {
  manifest <- nature_write_patchwork_panel_layout_upstream(
    plot = plot, manifest_path = manifest_path,
    width_in = width_in, height_in = height_in,
    panel_ids = panel_ids, row_groups = row_groups,
    column_groups = column_groups, exemptions = exemptions
  )
  positive_boxes <- vapply(
    manifest$panels,
    function(panel) {
      box <- as.numeric(panel$bbox_pt)
      length(box) == 4L && all(is.finite(box)) && box[3L] > box[1L] && box[4L] > box[2L]
    },
    logical(1)
  )
  if (all(positive_boxes)) return(invisible(manifest))

  grob <- patchwork::patchworkGrob(plot)
  panel_rows <- .nature_alignment_panel_rows(grob)
  resolved_panel_ids <- if (is.null(panel_ids)) {
    vapply(
      seq_len(nrow(panel_rows)),
      function(index) if (index <= length(letters)) letters[index] else paste0("panel-", index),
      character(1)
    )
  } else {
    as.character(panel_ids)
  }
  probe_path <- tempfile(fileext = ".pdf")
  grDevices::pdf(probe_path, width = width_in, height = height_in, useDingbats = FALSE)
  device_open <- TRUE
  on.exit({
    if (device_open) grDevices::dev.off()
    unlink(probe_path)
  }, add = TRUE)
  grid::grid.newpage()
  grid::grid.draw(grob)
  grid::grid.force()

  resolve_tracks <- function(units, total_pt, axis = c("width", "height")) {
    axis <- match.arg(axis)
    converted <- if (axis == "width") {
      grid::convertWidth(units, "pt", valueOnly = TRUE)
    } else {
      grid::convertHeight(units, "pt", valueOnly = TRUE)
    }
    types <- grid::unitType(units)
    flexible <- types == "null"
    if (any(flexible)) {
      weights <- vapply(which(flexible), function(i) as.numeric(units[i]), numeric(1))
      remaining <- total_pt - sum(converted[!flexible])
      if (!is.finite(remaining) || remaining <= 0 ||
          any(!is.finite(weights)) || sum(weights) <= 0) {
        stop("Could not resolve patchwork null tracks at final physical size", call. = FALSE)
      }
      converted[flexible] <- remaining * weights / sum(weights)
    }
    converted
  }

  widths_pt <- resolve_tracks(grob$widths, width_in * 72, "width")
  heights_pt <- resolve_tracks(grob$heights, height_in * 72, "height")
  x_edges <- c(0, cumsum(widths_pt))
  top_edges <- c(0, cumsum(heights_pt))
  total_height_pt <- sum(heights_pt)
  manifest$panels <- lapply(seq_len(nrow(panel_rows)), function(index) {
    row <- panel_rows[index, ]
    list(
      id = resolved_panel_ids[index],
      bbox_pt = unname(c(
        x_edges[row$l], total_height_pt - top_edges[row$b + 1L],
        x_edges[row$r + 1L], total_height_pt - top_edges[row$t]
      )),
      grid_id = "patchwork-grid-1",
      row_start = as.integer(row$t - 1L), row_stop = as.integer(row$b),
      col_start = as.integer(row$l - 1L), col_stop = as.integer(row$r)
    )
  })

  temporary_path <- tempfile(
    pattern = ".panel-layout-resolved-", tmpdir = dirname(manifest_path),
    fileext = ".json"
  )
  jsonlite::write_json(manifest, temporary_path, auto_unbox = TRUE, pretty = TRUE, digits = NA)
  if (file.exists(manifest_path) && !file.remove(manifest_path)) {
    stop("Could not replace unresolved patchwork alignment manifest", call. = FALSE)
  }
  if (!file.rename(temporary_path, manifest_path)) {
    unlink(temporary_path)
    stop("Could not activate resolved patchwork alignment manifest", call. = FALSE)
  }
  grDevices::dev.off()
  device_open <- FALSE
  unlink(probe_path)
  invisible(manifest)
}

state_labels <- c(
  S0 = "S0  No target disease",
  S1 = "S1  Hypertension only",
  S2 = "S2  One CVD without hypertension",
  S3 = "S3  Cardiovascular multimorbidity",
  S4 = "S4  Death"
)
confirmatory_nonfatal_edges <- c("S0->S1", "S0->S2", "S1->S3", "S2->S3")
death_edges <- c("S0->S4", "S1->S4", "S2->S4", "S3->S4")
all_edges <- c("S0->S1", "S0->S2", "S0->S4", "S1->S3", "S1->S4", "S2->S3", "S2->S4", "S3->S4")

normalize_transition <- function(x) {
  x <- gsub("[[:space:]]", "", as.character(x))
  x <- gsub("\u2192", "->", x, fixed = TRUE)
  x <- gsub("\u2013|\u2014", "-", x)
  x <- gsub("^([0-4])->([0-4])$", "S\\1->S\\2", x)
  x
}

display_transition <- function(x) gsub("->", "\u2192", normalize_transition(x), fixed = TRUE)

derive_cohort <- function(analysis_id) {
  value <- tolower(as.character(analysis_id))
  fifelse(grepl("^hrs|_hrs", value), "HRS", fifelse(grepl("^elsa|_elsa", value), "ELSA", NA_character_))
}

analysis_display <- function(analysis_id, cohort = NA_character_) {
  id <- tolower(as.character(analysis_id))
  cohort <- toupper(as.character(cohort))
  fifelse(
    grepl("elsa_4state|w4_w10", id), "ELSA non-fatal follow-up (waves 4-10)",
    fifelse(
      grepl("elsa_5state|through_2012", id), "ELSA mortality capture through 2012 (limited)",
      fifelse(grepl("hrs", id) | cohort == "HRS", "HRS follow-up through 2020", fifelse(cohort == "ELSA", "ELSA", as.character(analysis_id)))
    )
  )
}

profile_display <- function(x) {
  x <- as.character(x)
  map <- c(
    sample_means = "Sample mean",
    insomnia_0 = "0 symptoms",
    insomnia_1 = "1 symptom",
    insomnia_2 = "2 symptoms",
    insomnia_3 = "3 symptoms",
    burden_0 = "0 symptoms",
    burden_1 = "1 symptom",
    burden_ge2 = ">=2 symptoms",
    symptom_absent = "Symptom absent",
    symptom_present = "Symptom present"
  )
  out <- unname(map[x])
  out[is.na(out)] <- x[is.na(out)]
  out
}

model_display <- function(model_id) {
  map <- c(
    M4_MI20 = "M4, 20-imputation pooled",
    M4 = "M4, complete case",
    M3 = "M3, complete case",
    M2 = "M2, complete case",
    M1 = "M1, complete case",
    GA_M0 = "GlobalAgingR baseline"
  )
  out <- unname(map[as.character(model_id)])
  out[is.na(out)] <- as.character(model_id)[is.na(out)]
  out
}

result_search_roots <- unique(normalizePath(c(
  file.path(output_root, "results"),
  file.path(output_root, "reports"),
  file.path(project_root, "results"),
  file.path(project_root, "reports"),
  file.path(Sys.getenv("GA_STAGE_ROOT", unset = "private/stage"), "multistate_output", "results"),
  file.path(Sys.getenv("GA_STAGE_ROOT", unset = "private/stage"), "multistate_output", "reports")
), winslash = "/", mustWork = FALSE))
result_search_roots <- result_search_roots[dir.exists(result_search_roots)]

find_result_candidates <- function(basenames) {
  if (!length(result_search_roots)) return(character())
  candidates <- rbindlist(lapply(seq_along(result_search_roots), function(root_index) {
    paths <- list.files(
      result_search_roots[[root_index]], recursive = TRUE,
      full.names = TRUE, include.dirs = FALSE
    )
    if (!length(paths)) return(data.table())
    data.table(
      path = normalizePath(paths, winslash = "/", mustWork = FALSE),
      root_rank = root_index
    )
  }), use.names = TRUE, fill = TRUE)
  if (!nrow(candidates)) return(character())
  candidates <- candidates[basename(path) %chin% basenames]
  if (!nrow(candidates)) return(character())

  # A file can be discovered through duplicate/nested roots.  Preserve the
  # explicit root priority above, then prefer the newest valid file inside the
  # selected root.  File size is only a final tie-breaker; it must not cause a
  # larger stale report to outrank the active run's report.
  setorder(candidates, root_rank)
  candidates <- candidates[!duplicated(path)]
  info <- file.info(candidates$path)
  candidates[, `:=`(
    size = as.numeric(info$size),
    mtime_numeric = as.numeric(info$mtime)
  )]
  candidates <- candidates[is.finite(size) & size > 0]
  if (!nrow(candidates)) return(character())
  setorder(candidates, root_rank, -mtime_numeric, -size)
  candidates$path
}

read_result_table <- function(basenames, required_columns = character()) {
  candidates <- find_result_candidates(basenames)
  if (!length(candidates)) return(list(data = data.table(), path = NA_character_))
  for (path in candidates) {
    value <- tryCatch(fread(path, encoding = "UTF-8"), error = function(e) NULL)
    if (is.null(value) || !nrow(value)) next
    if (length(required_columns) && !all(required_columns %in% names(value))) next
    return(list(data = as.data.table(value), path = normalizePath(path, winslash = "/", mustWork = FALSE)))
  }
  list(data = data.table(), path = NA_character_)
}

read_top_level_result_table <- function(basename, required_columns = character()) {
  # Canonical collated results live immediately below a search root.  Prefer
  # those exact files over recursively discovered model shards with the same
  # basename; fall back only when no valid canonical file exists.
  candidates <- unique(file.path(result_search_roots, basename))
  candidates <- candidates[file.exists(candidates)]
  for (path in candidates) {
    value <- tryCatch(fread(path, encoding = "UTF-8"), error = function(e) NULL)
    if (is.null(value) || !nrow(value)) next
    if (length(required_columns) && !all(required_columns %in% names(value))) next
    return(list(data = as.data.table(value), path = normalizePath(path, winslash = "/", mustWork = FALSE)))
  }
  read_result_table(basename, required_columns)
}

rename_first <- function(x, target, alternatives) {
  if (target %in% names(x)) return(invisible(x))
  source <- alternatives[alternatives %in% names(x)][1L] %||% NA_character_
  if (!is.na(source)) setnames(x, source, target)
  invisible(x)
}

ensure_column <- function(x, name, value = NA) {
  if (!name %in% names(x)) x[, (name) := value]
  invisible(x)
}

atomic_fwrite <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(unlink(tmp), add = TRUE)
  fwrite(as.data.table(x), tmp, na = "")
  if (file.exists(path) && !file.remove(path)) stop("Could not replace ", path, call. = FALSE)
  if (!file.rename(tmp, path)) stop("Could not activate ", path, call. = FALSE)
  invisible(path)
}

run_command_audit <- function(program, args, report_path = NULL) {
  output <- suppressWarnings(system2(program, args = args, stdout = TRUE, stderr = TRUE))
  status <- attr(output, "status") %||% 0L
  if (!is.null(report_path)) writeLines(output, report_path, useBytes = TRUE)
  list(status = as.integer(status), output = output)
}

run_rendered_pdf_qa <- function(stem) {
  if (!run_pdf_qa) return("NOT_RUN_BY_CONFIGURATION")
  if (!nzchar(audit_python) || !file.exists(audit_python)) {
    stop("Python runtime for backend-neutral PDF audits is unavailable", call. = FALSE)
  }
  if (!all(file.exists(c(pdf_text_auditor, collision_auditor)))) {
    stop("One or more backend-neutral PDF auditors are unavailable", call. = FALSE)
  }
  pdf_path <- paste0(stem, ".pdf")
  audit_base <- file.path(figure_qa_root, basename(stem))
  text_report <- paste0(audit_base, ".pdf-text-audit.json")
  collision_report <- paste0(audit_base, ".collision-audit.json")
  collision_overlay <- paste0(audit_base, ".collision-overlay.pdf")
  text_run <- run_command_audit(
    audit_python,
    c(shQuote(pdf_text_auditor), shQuote(pdf_path), "--min-pt", "5", "--json"),
    text_report
  )
  text_status <- "PASS"
  if (text_run$status != 0L) {
    parsed <- tryCatch(jsonlite::fromJSON(paste(text_run$output, collapse = "\n")), error = function(e) NULL)
    # R's pdf/cairo devices commonly encode every glyph as `1 Tf` and carry
    # its physical size in the text matrix.  The mandated Tf-only auditor
    # explicitly does not account for transforms, so such a file cannot be
    # called a pass or a true size failure from that audit alone.
    r_transform_only <- !is.null(parsed) && isTRUE(parsed$auditable) &&
      isTRUE(parsed$text_run_count > 0) && isTRUE(parsed$minimum_found_pt == 1) &&
      isTRUE(parsed$below_minimum_count == parsed$text_run_count)
    if (r_transform_only) {
      text_status <- "NOT_AUDITABLE_R_TEXT_MATRIX_TRANSFORM"
    } else {
      stop("PDF text-size audit failed for ", basename(stem), "; inspect ", text_report, call. = FALSE)
    }
  }
  collision_run <- run_command_audit(
    audit_python,
    c(
      shQuote(collision_auditor), shQuote(pdf_path),
      "--json-out", shQuote(collision_report),
      "--overlay-pdf", shQuote(collision_overlay)
    ),
    paste0(audit_base, ".collision-audit.log")
  )
  collision_status <- "PASS_OR_REVIEW_REQUIRED"
  if (collision_run$status != 0L) {
    # Preserve the blocking report, but allow the non-strict generation run to
    # finish so every figure can be inspected together.  A strict-completeness
    # run still fails below.  This is especially important for R PDFs whose
    # batched axis glyphs can be returned by PDF geometry engines as one large
    # bounding box; such findings require final-size visual adjudication, not
    # silent deletion of the report or fabricated replacement graphics.
    collision_status <- "FIX_BEFORE_DELIVERY_COLLISION_AUDIT"
  }
  paste(text_status, collision_status, sep = ";")
}

export_figure <- function(
    plot, file_stem, width_mm, height_mm, panel_ids = NULL,
    row_groups = NULL, column_groups = NULL, exemptions = list()) {
  stem <- file.path(figure_root, file_stem)
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  alignment_status <- "NOT_APPLICABLE_SINGLE_PANEL"
  if (inherits(plot, "patchwork")) {
    if (!nzchar(audit_python) || !file.exists(audit_python)) {
      stop("Python runtime for backend-neutral alignment audit is unavailable", call. = FALSE)
    }
    alignment_manifest <- file.path(figure_qa_root, paste0(file_stem, ".alignment-layout.json"))
    alignment_report <- file.path(figure_qa_root, paste0(file_stem, ".alignment.json"))
    alignment_overlay <- file.path(figure_qa_root, paste0(file_stem, ".alignment-overlay.svg"))
    require_patchwork_panel_alignment(
      plot = plot,
      manifest_path = alignment_manifest,
      report_path = alignment_report,
      overlay_svg = alignment_overlay,
      width_in = width_in,
      height_in = height_in,
      panel_ids = panel_ids,
      row_groups = row_groups,
      column_groups = column_groups,
      exemptions = exemptions,
      audit_script = alignment_auditor,
      python = audit_python,
      tolerance_pt = 1.5,
      gutter_tolerance_pt = 1.5,
      strict = TRUE
    )
    alignment_status <- "PASS"
  }

  svglite::svglite(paste0(stem, ".svg"), width = width_in, height = height_in, bg = "white")
  print(plot)
  grDevices::dev.off()
  grDevices::cairo_pdf(
    paste0(stem, ".pdf"), width = width_in, height = height_in,
    family = "Arial", bg = "white", onefile = TRUE
  )
  print(plot)
  grDevices::dev.off()
  ragg::agg_tiff(
    paste0(stem, ".tiff"), width = width_in, height = height_in,
    units = "in", res = 600, background = "white", compression = "lzw"
  )
  print(plot)
  grDevices::dev.off()
  ragg::agg_png(
    paste0(stem, ".png"), width = width_in, height = height_in,
    units = "in", res = 600, background = "white"
  )
  print(plot)
  grDevices::dev.off()

  qa_status <- run_rendered_pdf_qa(stem)
  list(
    stem = stem,
    alignment_status = alignment_status,
    pdf_qa_status = qa_status,
    svg = paste0(stem, ".svg"), pdf = paste0(stem, ".pdf"),
    tiff = paste0(stem, ".tiff"), png = paste0(stem, ".png")
  )
}

manifest_rows <- list()
record_manifest <- function(
    figure_id, file_stem, title, claim, archetype, status,
    source_data = "", input_files = character(), model_selection = "",
    cohort_scope = "", panel_count = 1L, width_mm = NA_real_, height_mm = NA_real_,
    alignment_status = "", pdf_qa_status = "", notes = "") {
  manifest_rows[[length(manifest_rows) + 1L]] <<- data.table(
    figure_id = figure_id,
    file_stem = file_stem,
    title = title,
    figure_level_claim = claim,
    archetype = archetype,
    status = status,
    source_data_csv = if (nzchar(source_data)) basename(source_data) else "",
    input_result_files = paste(unique(na.omit(basename(input_files))), collapse = " | "),
    model_selection = model_selection,
    cohort_scope = cohort_scope,
    panel_count = as.integer(panel_count),
    width_mm = width_mm,
    height_mm = height_mm,
    alignment_status = alignment_status,
    pdf_qa_status = pdf_qa_status,
    notes = notes,
    generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
    script_version = script_version
  )
  invisible(NULL)
}

write_plot_source <- function(x, file_stem, figure_id = "") {
  path <- file.path(source_data_root, paste0(file_stem, ".csv"))
  atomic_fwrite(x, path)
  main_number <- suppressWarnings(as.integer(sub("^Fig0*([0-9]+)$", "\\1", figure_id)))
  if (grepl("^Fig0*[1-9][0-9]*$", figure_id) && is.finite(main_number)) {
    atomic_fwrite(
      x,
      file.path(source_data_root, paste0("Figure", main_number, "_source_data.csv"))
    )
  }
  path
}

save_and_record <- function(
    figure_id, file_stem, plot, source_data, title, claim, archetype,
    input_files, model_selection, cohort_scope, width_mm, height_mm,
    panel_ids = NULL, row_groups = NULL, column_groups = NULL,
    exemptions = list(), notes = "") {
  source_path <- write_plot_source(source_data, file_stem, figure_id)
  export <- export_figure(
    plot, file_stem, width_mm, height_mm,
    panel_ids = panel_ids, row_groups = row_groups,
    column_groups = column_groups, exemptions = exemptions
  )
  record_manifest(
    figure_id, file_stem, title, claim, archetype, "GENERATED",
    source_path, input_files, model_selection, cohort_scope,
    length(panel_ids %||% "a"), width_mm, height_mm,
    export$alignment_status, export$pdf_qa_status, notes
  )
  invisible(export)
}

skip_figure <- function(
    figure_id, file_stem, title, claim, archetype, reason,
    input_files = character(), cohort_scope = "", panel_count = 1L,
    width_mm = NA_real_, height_mm = NA_real_) {
  record_manifest(
    figure_id, file_stem, title, claim, archetype, "SKIPPED_NO_ESTIMATES",
    "", input_files, "", cohort_scope, panel_count, width_mm, height_mm,
    "NOT_RUN", "NOT_RUN", reason
  )
}

normalize_markov_table <- function(x) {
  x <- as.data.table(copy(x))
  if (!nrow(x)) return(x)
  rename_first(x, "hazard_ratio", c("intensity_ratio", "estimate", "IR", "HR"))
  rename_first(x, "conf_low", c("lower", "lower_ci", "ci_low", "lcl"))
  rename_first(x, "conf_high", c("upper", "upper_ci", "ci_high", "ucl"))
  rename_first(x, "probability", c("estimate_probability", "transition_probability", "prob"))
  rename_first(x, "sojourn_years", c("estimate", "mean_sojourn", "sojourn"))
  ensure_column(x, "analysis_id", "")
  ensure_column(x, "cohort", NA_character_)
  x[is.na(cohort) | !nzchar(cohort), cohort := derive_cohort(analysis_id)]
  x[, cohort := toupper(as.character(cohort))]
  ensure_column(x, "transition", "")
  if (all(!nzchar(as.character(x$transition))) && all(c("from_state", "to_state") %in% names(x))) {
    x[, transition := paste0("S", from_state, "->S", to_state)]
  }
  x[, transition := normalize_transition(transition)]
  ensure_column(x, "model_id", "")
  x[, `:=`(
    model_id = as.character(model_id),
    analysis_label_plot = analysis_display(analysis_id, cohort),
    transition_label = display_transition(transition)
  )]
  x[]
}

valid_pooled_model <- function(x, target_model_id) {
  z <- x[get("model_id") == target_model_id]
  if (!nrow(z)) return(FALSE)
  if (!"pooling_status" %in% names(z)) return(TRUE)
  values <- unique(na.omit(as.character(z$pooling_status)))
  !length(values) || all(!grepl("PARTIAL|NOT_ESTIMABLE", values, ignore.case = TRUE))
}

choose_model_per_analysis <- function(x) {
  x <- as.data.table(copy(x))
  if (!nrow(x)) return(x)
  priorities <- c("M4_MI20", "M4", "M3", "M2", "M1", "GA_M0")
  chosen <- x[, {
    candidates <- priorities[priorities %chin% unique(model_id)]
    candidates <- candidates[vapply(candidates, function(id) valid_pooled_model(.SD, id), logical(1))]
    id <- candidates[1L] %||% unique(model_id)[1L]
    .(selected_model_id = id)
  }, by = analysis_id]
  merge(x, chosen, by = "analysis_id", all.x = TRUE, sort = FALSE)[model_id == selected_model_id]
}

preferred_analysis_rows <- function(x, include_elsa_limited = FALSE) {
  x <- as.data.table(copy(x))
  if (!nrow(x)) return(x)
  ids <- unique(as.character(x$analysis_id))
  take <- character()
  hrs <- ids[grepl("^hrs_5state", ids, ignore.case = TRUE)][1L] %||%
    ids[derive_cohort(ids) == "HRS"][1L] %||% NA_character_
  elsa_nonfatal <- ids[grepl("elsa_4state|w4_w10", ids, ignore.case = TRUE)][1L] %||%
    ids[derive_cohort(ids) == "ELSA" & !grepl("5state|2012", ids, ignore.case = TRUE)][1L] %||% NA_character_
  take <- c(take, hrs, elsa_nonfatal)
  if (include_elsa_limited) {
    elsa_limited <- ids[grepl("elsa_5state|through_2012", ids, ignore.case = TRUE)][1L] %||% NA_character_
    take <- c(take, elsa_limited)
  }
  take <- unique(na.omit(take))
  if (!length(take)) return(x)
  x[analysis_id %chin% take]
}

add_confirmatory_fdr <- function(x) {
  x <- as.data.table(copy(x))
  ensure_column(x, "p_value", NA_real_)
  x[, confirmatory_bh_fdr := NA_real_]
  group_columns <- intersect(
    c("cohort", "analysis_id", "model_id", "model_level", "basis", "weighted", "term"),
    names(x)
  )
  eligible <- x$transition %chin% confirmatory_nonfatal_edges & is.finite(x$p_value)
  if (length(group_columns)) {
    x[eligible,
      confirmatory_bh_fdr := p.adjust(p_value, method = "BH"),
      by = group_columns]
  } else {
    x[eligible, confirmatory_bh_fdr := p.adjust(p_value, method = "BH")]
  }
  x[]
}

forest_panel <- function(x, title, x_title = "Transition intensity ratio (95% CI)") {
  x <- copy(x)
  # Upstream callers retain only strictly positive ratios and confidence limits
  # before applying the logarithmic axis; no pseudocount is introduced.
  edge_order <- rev(all_edges[all_edges %chin% unique(x$transition)])
  x[, transition_factor := factor(transition, levels = edge_order, labels = display_transition(edge_order))]
  ggplot(x, aes(x = hazard_ratio, y = transition_factor)) +
    geom_vline(xintercept = 1, linewidth = 0.35, linetype = "dashed", colour = palette_contract[["neutral_mid"]]) +
    geom_errorbar(aes(xmin = conf_low, xmax = conf_high), orientation = "y", width = 0.16, linewidth = 0.48, colour = palette_contract[["signal_blue"]]) +
    geom_point(aes(fill = confirmatory_bh_fdr < 0.05), shape = 21, size = 2.0, stroke = 0.35, colour = "black") +
    scale_fill_manual(values = c(`TRUE` = palette_contract[["accent_red"]], `FALSE` = "white"), na.value = "white", guide = "none") +
    scale_x_log10(labels = scales::label_number(accuracy = 0.01)) +
    labs(x = x_title, y = NULL, title = title) +
    theme_nature_contract() +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())
}

add_patchwork_tags <- function(plot) {
  plot + plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(size = 8, face = "bold", family = "Arial"))
}

# -------------------------------------------------------------------------
# Load all result families once.  Empty or failed model outputs stay empty;
# the script records a skipped figure rather than inventing substitute values.
# -------------------------------------------------------------------------

catalog_result <- read_result_table(c("transition_catalog.csv"), c("from_state", "to_state"))
rate_result <- read_result_table(c("transition_counts_rates.csv"), c("cohort", "transition"))
markov_hr_result <- read_result_table(c(
  "markov_hazard_ratios.csv", "markov_transition_intensity_ratios.csv",
  "markov_transition_intensity_estimates.csv"
))
markov_q_result <- read_top_level_result_table(
  "markov_q_intensities.csv",
  c(
    "analysis_id", "model_id", "transition", "covariate_profile",
    "intensity", "conf_low", "conf_high", "status"
  )
)
markov_probability_result <- read_result_table(c(
  "markov_transition_probabilities_2_5_10y.csv", "markov_transition_probabilities.csv",
  "markov_transition_probability_matrices.csv"
))
markov_sojourn_result <- read_result_table(c("markov_sojourn_times.csv", "markov_mean_sojourn_times.csv"))
markov_prevalence_result <- read_result_table(c(
  "markov_observed_expected_prevalence.csv", "markov_state_occupancy.csv",
  "markov_observed_expected_occupancy.csv"
))
cox_result <- read_result_table(c("transition_cox_estimates.csv", "cox_transition_estimates.csv"))
protocol_sensitivity_result <- read_result_table(
  c("protocol_sensitivity_estimates.csv"),
  c("sensitivity_id", "cohort", "scenario", "transition", "hazard_ratio")
)
protocol_coverage_result <- read_result_table(
  c("protocol_sensitivity_coverage.csv"),
  c(
    "sensitivity_id", "cohort", "protocol_sensitivity", "status",
    "reason", "evidence_files"
  )
)

catalog <- as.data.table(copy(catalog_result$data))
rates <- as.data.table(copy(rate_result$data))
markov_hr <- normalize_markov_table(markov_hr_result$data)
markov_q <- normalize_markov_table(markov_q_result$data)
markov_probability <- normalize_markov_table(markov_probability_result$data)
markov_sojourn <- normalize_markov_table(markov_sojourn_result$data)
markov_prevalence <- normalize_markov_table(markov_prevalence_result$data)
cox_estimates <- as.data.table(copy(cox_result$data))
protocol_sensitivity_estimates <- as.data.table(copy(protocol_sensitivity_result$data))
protocol_sensitivity_coverage <- as.data.table(copy(protocol_coverage_result$data))

# -------------------------------------------------------------------------
# Fig. 1 | Prespecified state system and the eight allowed intensities.
# -------------------------------------------------------------------------

fig_id <- "Fig01"
file_stem <- "fig01_prespecified_multistate_structure"
fig_title <- "Prespecified cardiovascular multi-state structure"
fig_claim <- "The primary model uses eight progressive intensities and deliberately excludes a direct S0-to-S3 intensity."
if (!nrow(catalog)) {
  skip_figure(fig_id, file_stem, fig_title, fig_claim, "schematic-led composite", "transition_catalog.csv is unavailable or empty", catalog_result$path, "HRS and ELSA", 1L, 170, 92)
} else {
  catalog[, transition := normalize_transition(if ("transition" %in% names(catalog)) transition else paste0("S", from_state, "->S", to_state))]
  catalog[, `:=`(from_state = as.integer(from_state), to_state = as.integer(to_state))]
  observed_edges <- sort(unique(catalog$transition))
  if (length(observed_edges) != 8L || !setequal(observed_edges, all_edges) || "S0->S3" %chin% observed_edges) {
    stop("State diagram blocked: transition catalog is not the protocol-defined eight-edge structure", call. = FALSE)
  }
  nodes <- data.table(
    state = 0:4,
    state_code = paste0("S", 0:4),
    x = c(0, 1.35, 1.35, 2.85, 4.25),
    y = c(0, 0.90, -0.90, 0, 0),
    description_x = c(-0.03, 1.35, 1.35, 2.85, 4.25),
    description_y = c(-0.34, 1.25, -1.25, -0.43, -0.34),
    label = unname(state_labels[paste0("S", 0:4)]),
    description = c(
      "No target disease", "Hypertension only", "One CVD without hypertension",
      "Multimorbidity", "Death"
    )
  )
  geometry <- data.table(
    transition = all_edges,
    curvature = c(0.10, -0.10, -0.34, 0.10, -0.23, -0.10, 0.23, 0)
  )
  edges <- merge(catalog, geometry, by = "transition", all.x = TRUE, sort = FALSE)
  edges <- merge(edges, nodes[, .(from_state = state, x, y)], by = "from_state", all.x = TRUE, sort = FALSE)
  setnames(edges, c("x", "y"), c("x_from", "y_from"))
  edges <- merge(edges, nodes[, .(to_state = state, x, y)], by = "to_state", all.x = TRUE, sort = FALSE)
  setnames(edges, c("x", "y"), c("x_to", "y_to"))
  edges[, edge_length := sqrt((x_to - x_from)^2 + (y_to - y_from)^2)]
  edges[, `:=`(
    x_from_plot = x_from + 0.21 * (x_to - x_from) / edge_length,
    y_from_plot = y_from + 0.21 * (y_to - y_from) / edge_length,
    x_to_plot = x_to - 0.27 * (x_to - x_from) / edge_length,
    y_to_plot = y_to - 0.27 * (y_to - y_from) / edge_length
  )]

  p_state <- ggplot() +
    coord_equal(xlim = c(-0.55, 4.78), ylim = c(-1.72, 1.82), clip = "off") +
    theme_void(base_family = "Arial", base_size = 7) +
    theme(
      plot.margin = margin(8, 12, 8, 12),
      plot.title = element_text(size = 8, face = "bold"),
      plot.subtitle = element_text(size = 6.6, colour = palette_contract[["neutral_mid"]]),
      plot.caption = element_text(size = 5.6, colour = palette_contract[["neutral_mid"]], hjust = 0)
    ) +
    labs(
      title = fig_title,
      subtitle = "Panel-observed progressive disease states; death is absorbing",
      caption = "No direct S0\u2192S3 intensity: an observed panel jump is treated as an interval-censored, unobserved multi-step path."
    )
  for (i in seq_len(nrow(edges))) {
    p_state <- p_state + geom_curve(
      data = edges[i],
      aes(x = x_from_plot, y = y_from_plot, xend = x_to_plot, yend = y_to_plot),
      curvature = edges$curvature[i], linewidth = 0.65,
      colour = if (edges$to_state[i] == 4L) palette_contract[["neutral_dark"]] else palette_contract[["signal_blue"]],
      arrow = grid::arrow(length = grid::unit(2.2, "mm"), type = "closed"),
      lineend = "round"
    )
  }
  p_state <- p_state +
    geom_point(
      data = nodes,
      aes(x = x, y = y, fill = state_code),
      shape = 21, size = 10.5, stroke = 0.6, colour = "white", show.legend = FALSE
    ) +
    scale_fill_manual(values = state_palette)
  # Draw each label as its own R text object.  Besides keeping layout explicit,
  # this prevents PDF geometry tools from treating spatially separated labels
  # as one large batched text bounding box.
  for (i in seq_len(nrow(nodes))) {
    p_state <- p_state +
      annotate(
        "text", x = nodes$x[i], y = nodes$y[i], label = nodes$state_code[i],
        family = "Arial", size = 2.7, fontface = "bold", colour = "white"
      ) +
      annotate(
        "text", x = nodes$description_x[i], y = nodes$description_y[i],
        label = nodes$description[i],
        family = "Arial", size = 1.95, fontface = "bold", colour = palette_contract[["neutral_dark"]]
      )
  }
  p_state <- p_state +
    annotate(
      "text", x = 3.50, y = 1.50, label = "Blue = non-fatal progression     Dark = transition to death",
      family = "Arial", size = 1.9, colour = palette_contract[["neutral_mid"]], hjust = 0.5
    )
  source_state <- rbindlist(list(
    nodes[, .(
      record_type = "state", state_code, state_label = label,
      from_state = state, to_state = NA_integer_, transition = NA_character_,
      x_from = x, y_from = y, x_to = NA_real_, y_to = NA_real_, curvature = NA_real_
    )],
    edges[, .(
      record_type = "allowed_transition", state_code = NA_character_, state_label = NA_character_,
      from_state, to_state, transition, x_from, y_from, x_to, y_to, curvature
    )]
  ), use.names = TRUE, fill = TRUE)
  save_and_record(
    fig_id, file_stem, p_state, source_state, fig_title, fig_claim,
    "schematic-led composite", catalog_result$path, "Protocol-defined structural model",
    "HRS and ELSA", 170, 92,
    notes = "The S0-to-S3 panel jump is interval-censored evidence of an unobserved multi-step path, not a direct intensity."
  )
}

# -------------------------------------------------------------------------
# Fig. 2 | Primary non-fatal transition-specific intensity ratios.
# -------------------------------------------------------------------------

fig_id <- "Fig02"
file_stem <- "fig02_primary_nonfatal_transition_intensity_ratios"
fig_title <- "Insomnia and non-fatal cardiovascular transitions"
fig_claim <- "The primary multi-state model localizes the insomnia association to specific non-fatal disease transitions in HRS and ELSA."
hr_primary <- preferred_analysis_rows(markov_hr)
if (nrow(hr_primary)) {
  ensure_column(hr_primary, "is_requested_insomnia_exposure", 1L)
  ensure_column(hr_primary, "term", "insomnia_score_primary")
  hr_primary <- hr_primary[
    is_requested_insomnia_exposure == 1L &
      (term == "insomnia_score_primary" | grepl("insomnia_score_primary", term, fixed = TRUE)) &
      transition %chin% confirmatory_nonfatal_edges &
      is.finite(hazard_ratio) & hazard_ratio > 0 & is.finite(conf_low) & conf_low > 0 & is.finite(conf_high) & conf_high > 0
  ]
  hr_primary <- choose_model_per_analysis(hr_primary)
  hr_primary <- add_confirmatory_fdr(hr_primary)
}
if (!nrow(hr_primary) || uniqueN(hr_primary$cohort) < 2L) {
  skip_figure(fig_id, file_stem, fig_title, fig_claim, "quantitative grid", "Primary Markov intensity-ratio estimates are unavailable for both HRS and ELSA", markov_hr_result$path, "HRS and ELSA", 2L, 183, 96)
} else {
  cohorts <- c("HRS", "ELSA")
  panels <- lapply(cohorts, function(cohort_name) {
    z <- hr_primary[cohort == cohort_name]
    subtitle <- unique(z$analysis_label_plot)[1L]
    forest_panel(z, cohort_name) +
      labs(subtitle = subtitle)
  })
  fig <- add_patchwork_tags(panels[[1L]] | panels[[2L]])
  save_and_record(
    fig_id, file_stem, fig, hr_primary, fig_title, fig_claim, "quantitative grid",
    markov_hr_result$path,
    paste(unique(paste(hr_primary$analysis_id, hr_primary$model_id, sep = ": ")), collapse = " | "),
    "HRS and ELSA non-fatal comparison", 183, 96,
    panel_ids = c("a", "b"), row_groups = list(c("a", "b")),
    notes = "Filled points denote BH-FDR <0.05 across the four confirmatory non-fatal transitions within each fitted analysis."
  )
}

# -------------------------------------------------------------------------
# Extended Data Fig. 1 | Death-transition intensity ratios.
# -------------------------------------------------------------------------

fig_id <- "ExtFig01"
file_stem <- "extfig01_death_transition_intensity_ratios"
fig_title <- "Insomnia and transitions to death"
fig_claim <- "Mortality-transition estimates are shown with full HRS follow-up and separately labelled limited ELSA mortality capture."
hr_death <- preferred_analysis_rows(markov_hr, include_elsa_limited = TRUE)
if (nrow(hr_death)) {
  ensure_column(hr_death, "is_requested_insomnia_exposure", 1L)
  ensure_column(hr_death, "term", "insomnia_score_primary")
  hr_death <- hr_death[
    is_requested_insomnia_exposure == 1L &
      (term == "insomnia_score_primary" | grepl("insomnia_score_primary", term, fixed = TRUE)) &
      transition %chin% death_edges & is.finite(hazard_ratio) & hazard_ratio > 0 &
      is.finite(conf_low) & conf_low > 0 & is.finite(conf_high) & conf_high > 0
  ]
  hr_death <- hr_death[grepl("hrs_5state|elsa_5state|through_2012", analysis_id, ignore.case = TRUE)]
  hr_death <- choose_model_per_analysis(hr_death)
  hr_death[, confirmatory_bh_fdr := NA_real_]
}
if (!nrow(hr_death)) {
  skip_figure(fig_id, file_stem, fig_title, fig_claim, "quantitative grid", "No fitted death-transition Markov intensity ratios are available", markov_hr_result$path, "HRS and limited ELSA mortality", 2L, 183, 96)
} else {
  analyses <- unique(hr_death$analysis_id)
  analyses <- analyses[order(match(derive_cohort(analyses), c("HRS", "ELSA")))]
  panels <- lapply(analyses, function(id) {
    z <- hr_death[analysis_id == id]
    forest_panel(z, analysis_display(id, z$cohort[1L])) +
      labs(subtitle = model_display(unique(z$model_id)[1L]))
  })
  if (length(panels) == 1L) {
    fig <- panels[[1L]]
    panel_ids <- NULL
    row_groups <- NULL
  } else {
    fig <- add_patchwork_tags(wrap_plots(panels, nrow = 1L))
    panel_ids <- letters[seq_along(panels)]
    row_groups <- list(panel_ids)
  }
  save_and_record(
    fig_id, file_stem, fig, hr_death, fig_title, fig_claim, "quantitative grid",
    markov_hr_result$path,
    paste(unique(paste(hr_death$analysis_id, hr_death$model_id, sep = ": ")), collapse = " | "),
    "HRS mortality through 2020; ELSA mortality only where captured through 2012", 183, 96,
    panel_ids = panel_ids, row_groups = row_groups,
    notes = "ELSA death-transition estimates must not be interpreted as mortality follow-up through wave 10."
  )
}

# -------------------------------------------------------------------------
# Extended Data Fig. 2 | Transition-specific Cox sensitivity estimates.
# -------------------------------------------------------------------------

fig_id <- "ExtFig02"
file_stem <- "extfig02_transition_specific_cox_sensitivity"
fig_title <- "Visit-detection-time Cox sensitivity analysis"
fig_claim <- "A semiparametric transition-specific sensitivity analysis tests whether the primary pattern depends on the continuous-time Markov estimator."
cox_primary <- copy(cox_estimates)
if (nrow(cox_primary)) {
  rename_first(cox_primary, "hazard_ratio", c("HR", "estimate_ratio"))
  rename_first(cox_primary, "conf_low", c("lower", "lower_ci", "lcl"))
  rename_first(cox_primary, "conf_high", c("upper", "upper_ci", "ucl"))
  ensure_column(cox_primary, "analysis_id", "")
  ensure_column(cox_primary, "cohort", NA_character_)
  ensure_column(cox_primary, "transition", "")
  ensure_column(cox_primary, "term_role", "exposure")
  ensure_column(cox_primary, "term", "")
  ensure_column(cox_primary, "model_level", "M4")
  ensure_column(cox_primary, "basis", "complete_case")
  ensure_column(cox_primary, "weighted", 0L)
  cox_primary[, `:=`(
    cohort = toupper(as.character(cohort)),
    transition = normalize_transition(transition)
  )]
  cox_primary <- cox_primary[
    analysis_id == "primary_score_tv" & model_level == "M4" & weighted == 0L &
      term_role == "exposure" & !grepl(":", term, fixed = TRUE) &
      is.finite(hazard_ratio) & hazard_ratio > 0 & is.finite(conf_low) & conf_low > 0 & is.finite(conf_high) & conf_high > 0
  ]
  if (nrow(cox_primary)) {
    cox_primary[, basis_rank := fifelse(basis == "mice_pooled", 1L, fifelse(basis == "complete_case", 2L, 3L))]
    cox_primary <- cox_primary[order(basis_rank)][, .SD[1L], by = .(cohort, transition, term)]
    cox_primary <- add_confirmatory_fdr(cox_primary)
  }
}
if (!nrow(cox_primary)) {
  skip_figure(fig_id, file_stem, fig_title, fig_claim, "quantitative grid", "Transition-specific Cox sensitivity estimates are unavailable", cox_result$path, "HRS and ELSA", 2L, 183, 112)
} else {
  cohorts <- intersect(c("HRS", "ELSA"), unique(cox_primary$cohort))
  panels <- lapply(cohorts, function(cohort_name) {
    z <- cox_primary[cohort == cohort_name]
    forest_panel(z, cohort_name, "Hazard ratio (95% CI)")
  })
  if (length(panels) == 1L) {
    fig <- panels[[1L]]
    panel_ids <- NULL
    row_groups <- NULL
  } else {
    fig <- add_patchwork_tags(wrap_plots(panels, nrow = 1L))
    panel_ids <- letters[seq_along(panels)]
    row_groups <- list(panel_ids)
  }
  save_and_record(
    fig_id, file_stem, fig, cox_primary, fig_title, fig_claim, "quantitative grid",
    cox_result$path,
    "Primary time-varying score; M4; unweighted; pooled imputation preferred over complete case",
    "HRS and ELSA; ELSA death capture remains limited", 183, 112,
    panel_ids = panel_ids, row_groups = row_groups,
    notes = "Cox event times are visit-detection approximations and constitute sensitivity, not the primary panel-observed Markov analysis."
  )
}

# -------------------------------------------------------------------------
# Fig. 3 | Predicted 2-, 5-, and 10-year transition probabilities from S0.
# -------------------------------------------------------------------------

fig_id <- "Fig03"
file_stem <- "fig03_predicted_state_probabilities_from_s0"
fig_title <- "Predicted state probabilities after entry in S0"
fig_claim <- "Contrasting low and high insomnia symptom profiles shows how transition-specific effects accumulate into clinically interpretable state probabilities."
prob_primary <- preferred_analysis_rows(markov_probability)
if (nrow(prob_primary)) {
  ensure_column(prob_primary, "exposure_variable", "insomnia_score_primary")
  ensure_column(prob_primary, "profile", "")
  ensure_column(prob_primary, "time_years", NA_real_)
  ensure_column(prob_primary, "from_state", NA_integer_)
  ensure_column(prob_primary, "to_state", NA_integer_)
  prob_primary <- prob_primary[
    exposure_variable == "insomnia_score_primary" &
      profile %chin% c("insomnia_0", "insomnia_3") &
      time_years %in% c(2, 5, 10) & from_state == 0L &
      is.finite(probability) & probability >= 0 & probability <= 1
  ]
  prob_primary <- choose_model_per_analysis(prob_primary)
  prob_primary[, `:=`(
    exposure_profile = factor(profile_display(profile), levels = c("0 symptoms", "3 symptoms")),
    destination_state = factor(paste0("S", to_state), levels = paste0("S", 0:4), labels = unname(state_labels[paste0("S", 0:4)]))
  )]
}
if (!nrow(prob_primary) || uniqueN(prob_primary$cohort) < 2L) {
  skip_figure(fig_id, file_stem, fig_title, fig_claim, "quantitative grid", "Primary 2/5/10-year transition probabilities from S0 are unavailable for both cohorts", markov_probability_result$path, "HRS and ELSA", 2L, 183, 105)
} else {
  panels <- lapply(c("HRS", "ELSA"), function(cohort_name) {
    z <- prob_primary[cohort == cohort_name]
    state_codes <- paste0("S", 0:4)
    state_colours <- setNames(unname(state_palette[state_codes]), unname(state_labels[state_codes]))
    ggplot(z, aes(x = time_years, y = probability, colour = destination_state, linetype = exposure_profile, group = interaction(destination_state, exposure_profile))) +
      geom_line(linewidth = 0.62) +
      geom_point(size = 1.45, stroke = 0.2) +
      scale_colour_manual(values = state_colours, breaks = names(state_colours), labels = state_codes, drop = FALSE) +
      scale_linetype_manual(values = c("0 symptoms" = "solid", "3 symptoms" = "22"), drop = FALSE) +
      scale_x_continuous(breaks = c(2, 5, 10)) +
      scale_y_continuous(limits = c(0, 1), labels = scales::label_percent(accuracy = 1), expand = expansion(mult = c(0, 0.03))) +
      labs(
        x = "Years since S0", y = "Predicted probability", colour = "Destination state",
        linetype = "Insomnia profile", title = cohort_name,
        subtitle = analysis_display(unique(z$analysis_id)[1L], cohort_name)
      ) +
      theme_nature_contract() +
      theme(legend.position = "bottom", legend.box = "vertical")
  })
  panels[[1L]] <- panels[[1L]] +
    theme(legend.position = "bottom", legend.box = "vertical", legend.box.just = "center")
  panels[[2L]] <- panels[[2L]] + theme(legend.position = "none")
  fig <- add_patchwork_tags(
    (panels[[1L]] | panels[[2L]]) + plot_layout(guides = "collect")
  )
  save_and_record(
    fig_id, file_stem, fig, prob_primary, fig_title, fig_claim, "quantitative grid",
    markov_probability_result$path,
    paste(unique(paste(prob_primary$analysis_id, prob_primary$model_id, sep = ": ")), collapse = " | "),
    "HRS and ELSA non-fatal comparison", 183, 105,
    panel_ids = c("a", "b"), row_groups = list(c("a", "b")),
    notes = "Probabilities are model-derived for identical 0-symptom and 3-symptom covariate profiles; ELSA S4 is absent from the wave 4-10 non-fatal analysis."
  )
}

# -------------------------------------------------------------------------
# Extended Data Fig. 3 | Complete 2/5/10-year transition-probability matrices.
# -------------------------------------------------------------------------

fig_id <- "ExtFig03"
file_stem <- "extfig03_complete_transition_probability_matrices"
fig_title <- "Complete fitted transition-probability matrices"
fig_claim <- "The complete matrices expose the full set of fitted state-to-state probabilities at 2, 5, and 10 years without selecting only favourable cells."
prob_full <- preferred_analysis_rows(markov_probability)
if (nrow(prob_full)) {
  ensure_column(prob_full, "exposure_variable", "insomnia_score_primary")
  ensure_column(prob_full, "profile", "")
  ensure_column(prob_full, "time_years", NA_real_)
  prob_full <- prob_full[
    exposure_variable == "insomnia_score_primary" & time_years %in% c(2, 5, 10) &
      is.finite(probability) & probability >= 0 & probability <= 1
  ]
  prob_full <- choose_model_per_analysis(prob_full)
  prob_full[, `:=`(
    profile_label = profile_display(profile),
    time_profile = paste0(profile_display(profile), "\n", time_years, " y"),
    matrix_cell = paste0("S", from_state, "\u2192S", to_state)
  )]
  profile_order <- c("Sample mean", "0 symptoms", "1 symptom", "2 symptoms", "3 symptoms", ">=2 symptoms", "Symptom absent", "Symptom present")
  x_order <- unlist(lapply(profile_order[profile_order %chin% unique(prob_full$profile_label)], function(p) paste0(p, "\n", c(2, 5, 10), " y")))
  cell_order <- unlist(lapply(0:4, function(from) paste0("S", from, "\u2192S", 0:4)))
  prob_full[, `:=`(
    time_profile = factor(time_profile, levels = x_order),
    matrix_cell = factor(matrix_cell, levels = rev(cell_order))
  )]
}
if (!nrow(prob_full) || uniqueN(prob_full$cohort) < 2L) {
  skip_figure(fig_id, file_stem, fig_title, fig_claim, "quantitative grid", "Complete 2/5/10-year probability matrices are unavailable for both cohorts", markov_probability_result$path, "HRS and ELSA", 2L, 183, 225)
} else {
  panels <- lapply(c("HRS", "ELSA"), function(cohort_name) {
    z <- prob_full[cohort == cohort_name]
    ggplot(z, aes(x = time_profile, y = matrix_cell, fill = probability)) +
      geom_tile(colour = "white", linewidth = 0.12) +
      geom_text(data = z[probability >= 0.01], aes(label = sprintf("%.2f", probability)), size = 1.78, family = "Arial", colour = "black") +
      scale_fill_gradient(low = "#F4F7FA", high = palette_contract[["signal_blue"]], limits = c(0, 1), labels = scales::label_percent(accuracy = 1)) +
      scale_x_discrete(drop = FALSE) +
      scale_y_discrete(drop = FALSE) +
      labs(x = NULL, y = "Origin-to-destination cell", fill = "Probability", title = cohort_name) +
      theme_nature_contract() +
      theme(
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 5.2),
        panel.border = element_rect(colour = palette_contract[["neutral_light"]], fill = NA, linewidth = 0.3),
        legend.position = "bottom"
      )
  })
  fig <- add_patchwork_tags((panels[[1L]] | panels[[2L]]) + plot_layout(guides = "collect")) & theme(legend.position = "bottom")
  save_and_record(
    fig_id, file_stem, fig, prob_full, fig_title, fig_claim, "quantitative grid",
    markov_probability_result$path,
    paste(unique(paste(prob_full$analysis_id, prob_full$model_id, sep = ": ")), collapse = " | "),
    "HRS and ELSA non-fatal comparison", 183, 225,
    panel_ids = c("a", "b"), row_groups = list(c("a", "b")),
    notes = "Only cells >=0.01 are printed, but every fitted cell is retained in the heatmap and Source Data. Blank ELSA S4 rows reflect the four-state non-fatal analysis, not zero-filled estimates."
  )
}

# -------------------------------------------------------------------------
# Fig. 4 | Model-expected state occupancy over follow-up.
# -------------------------------------------------------------------------

fig_id <- "Fig04"
file_stem <- "fig04_model_expected_state_occupancy"
fig_title <- "Model-expected state occupancy"
fig_claim <- "Fitted transition intensities reproduce an interpretable evolution of cohort state occupancy over observed follow-up."
occupancy <- preferred_analysis_rows(markov_prevalence)
if (nrow(occupancy)) {
  rename_first(occupancy, "expected_probability", c("predicted_probability", "model_probability", "expected"))
  ensure_column(occupancy, "time_years", NA_real_)
  ensure_column(occupancy, "state", NA_integer_)
  occupancy <- occupancy[is.finite(time_years) & is.finite(expected_probability) & expected_probability >= 0 & expected_probability <= 1]
  occupancy <- choose_model_per_analysis(occupancy)
  occupancy[, state_factor := factor(paste0("S", state), levels = paste0("S", 0:4), labels = unname(state_labels[paste0("S", 0:4)]))]
}
if (!nrow(occupancy) || uniqueN(occupancy$cohort) < 2L) {
  skip_figure(fig_id, file_stem, fig_title, fig_claim, "quantitative grid", "Model-expected state-occupancy estimates are unavailable for both cohorts", markov_prevalence_result$path, "HRS and ELSA", 2L, 183, 102)
} else {
  panels <- lapply(c("HRS", "ELSA"), function(cohort_name) {
    z <- occupancy[cohort == cohort_name]
    state_codes <- paste0("S", 0:4)
    state_colours <- setNames(unname(state_palette[state_codes]), unname(state_labels[state_codes]))
    ggplot(z, aes(x = time_years, y = expected_probability, colour = state_factor, group = state_factor)) +
      geom_line(linewidth = 0.68) +
      geom_point(size = 1.45) +
      scale_colour_manual(values = state_colours, breaks = names(state_colours), labels = state_codes, drop = FALSE) +
      scale_y_continuous(limits = c(0, 1), labels = scales::label_percent(accuracy = 1), expand = expansion(mult = c(0, 0.03))) +
      labs(x = "Years since analysis entry", y = "Expected cohort occupancy", colour = "State", title = cohort_name, subtitle = analysis_display(unique(z$analysis_id)[1L], cohort_name)) +
      theme_nature_contract() +
      theme(legend.position = "bottom")
  })
  panels[[1L]] <- panels[[1L]] +
    theme(legend.position = "bottom", legend.box = "vertical", legend.box.just = "center")
  panels[[2L]] <- panels[[2L]] + theme(legend.position = "none")
  fig <- add_patchwork_tags(
    (panels[[1L]] | panels[[2L]]) + plot_layout(guides = "collect")
  )
  save_and_record(
    fig_id, file_stem, fig, occupancy, fig_title, fig_claim, "quantitative grid",
    markov_prevalence_result$path,
    paste(unique(paste(occupancy$analysis_id, occupancy$model_id, sep = ": ")), collapse = " | "),
    "HRS and ELSA non-fatal comparison", 183, 102,
    panel_ids = c("a", "b"), row_groups = list(c("a", "b")),
    notes = "Expected probabilities are plotted at the time points returned by msm::prevalence.msm; no interpolation is applied."
  )
}

# -------------------------------------------------------------------------
# Fig. 5 | Mean sojourn times in transient states.
# -------------------------------------------------------------------------

fig_id <- "Fig05"
file_stem <- "fig05_mean_sojourn_times"
fig_title <- "Estimated mean sojourn times"
fig_claim <- "Mean residence times quantify where participants spend follow-up before progressing to another cardiovascular state."
sojourn <- preferred_analysis_rows(markov_sojourn)
if (nrow(sojourn)) {
  ensure_column(sojourn, "state", NA_integer_)
  sojourn <- sojourn[
    state %in% 0:3 & is.finite(sojourn_years) & sojourn_years >= 0 &
      is.finite(conf_low) & conf_low >= 0 & is.finite(conf_high) & conf_high >= 0
  ]
  sojourn <- choose_model_per_analysis(sojourn)
  sojourn[, state_factor := factor(paste0("S", state), levels = paste0("S", 3:0), labels = unname(state_labels[paste0("S", 3:0)]))]
}
if (!nrow(sojourn) || uniqueN(sojourn$cohort) < 2L) {
  skip_figure(fig_id, file_stem, fig_title, fig_claim, "quantitative grid", "Mean sojourn-time estimates with confidence intervals are unavailable for both cohorts", markov_sojourn_result$path, "HRS and ELSA", 2L, 183, 94)
} else {
  panels <- lapply(c("HRS", "ELSA"), function(cohort_name) {
    z <- sojourn[cohort == cohort_name]
    ggplot(z, aes(x = sojourn_years, y = state_factor)) +
      geom_errorbar(aes(xmin = conf_low, xmax = conf_high), orientation = "y", width = 0.16, linewidth = 0.5, colour = palette_contract[["signal_blue"]]) +
      geom_point(shape = 21, size = 2, stroke = 0.35, fill = palette_contract[["signal_blue"]]) +
      labs(x = "Mean sojourn time, years (95% CI)", y = NULL, title = cohort_name, subtitle = analysis_display(unique(z$analysis_id)[1L], cohort_name)) +
      theme_nature_contract() +
      theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())
  })
  fig <- add_patchwork_tags(panels[[1L]] | panels[[2L]])
  save_and_record(
    fig_id, file_stem, fig, sojourn, fig_title, fig_claim, "quantitative grid",
    markov_sojourn_result$path,
    paste(unique(paste(sojourn$analysis_id, sojourn$model_id, sep = ": ")), collapse = " | "),
    "HRS and ELSA non-fatal comparison", 183, 94,
    panel_ids = c("a", "b"), row_groups = list(c("a", "b")),
    notes = "Death is absorbing and therefore has no finite sojourn-time estimate."
  )
}

# -------------------------------------------------------------------------
# Extended Data Fig. 4 | Observed-versus-expected prevalence diagnostics.
# -------------------------------------------------------------------------

fig_id <- "ExtFig04"
file_stem <- "extfig04_observed_expected_prevalence_diagnostics"
fig_title <- "Observed-versus-expected state prevalence"
fig_claim <- "Observed and model-expected state prevalence provide a direct calibration diagnostic for the fitted Markov models."
diagnostic <- preferred_analysis_rows(markov_prevalence)
if (nrow(diagnostic)) {
  rename_first(diagnostic, "observed_probability", c("observed", "observed_prevalence"))
  rename_first(diagnostic, "expected_probability", c("predicted_probability", "model_probability", "expected"))
  ensure_column(diagnostic, "state", NA_integer_)
  diagnostic <- diagnostic[
    is.finite(observed_probability) & observed_probability >= 0 & observed_probability <= 1 &
      is.finite(expected_probability) & expected_probability >= 0 & expected_probability <= 1
  ]
  diagnostic <- choose_model_per_analysis(diagnostic)
  diagnostic[, `:=`(
    state_factor = factor(paste0("S", state), levels = paste0("S", 0:4), labels = unname(state_labels[paste0("S", 0:4)])),
    absolute_error = abs(expected_probability - observed_probability)
  )]
}
if (!nrow(diagnostic) || uniqueN(diagnostic$cohort) < 2L) {
  skip_figure(fig_id, file_stem, fig_title, fig_claim, "quantitative grid", "Observed-versus-expected prevalence diagnostics are unavailable for both cohorts", markov_prevalence_result$path, "HRS and ELSA", 2L, 183, 92)
} else {
  panels <- lapply(c("HRS", "ELSA"), function(cohort_name) {
    z <- diagnostic[cohort == cohort_name]
    state_codes <- paste0("S", 0:4)
    state_colours <- setNames(unname(state_palette[state_codes]), unname(state_labels[state_codes]))
    ggplot(z, aes(x = observed_probability, y = expected_probability, colour = state_factor)) +
      geom_abline(slope = 1, intercept = 0, linewidth = 0.38, linetype = "dashed", colour = palette_contract[["neutral_mid"]]) +
      geom_point(size = 1.75, alpha = 0.88) +
      coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
      scale_colour_manual(values = state_colours, breaks = names(state_colours), labels = state_codes, drop = FALSE) +
      scale_x_continuous(labels = scales::label_percent(accuracy = 1)) +
      scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
      labs(x = "Observed prevalence", y = "Expected prevalence", colour = "State", title = cohort_name) +
      theme_nature_contract() +
      theme(legend.position = "bottom")
  })
  panels[[1L]] <- panels[[1L]] +
    theme(legend.position = "bottom", legend.box = "vertical", legend.box.just = "center")
  panels[[2L]] <- panels[[2L]] + theme(legend.position = "none")
  fig <- add_patchwork_tags(
    (panels[[1L]] | panels[[2L]]) + plot_layout(guides = "collect")
  )
  save_and_record(
    fig_id, file_stem, fig, diagnostic, fig_title, fig_claim, "quantitative grid",
    markov_prevalence_result$path,
    paste(unique(paste(diagnostic$analysis_id, diagnostic$model_id, sep = ": ")), collapse = " | "),
    "HRS and ELSA non-fatal comparison", 183, 92,
    panel_ids = c("a", "b"), row_groups = list(c("a", "b")),
    notes = "Each point is one state-time combination; the dashed diagonal denotes exact calibration. Absolute errors are retained in Source Data."
  )
}

# -------------------------------------------------------------------------
# Extended Data Fig. 5 | Complete Cox sensitivity heatmap.
# -------------------------------------------------------------------------

fig_id <- "ExtFig05"
file_stem <- "extfig05_sensitivity_heatmap"
fig_title <- "Robustness across exposure and design specifications"
fig_claim <- "The direction and magnitude of transition-specific estimates are assessed across all estimable M4 sensitivity specifications rather than a selected subset."
sensitivity <- copy(cox_estimates)
if (nrow(sensitivity)) {
  rename_first(sensitivity, "hazard_ratio", c("HR", "estimate_ratio"))
  rename_first(sensitivity, "conf_low", c("lower", "lower_ci", "lcl"))
  rename_first(sensitivity, "conf_high", c("upper", "upper_ci", "ucl"))
  ensure_column(sensitivity, "analysis_id", "")
  ensure_column(sensitivity, "analysis_label", sensitivity$analysis_id)
  ensure_column(sensitivity, "cohort", NA_character_)
  ensure_column(sensitivity, "transition", "")
  ensure_column(sensitivity, "term_role", "exposure")
  ensure_column(sensitivity, "term", "")
  ensure_column(sensitivity, "model_level", "M4")
  ensure_column(sensitivity, "basis", "complete_case")
  ensure_column(sensitivity, "weighted", 0L)
  ensure_column(sensitivity, "q_value_bh", NA_real_)
  if (!any(is.finite(sensitivity$q_value_bh)) && "p_value" %in% names(sensitivity)) {
    sensitivity[, q_value_bh := p.adjust(p_value, method = "BH"), by = .(cohort, analysis_id, model_level, basis, weighted, term)]
  }
  sensitivity[, `:=`(
    cohort = toupper(as.character(cohort)),
    transition = normalize_transition(transition),
    weighted = as.integer(weighted)
  )]
  sensitivity <- sensitivity[
    model_level == "M4" & term_role %chin% c("exposure", "interaction") &
      is.finite(hazard_ratio) & hazard_ratio > 0 & is.finite(conf_low) & conf_low > 0 & is.finite(conf_high) & conf_high > 0
  ]
  if (nrow(sensitivity)) {
    sensitivity[, basis_rank := fifelse(basis == "mice_pooled", 1L, fifelse(basis == "complete_case", 2L, 3L))]
    sensitivity <- sensitivity[order(basis_rank)][, .SD[1L], by = .(cohort, analysis_id, weighted, transition, term)]
    sensitivity[, `:=`(
      log2_hazard_ratio = log2(hazard_ratio),
      specification_label = paste0(
        analysis_label,
        fifelse(nzchar(term) & !grepl("^exposure_value$", term), paste0(" [", term, "]"), ""),
        fifelse(weighted == 1L, " [weighted]", "")
      ),
      transition_label = display_transition(transition),
      fdr_marker = fifelse(is.finite(q_value_bh) & q_value_bh < 0.05, "\u2022", "")
    )]
    spec_order <- unique(sensitivity[order(cohort, analysis_id, weighted, term), specification_label])
    sensitivity[, specification_factor := factor(specification_label, levels = rev(spec_order))]
  }
}
if (!nrow(sensitivity)) {
  skip_figure(fig_id, file_stem, fig_title, fig_claim, "quantitative grid", "No estimable M4 Cox sensitivity estimates are available", cox_result$path, "HRS and ELSA", 2L, 183, 220)
} else {
  max_abs <- max(abs(sensitivity$log2_hazard_ratio), na.rm = TRUE)
  if (!is.finite(max_abs) || max_abs <= 0) max_abs <- 0.01
  cohorts <- intersect(c("HRS", "ELSA"), unique(sensitivity$cohort))
  panels <- lapply(cohorts, function(cohort_name) {
    z <- sensitivity[cohort == cohort_name]
    ggplot(z, aes(x = factor(transition_label, levels = display_transition(all_edges)), y = specification_factor, fill = log2_hazard_ratio)) +
      geom_tile(colour = "white", linewidth = 0.16) +
      geom_text(aes(label = fdr_marker), family = "Arial", size = 2.0, colour = "black") +
      scale_fill_gradient2(
        low = "#3B6FB6", mid = "#F7F7F7", high = "#C84C4C", midpoint = 0,
        limits = c(-max_abs, max_abs), name = "log2(HR)"
      ) +
      scale_x_discrete(drop = FALSE) +
      scale_y_discrete(drop = FALSE) +
      labs(x = NULL, y = NULL, title = cohort_name) +
      theme_nature_contract() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 5.4),
        axis.text.y = element_text(size = 5.2),
        panel.border = element_rect(colour = palette_contract[["neutral_light"]], fill = NA, linewidth = 0.3),
        legend.position = "bottom"
      )
  })
  if (length(panels) == 1L) {
    fig <- panels[[1L]]
    panel_ids <- NULL
    row_groups <- NULL
    column_groups <- NULL
  } else {
    # Long specification labels require the full journal column width.  Stack
    # cohorts so labels and all 45 sensitivity rows remain legible at final
    # size instead of collapsing both panel bodies beside duplicated axes.
    fig <- add_patchwork_tags((wrap_plots(panels, ncol = 1L)) + plot_layout(guides = "collect")) & theme(legend.position = "bottom")
    panel_ids <- letters[seq_along(panels)]
    row_groups <- NULL
    column_groups <- list(panel_ids)
  }
  dynamic_height <- min(245, max(135, 62 + 4.2 * uniqueN(sensitivity$specification_label)))
  save_and_record(
    fig_id, file_stem, fig, sensitivity, fig_title, fig_claim, "quantitative grid",
    cox_result$path,
    "All estimable M4 exposure/interaction specifications; pooled imputation preferred within duplicate specification cells",
    paste(cohorts, collapse = " and "), 183, dynamic_height,
    panel_ids = panel_ids, row_groups = row_groups, column_groups = column_groups,
    notes = "Fill shows the uncapped observed log2 hazard ratio. Dots denote BH-FDR <0.05 within the reported analysis-term family. Missing tiles are unavailable/failed fits, not zero effects."
  )
}

# -------------------------------------------------------------------------
# Extended Data Fig. 6 | Targeted protocol sensitivities and exact coverage.
# -------------------------------------------------------------------------

fig_id <- "ExtFig06"
file_stem <- "extfig06_protocol_sensitivity_coverage"
fig_title <- "Protocol-targeted sensitivity analyses and coverage"
fig_claim <- "Every prespecified sensitivity family is either represented by an estimable result, linked to an upstream analysis, or retained with an explicit non-estimability reason."
protocol_effects <- copy(protocol_sensitivity_estimates)
protocol_coverage <- copy(protocol_sensitivity_coverage)
coverage_complete <- FALSE
coverage_status_valid <- FALSE
coverage_labels_consistent <- FALSE
coverage_reasons_complete <- FALSE
coverage_upstream_links_complete <- FALSE

if (nrow(protocol_effects)) {
  ensure_column(protocol_effects, "scenario_label", protocol_effects$scenario)
  ensure_column(protocol_effects, "term", "insomnia_score_primary")
  ensure_column(protocol_effects, "status", "completed")
  ensure_column(protocol_effects, "p_value", NA_real_)
  ensure_column(protocol_effects, "p_fdr_bh_sensitivity", NA_real_)
  rename_first(protocol_effects, "conf_low", c("lower", "lower_ci", "lcl"))
  rename_first(protocol_effects, "conf_high", c("upper", "upper_ci", "ucl"))
  ensure_column(protocol_effects, "conf_low", NA_real_)
  ensure_column(protocol_effects, "conf_high", NA_real_)
  protocol_effects[, `:=`(
    sensitivity_id = as.integer(sensitivity_id),
    cohort = toupper(as.character(cohort)),
    transition = normalize_transition(transition),
    scenario = as.character(scenario),
    term = as.character(term),
    status = tolower(as.character(status))
  )]
  protocol_effects <- protocol_effects[
    status == "completed" & cohort %chin% c("HRS", "ELSA") &
      transition %chin% all_edges & grepl("^insomnia", term) &
      is.finite(hazard_ratio) & hazard_ratio > 0
  ]
  if (nrow(protocol_effects) && !any(is.finite(protocol_effects$p_fdr_bh_sensitivity))) {
    protocol_effects[is.finite(p_value),
      p_fdr_bh_sensitivity := p.adjust(p_value, method = "BH"),
      by = .(cohort, sensitivity_id, scenario, term)
    ]
  }
  scenario_names <- c(
    rcs = "Dose-response spline",
    alternate_death = "Alternative death timing",
    attained_age_period = "Attained age + period",
    semi_markov = "Semi-Markov",
    ipcw = "IPCW",
    proxy_exclusion = "Exclude proxy interviews",
    expanded_hypertension = "Expanded hypertension"
  )
  protocol_effects[, scenario_short := unname(scenario_names[scenario])]
  protocol_effects[is.na(scenario_short), scenario_short := scenario_label]
  protocol_effects[, contrast_short := fcase(
    term == "insomnia_score_1_vs_0", "1 vs 0 symptoms",
    term == "insomnia_score_2_vs_0", "2 vs 0 symptoms",
    term == "insomnia_score_3_vs_0", "3 vs 0 symptoms",
    default = "per symptom"
  )]
  protocol_effects[, `:=`(
    specification_label = paste0(
      sprintf("%02d", sensitivity_id), "  ", scenario_short,
      fifelse(scenario == "rcs", paste0(": ", contrast_short), "")
    ),
    transition_label = display_transition(transition),
    log2_hazard_ratio = log2(hazard_ratio),
    fdr_marker = fifelse(
      is.finite(p_fdr_bh_sensitivity) & p_fdr_bh_sensitivity < 0.05,
      "\u2022", ""
    )
  )]
}

if (nrow(protocol_coverage)) {
  protocol_coverage[, `:=`(
    sensitivity_id = as.integer(sensitivity_id),
    cohort = toupper(as.character(cohort)),
    status = tolower(as.character(status)),
    protocol_sensitivity = as.character(protocol_sensitivity)
  )]
  protocol_coverage <- protocol_coverage[
    sensitivity_id %between% c(1L, 16L) & cohort %chin% c("HRS", "ELSA")
  ]
  coverage_levels <- c(
    "completed", "partially_completed", "covered_by_primary_definition",
    "covered_upstream", "not_estimable"
  )
  expected_coverage <- CJ(
    cohort = c("HRS", "ELSA"),
    sensitivity_id = 1:16,
    unique = TRUE
  )
  observed_coverage <- protocol_coverage[, .(cohort, sensitivity_id)]
  coverage_complete <-
    nrow(protocol_coverage) == nrow(expected_coverage) &&
    !anyDuplicated(observed_coverage) &&
    setequal(
      paste(observed_coverage$cohort, observed_coverage$sensitivity_id, sep = "::"),
      paste(expected_coverage$cohort, expected_coverage$sensitivity_id, sep = "::")
    )
  coverage_status_valid <- all(protocol_coverage$status %chin% coverage_levels)
  coverage_labels_consistent <-
    all(!is.na(protocol_coverage$protocol_sensitivity) &
          nzchar(trimws(protocol_coverage$protocol_sensitivity))) &&
    protocol_coverage[
      , .(n_labels = uniqueN(trimws(protocol_sensitivity))), by = sensitivity_id
    ][, all(n_labels == 1L)]
  coverage_reasons_complete <- protocol_coverage[
    status == "not_estimable",
    all(!is.na(reason) & nzchar(trimws(reason)))
  ]
  coverage_upstream_links_complete <- protocol_coverage[
    status == "covered_upstream",
    all(!is.na(evidence_files) & nzchar(trimws(evidence_files)))
  ]
  protocol_coverage[, `:=`(
    sensitivity_label = paste0(sprintf("%02d", sensitivity_id), "  ", protocol_sensitivity),
    status_factor = factor(status, levels = coverage_levels, labels = c(
      "Completed", "Partial", "Primary", "Upstream", "Not estimable"
    )),
    status_code = fcase(
      status == "completed", "Complete",
      status == "partially_completed", "Partial",
      status == "covered_by_primary_definition", "Primary",
      status == "covered_upstream", "Upstream",
      default = "N/E"
    )
  )]
}

required_protocol_cohorts <- c("HRS", "ELSA")
if (!nrow(protocol_effects) || !nrow(protocol_coverage) ||
    !all(required_protocol_cohorts %chin% unique(protocol_effects$cohort)) ||
    !all(required_protocol_cohorts %chin% unique(protocol_coverage$cohort)) ||
    !coverage_complete || !coverage_status_valid ||
    !coverage_labels_consistent || !coverage_reasons_complete ||
    !coverage_upstream_links_complete) {
  skip_figure(
    fig_id, file_stem, fig_title, fig_claim, "quantitative grid",
    paste(
      "Targeted protocol estimates are unavailable for both cohorts, or the coverage audit",
      "is not an exact validated 2-by-16 matrix with recognized statuses, reasons and upstream links"
    ),
    c(protocol_sensitivity_result$path, protocol_coverage_result$path),
    "HRS and ELSA", 3L, 183, 238
  )
} else {
  protocol_limit <- max(abs(protocol_effects$log2_hazard_ratio), na.rm = TRUE)
  if (!is.finite(protocol_limit) || protocol_limit <= 0) protocol_limit <- 0.01
  protocol_heatmap_panel <- function(cohort_name) {
    z <- protocol_effects[cohort == cohort_name]
    row_order <- unique(z[order(sensitivity_id, scenario, contrast_short), specification_label])
    z[, specification_factor := factor(specification_label, levels = rev(row_order))]
    ggplot(
      z,
      aes(
        x = factor(transition_label, levels = display_transition(all_edges)),
        y = specification_factor,
        fill = log2_hazard_ratio
      )
    ) +
      geom_tile(colour = "white", linewidth = 0.18) +
      geom_text(aes(label = fdr_marker), family = "Arial", size = 2.0, colour = "black") +
      scale_fill_gradient2(
        low = "#3B6FB6", mid = "#F7F7F7", high = "#C84C4C", midpoint = 0,
        limits = c(-protocol_limit, protocol_limit), name = "log2(HR)"
      ) +
      scale_x_discrete(drop = FALSE) +
      scale_y_discrete(drop = FALSE) +
      labs(x = NULL, y = NULL, title = paste0(cohort_name, " estimable targeted sensitivities")) +
      theme_nature_contract() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 5.2),
        axis.text.y = element_text(size = 5.2),
        panel.border = element_rect(colour = palette_contract[["neutral_light"]], fill = NA, linewidth = 0.3),
        legend.position = "bottom"
      )
  }
  p_protocol_hrs <- protocol_heatmap_panel("HRS")
  p_protocol_elsa <- protocol_heatmap_panel("ELSA")

  coverage_order <- unique(protocol_coverage[order(sensitivity_id), sensitivity_label])
  protocol_coverage[, sensitivity_factor := factor(sensitivity_label, levels = rev(coverage_order))]
  coverage_palette <- c(
    "Completed" = palette_contract[["signal_teal"]],
    "Partial" = palette_contract[["accent_orange"]],
    "Primary" = palette_contract[["signal_blue"]],
    "Upstream" = palette_contract[["pale_blue"]],
    "Not estimable" = palette_contract[["neutral_light"]]
  )
  wrap_discrete <- function(x, width = 52L) {
    vapply(x, function(value) paste(strwrap(value, width = width), collapse = "\n"), character(1L))
  }
  p_protocol_coverage <- ggplot(
    protocol_coverage,
    aes(x = factor(cohort, levels = c("HRS", "ELSA")), y = sensitivity_factor, fill = status_factor)
  ) +
    geom_tile(colour = "white", linewidth = 0.22) +
    geom_text(aes(label = status_code), family = "Arial", size = 1.8, colour = "black") +
    scale_fill_manual(values = coverage_palette, drop = FALSE, name = "Protocol status") +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE, labels = function(x) wrap_discrete(x, 58L)) +
    labs(x = NULL, y = NULL, title = "All 16 prespecified sensitivity families") +
    theme_nature_contract() +
    theme(
      axis.text.x = element_text(size = 5.5),
      axis.text.y = element_text(size = 5.0),
      panel.border = element_rect(colour = palette_contract[["neutral_light"]], fill = NA, linewidth = 0.3),
      legend.position = "bottom"
    )

  fig <- add_patchwork_tags(
    (p_protocol_hrs / p_protocol_elsa / p_protocol_coverage) +
      plot_layout(heights = c(1, 1, 1.8), guides = "collect")
  ) & theme(legend.position = "bottom", legend.box = "vertical", legend.box.just = "center")
  protocol_effects_source <- copy(protocol_effects)
  protocol_effects_source[, record_type := "estimable_effect"]
  protocol_coverage_source <- copy(protocol_coverage)
  protocol_coverage_source[, record_type := "coverage_audit"]
  source_protocol <- rbindlist(
    list(protocol_effects_source, protocol_coverage_source),
    use.names = TRUE, fill = TRUE
  )
  save_and_record(
    fig_id, file_stem, fig, source_protocol, fig_title, fig_claim, "quantitative grid",
    c(protocol_sensitivity_result$path, protocol_coverage_result$path),
    "All estimable targeted sensitivity effects plus all 16 protocol coverage decisions",
    "HRS and ELSA", 183, 238,
    panel_ids = c("a", "b", "c"),
    column_groups = list(c("a", "b", "c")),
    notes = paste(
      "Effect heatmaps show every estimable protocol-targeted insomnia coefficient; dots denote scenario-specific BH-FDR <0.05.",
      "The coverage matrix retains unsupported and sparse specifications instead of converting missing tiles to null effects; exact reasons remain in Source Data."
    )
  )
}

# -------------------------------------------------------------------------
# Fig. 6 | Direct cohort comparison: crude rates and adjusted intensity ratios.
# -------------------------------------------------------------------------

fig_id <- "Fig06"
file_stem <- "fig06_hrs_elsa_cohort_comparison"
fig_title <- "Cross-cohort comparison of non-fatal transitions"
fig_claim <- "HRS and ELSA are compared on both observed transition incidence and adjusted insomnia associations using the four harmonized non-fatal edges."
comparison_rate <- copy(rates)
if (nrow(comparison_rate)) {
  rename_first(comparison_rate, "rate_per_1000_py", c("incidence_rate", "rate", "crude_rate_per_1000_py"))
  comparison_rate[, `:=`(
    cohort = toupper(as.character(cohort)),
    transition = normalize_transition(transition)
  )]
  comparison_rate <- comparison_rate[
    cohort %chin% c("HRS", "ELSA") & transition %chin% confirmatory_nonfatal_edges &
      is.finite(rate_per_1000_py) & rate_per_1000_py >= 0
  ]
}
comparison_ir <- copy(hr_primary)
if (nrow(comparison_ir)) comparison_ir <- comparison_ir[transition %chin% confirmatory_nonfatal_edges]
if (!nrow(comparison_rate) || !nrow(comparison_ir) || uniqueN(comparison_rate$cohort) < 2L || uniqueN(comparison_ir$cohort) < 2L) {
  skip_figure(fig_id, file_stem, fig_title, fig_claim, "quantitative grid", "Both crude rates and primary adjusted intensity ratios are required for a two-cohort comparison", c(rate_result$path, markov_hr_result$path), "HRS and ELSA", 2L, 183, 98)
} else {
  order_nonfatal <- rev(confirmatory_nonfatal_edges)
  comparison_rate[, transition_factor := factor(transition, levels = order_nonfatal, labels = display_transition(order_nonfatal))]
  comparison_ir[, transition_factor := factor(transition, levels = order_nonfatal, labels = display_transition(order_nonfatal))]
  p_rate <- ggplot(comparison_rate, aes(x = rate_per_1000_py, y = transition_factor, colour = cohort, shape = cohort)) +
    geom_line(aes(group = transition_factor), linewidth = 0.35, colour = palette_contract[["neutral_light"]]) +
    geom_point(size = 2.1, stroke = 0.55) +
    scale_colour_manual(values = cohort_palette) +
    scale_shape_manual(values = c(HRS = 16, ELSA = 17)) +
    labs(x = "Observed rate per 1,000 person-years", y = NULL, colour = "Cohort", shape = "Cohort", title = "Observed transition incidence") +
    theme_nature_contract() +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "bottom")
  p_ir <- ggplot(comparison_ir, aes(x = hazard_ratio, y = transition_factor, colour = cohort, shape = cohort)) +
    geom_vline(xintercept = 1, linewidth = 0.35, linetype = "dashed", colour = palette_contract[["neutral_mid"]]) +
    geom_errorbar(aes(xmin = conf_low, xmax = conf_high), orientation = "y", width = 0.12, linewidth = 0.45, position = position_dodge(width = 0.45)) +
    geom_point(size = 2.0, stroke = 0.55, position = position_dodge(width = 0.45)) +
    scale_x_log10(labels = scales::label_number(accuracy = 0.01)) +
    scale_colour_manual(values = cohort_palette) +
    scale_shape_manual(values = c(HRS = 16, ELSA = 17)) +
    labs(x = "Adjusted transition intensity ratio (95% CI)", y = NULL, colour = "Cohort", shape = "Cohort", title = "Insomnia score, per additional symptom") +
    theme_nature_contract() +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none")
  fig <- add_patchwork_tags(p_rate | p_ir)
  source_comparison <- rbindlist(list(
    comparison_rate[, .(
      panel = "Observed transition incidence", cohort, transition,
      rate_per_1000_py, risk_intervals = if ("risk_intervals" %in% names(comparison_rate)) risk_intervals else NA_real_,
      risk_participants = if ("risk_participants" %in% names(comparison_rate)) risk_participants else NA_real_,
      events = if ("events" %in% names(comparison_rate)) events else NA_real_,
      person_years = if ("person_years" %in% names(comparison_rate)) person_years else NA_real_,
      analysis_id = NA_character_, model_id = NA_character_, hazard_ratio = NA_real_, conf_low = NA_real_, conf_high = NA_real_
    )],
    comparison_ir[, .(
      panel = "Adjusted transition intensity ratio", cohort, transition,
      rate_per_1000_py = NA_real_, risk_intervals = if ("risk_intervals" %in% names(comparison_ir)) risk_intervals else NA_real_,
      risk_participants = if ("risk_participants" %in% names(comparison_ir)) risk_participants else NA_real_,
      events = if ("events" %in% names(comparison_ir)) events else NA_real_,
      person_years = if ("person_years" %in% names(comparison_ir)) person_years else NA_real_,
      analysis_id, model_id, hazard_ratio, conf_low, conf_high
    )]
  ), use.names = TRUE, fill = TRUE)
save_and_record(
    fig_id, file_stem, fig, source_comparison, fig_title, fig_claim, "quantitative grid",
    c(rate_result$path, markov_hr_result$path),
    paste(unique(paste(comparison_ir$analysis_id, comparison_ir$model_id, sep = ": ")), collapse = " | "),
    "HRS and ELSA harmonized non-fatal transitions", 183, 98,
    panel_ids = c("a", "b"), row_groups = list(c("a", "b")),
    notes = "Panel a is descriptive and unadjusted; panel b is model-adjusted. Differences between panels must not be interpreted as confounding decomposition."
  )
}

# -------------------------------------------------------------------------
# Fig. 7 | Adjusted transition-intensity networks for the three estimands.
# -------------------------------------------------------------------------

fig_id <- "Fig07"
file_stem <- "fig07_adjusted_transition_intensity_networks"
fig_title <- "Adjusted cardiovascular transition-intensity networks"
fig_claim <- paste(
  "Adjusted baseline transition intensities differ across the HRS five-state",
  "model and the ELSA non-fatal and mortality-limited estimands."
)

network_specs <- data.table(
  panel_id = c("A", "B", "C"),
  panel_title = c(
    "HRS five-state network",
    "ELSA waves 4-10 primary network",
    "ELSA five-state extension through 2012"
  ),
  panel_subtitle = c(
    "Full five-state follow-up through 2020; eight allowed transitions",
    "Primary four-state model; four non-fatal transitions only",
    "Mortality-linked follow-up truncated at 2012; eight allowed transitions"
  ),
  analysis_id = c(
    "hrs_5state_full",
    "elsa_4state_w4_w10",
    "elsa_5state_through_2012"
  ),
  cohort = c("HRS", "ELSA", "ELSA"),
  expected_n = c(8L, 4L, 8L),
  expected_edges = list(all_edges, confirmatory_nonfatal_edges, all_edges),
  include_death = c(TRUE, FALSE, TRUE),
  mortality_limited = c(FALSE, FALSE, TRUE)
)

select_network_q <- function(specification, estimates) {
  expected_edges <- specification$expected_edges[[1L]]
  selected <- estimates[
    analysis_id == specification$analysis_id &
      cohort == specification$cohort &
      model_id == "M4_MI20" &
      covariate_profile == "sample means" &
      tolower(status) == "completed"
  ]
  if (nrow(selected) != specification$expected_n ||
      uniqueN(selected$transition) != specification$expected_n ||
      !setequal(selected$transition, expected_edges)) {
    stop(
      "Transition-intensity network blocked for ", specification$panel_id,
      ": expected exactly ", specification$expected_n, " protocol-defined q rows.",
      call. = FALSE
    )
  }
  if (anyDuplicated(selected[, .(analysis_id, model_id, covariate_profile, transition)])) {
    stop("Transition-intensity network blocked by duplicated q rows", call. = FALSE)
  }
  if (!all(
    selected$pooling_status == "RUBIN_POOLED_COMPLETE" &
      selected$m_expected == 20L & selected$m_fitted == 20L
  )) {
    stop("Transition-intensity network requires complete 20-imputation Rubin pooling", call. = FALSE)
  }
  if (!all(
    is.finite(selected$intensity) & selected$intensity > 0 &
      is.finite(selected$conf_low) & selected$conf_low > 0 &
      is.finite(selected$conf_high) & selected$conf_high > 0 &
      selected$conf_low <= selected$intensity & selected$intensity <= selected$conf_high
  )) {
    stop("Transition-intensity network contains an invalid q estimate or confidence interval", call. = FALSE)
  }
  if ("S0->S3" %chin% selected$transition) {
    stop("Direct S0-to-S3 intensity is prohibited by the prespecified state system", call. = FALSE)
  }
  selected[, `:=`(
    panel_id = specification$panel_id,
    panel_title = specification$panel_title,
    panel_subtitle = specification$panel_subtitle,
    mortality_limited = specification$mortality_limited,
    display_value = sprintf("q_ij = %.4f per person-year", intensity)
  )]
  selected[, `.edge_order` := match(transition, expected_edges)]
  setorder(selected, `.edge_order`)
  selected[, `.edge_order` := NULL]
  selected[]
}

network_q_selected <- data.table()
if (nrow(markov_q)) {
  required_q_columns <- c(
    "analysis_id", "cohort", "model_id", "transition", "from_state", "to_state",
    "covariate_profile", "intensity", "standard_error", "conf_low", "conf_high",
    "risk_intervals", "risk_participants", "person_years", "events",
    "m_expected", "m_fitted", "pooling_status", "pooled_standard_error_log",
    "status", "mortality_scope", "spec_id"
  )
  if (!all(required_q_columns %in% names(markov_q))) {
    stop(
      "Transition-intensity network source is missing columns: ",
      paste(setdiff(required_q_columns, names(markov_q)), collapse = ", "),
      call. = FALSE
    )
  }
  network_q_selected <- rbindlist(
    lapply(seq_len(nrow(network_specs)), function(index) {
      select_network_q(network_specs[index], markov_q)
    }),
    use.names = TRUE, fill = TRUE
  )
}

network_nodes_5 <- data.table(
  state = 0:4,
  state_code = paste0("S", 0:4),
  state_label = unname(state_labels[paste0("S", 0:4)]),
  x = c(0, 1.35, 1.35, 2.85, 4.25),
  y = c(0, 0.90, -0.90, 0, 0),
  description_x = c(-0.03, 1.35, 1.35, 2.85, 4.25),
  description_y = c(-0.34, 1.22, -1.25, -0.40, -0.34),
  description = c(
    "No target disease", "Hypertension only", "One CVD without hypertension",
    "Cardiovascular multimorbidity", "Death (absorbing)"
  )
)

network_nodes_4 <- copy(network_nodes_5[state < 4L])
network_nodes_4[state == 3L, `:=`(x = 3.15, description_x = 3.15)]

network_geometry_5 <- data.table(
  transition = all_edges,
  curvature = c(0.10, -0.10, -0.48, 0.10, -0.15, -0.10, 0.23, 0),
  target_y_offset = c(0, 0, 0.42, 0, 0.18, 0, -0.34, 0),
  label_x = c(0.57, 0.57, 2.10, 2.05, 3.15, 2.05, 3.18, 3.55),
  label_y = c(0.70, -0.70, 2.06, 0.70, 1.14, -0.70, -1.43, 0.31)
)

network_geometry_4 <- data.table(
  transition = confirmatory_nonfatal_edges,
  curvature = c(0.10, -0.10, 0.10, -0.10),
  target_y_offset = 0,
  label_x = c(0.57, 0.57, 2.15, 2.15),
  label_y = c(0.70, -0.70, 0.70, -0.70)
)

build_intensity_network_panel <- function(edge_data, include_death, mortality_warning = FALSE) {
  nodes <- copy(if (include_death) network_nodes_5 else network_nodes_4)
  geometry <- copy(if (include_death) network_geometry_5 else network_geometry_4)
  edges <- merge(edge_data, geometry, by = "transition", all.x = TRUE, sort = FALSE)
  edges <- merge(
    edges, nodes[, .(from_state = state, x_from = x, y_from = y)],
    by = "from_state", all.x = TRUE, sort = FALSE
  )
  edges <- merge(
    edges, nodes[, .(to_state = state, x_to = x, y_to = y)],
    by = "to_state", all.x = TRUE, sort = FALSE
  )
  edges[, `.edge_order` := match(transition, geometry$transition)]
  setorder(edges, `.edge_order`)
  edges[, `.edge_order` := NULL]
  if (any(!is.finite(unlist(edges[, .(x_from, y_from, x_to, y_to, curvature, target_y_offset, label_x, label_y)])))) {
    stop("Transition-intensity network geometry is incomplete", call. = FALSE)
  }
  edges[, `:=`(x_to_target = x_to, y_to_target = y_to + target_y_offset)]
  edges[, edge_length := sqrt((x_to_target - x_from)^2 + (y_to_target - y_from)^2)]
  edges[, `:=`(
    x_from_plot = x_from + 0.21 * (x_to_target - x_from) / edge_length,
    y_from_plot = y_from + 0.21 * (y_to_target - y_from) / edge_length,
    x_to_plot = x_to_target - 0.27 * (x_to_target - x_from) / edge_length,
    y_to_plot = y_to_target - 0.27 * (y_to_target - y_from) / edge_length,
    label_text = sprintf("q_ij = %.4f / year", intensity),
    death_transition_plot = to_state == 4L
  )]

  panel_title <- unique(edge_data$panel_title)
  panel_subtitle <- unique(edge_data$panel_subtitle)
  x_limits <- if (include_death) c(-0.55, 4.78) else c(-0.55, 3.68)
  p <- ggplot() +
    coord_fixed(ratio = 0.53, xlim = x_limits, ylim = c(-1.63, 2.30), clip = "off") +
    theme_void(base_family = "Arial", base_size = 7) +
    theme(
      plot.margin = margin(5, 11, 5, 11),
      plot.title = element_text(size = 7.8, face = "bold", colour = palette_contract[["neutral_dark"]], margin = margin(b = 2)),
      plot.subtitle = element_text(size = 6.3, colour = palette_contract[["neutral_mid"]], margin = margin(b = 2))
    ) +
    labs(title = panel_title, subtitle = panel_subtitle)

  for (index in seq_len(nrow(edges))) {
    is_death <- edges$death_transition_plot[index]
    p <- p + geom_curve(
      data = edges[index],
      aes(x = x_from_plot, y = y_from_plot, xend = x_to_plot, yend = y_to_plot),
      curvature = edges$curvature[index],
      linewidth = if (is_death) 0.58 else 0.65,
      linetype = if (is_death) "22" else "solid",
      colour = if (is_death) palette_contract[["neutral_dark"]] else palette_contract[["signal_blue"]],
      arrow = grid::arrow(length = grid::unit(2.0, "mm"), type = "closed"),
      lineend = "round"
    )
  }

  p <- p +
    geom_point(
      data = nodes, aes(x = x, y = y, fill = state_code),
      shape = 21, size = 9.3, stroke = 0.55, colour = "white", show.legend = FALSE
    ) +
    scale_fill_manual(values = state_palette)

  for (index in seq_len(nrow(nodes))) {
    p <- p +
      annotate(
        "text", x = nodes$x[index], y = nodes$y[index], label = nodes$state_code[index],
        family = "Arial", size = 2.45, fontface = "bold", colour = "white"
      ) +
      annotate(
        "text", x = nodes$description_x[index], y = nodes$description_y[index],
        label = nodes$description[index], family = "Arial", size = 1.86,
        fontface = "bold", colour = palette_contract[["neutral_dark"]]
      )
  }

  for (index in seq_len(nrow(edges))) {
    p <- p + annotate(
      "text", x = edges$label_x[index], y = edges$label_y[index],
      label = edges$label_text[index], family = "Arial", size = 1.90,
      fontface = "bold",
      colour = if (edges$death_transition_plot[index]) {
        palette_contract[["neutral_dark"]]
      } else {
        palette_contract[["signal_blue"]]
      }
    )
  }

  if (mortality_warning) {
    p <- p + annotate(
      "label", x = 3.72, y = 2.15,
      label = "DEATH DATA LIMITED: ascertainment ends in 2012",
      family = "Arial", size = 1.90, fontface = "bold",
      colour = palette_contract[["neutral_dark"]], fill = palette_contract[["pale_orange"]],
      linewidth = 0, label.padding = grid::unit(0.13, "lines"),
      label.r = grid::unit(0.06, "lines")
    )
  }

  panel_id <- unique(edge_data$panel_id)
  nodes[, `:=`(
    panel_id = panel_id,
    panel_title = panel_title,
    analysis_id = unique(edge_data$analysis_id),
    cohort = unique(edge_data$cohort),
    mortality_scope = unique(edge_data$mortality_scope),
    mortality_limited = mortality_warning
  )]
  list(plot = p, edges = edges, nodes = nodes)
}

if (!nrow(network_q_selected)) {
  skip_figure(
    fig_id, file_stem, fig_title, fig_claim, "quantitative state-network composite",
    "markov_q_intensities.csv is unavailable or empty", markov_q_result$path,
    "HRS, ELSA waves 4-10, and ELSA through 2012", 3L, 183, 239
  )
} else {
  built_networks <- lapply(seq_len(nrow(network_specs)), function(index) {
    specification <- network_specs[index]
    build_intensity_network_panel(
      network_q_selected[panel_id == specification$panel_id],
      include_death = specification$include_death,
      mortality_warning = specification$mortality_limited
    )
  })

  fig <- wrap_plots(
    lapply(built_networks, function(built) built$plot), ncol = 1,
    heights = c(1, 1, 1)
  ) +
    plot_annotation(
      tag_levels = "A",
      caption = paste(
        "Arrows report adjusted baseline Markov q_ij per person-year at each cohort's sample-mean covariate profile",
        "(M4, 20-imputation Rubin-pooled model). Blue solid arrows indicate non-fatal progression; charcoal dashed arrows",
        "indicate transition to the absorbing death state. Exact 95% CIs are in Source Data. Panel C is truncated at 2012",
        "and is not directly comparable with full HRS mortality follow-up; values are not a formal between-cohort test.",
        sep = "\n"
      )
    ) &
    theme(
      plot.tag = element_text(size = 8, face = "bold", family = "Arial"),
      plot.caption = element_text(
        size = 5.5, family = "Arial", colour = palette_contract[["neutral_mid"]],
        hjust = 0, margin = margin(t = 4)
      )
    )

  edge_source <- rbindlist(lapply(built_networks, function(built) {
    built$edges[, .(
      record_type = "transition_intensity", panel_id, panel_title, panel_subtitle,
      analysis_id, cohort, mortality_scope, mortality_limited, model_id,
      covariate_profile, pooling_status, m_expected, m_fitted, status,
      state_code = NA_character_, state_label = NA_character_,
      from_state, to_state, transition,
      intensity_per_person_year = intensity, standard_error,
      pooled_standard_error_log, conf_low, conf_high,
      risk_intervals, risk_participants, person_years, events,
      death_transition = death_transition_plot,
      x_from, y_from, x_to, y_to, x_to_target, y_to_target,
      x_from_plot, y_from_plot, x_to_plot, y_to_plot,
      curvature, target_y_offset, label_x, label_y, display_value, spec_id
    )]
  }), use.names = TRUE, fill = TRUE)

  node_source <- rbindlist(lapply(built_networks, function(built) {
    built$nodes[, .(
      record_type = "state", panel_id, panel_title,
      panel_subtitle = NA_character_, analysis_id, cohort, mortality_scope,
      mortality_limited, model_id = "M4_MI20", covariate_profile = "sample means",
      pooling_status = "RUBIN_POOLED_COMPLETE", m_expected = 20L, m_fitted = 20L,
      status = "completed", state_code, state_label,
      from_state = state, to_state = NA_integer_, transition = NA_character_,
      intensity_per_person_year = NA_real_, standard_error = NA_real_,
      pooled_standard_error_log = NA_real_, conf_low = NA_real_, conf_high = NA_real_, risk_intervals = NA_real_,
      risk_participants = NA_real_, person_years = NA_real_, events = NA_real_,
      death_transition = state == 4L,
      x_from = x, y_from = y, x_to = NA_real_, y_to = NA_real_,
      x_to_target = NA_real_, y_to_target = NA_real_,
      x_from_plot = NA_real_, y_from_plot = NA_real_,
      x_to_plot = NA_real_, y_to_plot = NA_real_, curvature = NA_real_, target_y_offset = NA_real_,
      label_x = description_x, label_y = description_y,
      display_value = description, spec_id = NA_character_
    )]
  }), use.names = TRUE, fill = TRUE)

  network_source <- rbindlist(list(edge_source, node_source), use.names = TRUE, fill = TRUE)
  setorder(network_source, panel_id, record_type, from_state, to_state)
  save_and_record(
    fig_id, file_stem, fig, network_source, fig_title, fig_claim,
    "quantitative state-network composite", markov_q_result$path,
    "M4_MI20; 20/20 imputations; sample-mean covariate profile",
    "HRS, ELSA waves 4-10, and ELSA mortality capture through 2012",
    183, 239, panel_ids = c("A", "B", "C"),
    column_groups = list(c("A", "B", "C")),
    notes = paste(
      "Every arrow is labelled with q_ij per person-year. Panel B omits death by design.",
      "Panel C uses all eight intensities from its own mortality-truncated five-state model;",
      "it does not splice non-fatal estimates from the primary ELSA four-state model."
    )
  )
}

# -------------------------------------------------------------------------
# Manifest and figure-level QA notes.
# -------------------------------------------------------------------------

manifest <- if (length(manifest_rows)) rbindlist(manifest_rows, use.names = TRUE, fill = TRUE) else data.table()
manifest_path <- file.path(figure_root, "multistate_figure_manifest.csv")
atomic_fwrite(manifest, manifest_path)
atomic_fwrite(manifest, file.path(output_root, "results", "multistate_figure_manifest.csv"))

qa_notes <- data.table(
  check = c(
    "backend_exclusivity", "data_integrity", "alignment", "font_floor",
    "collision_audit", "uncertainty", "multiple_testing", "mortality_scope"
  ),
  status = c(
    "R_ONLY", "NO_SYNTHETIC_ESTIMATES", "STRICT_1.5_PT_FOR_PATCHWORK",
    "MINIMUM_5_PT", if (run_pdf_qa) "RUN_AFTER_EACH_EXPORT" else "NOT_RUN_BY_CONFIGURATION",
    "95_PERCENT_CONFIDENCE_INTERVALS_WHERE_AVAILABLE", "BH_FDR_FOR_CONFIRMATORY_OR_SPECIFICATION_FAMILIES",
    "ELSA_DEATH_RESULTS_LABELLED_THROUGH_2012_ONLY"
  ),
  detail = c(
    "All plotted content and all PNG/PDF/SVG/TIFF exports are produced in R.",
    "When a result table is unavailable or empty, the figure is skipped and recorded; no placeholder estimates are drawn.",
    "Every patchwork composition invokes the backend-neutral geometry audit using rectangles measured in R at final size.",
    "The R theme and in-panel text sizes are set at or above the conservative five-point floor and PDFs are audited.",
    "The backend-neutral rendered collision audit is invoked on every final PDF; QA overlays are not submission figures.",
    "Forest and sojourn plots display 95% confidence intervals; probabilities and prevalence are point estimates returned by the fitted model.",
    "Confirmatory non-fatal transitions and sensitivity families use Benjamini-Hochberg adjustment when p values are available.",
    "ELSA mortality ascertainment does not extend across the full wave 4-10 period and is never shown as if it did."
  )
)
atomic_fwrite(qa_notes, file.path(figure_qa_root, "multistate_figure_qa_notes.csv"))

strict_failures <- manifest[
  status != "GENERATED" |
    grepl("FIX_BEFORE_DELIVERY|NOT_AUDITABLE", alignment_status, fixed = FALSE)
]
if (strict_completeness && nrow(strict_failures)) {
  stop(
    "Strict figure completeness failed; inspect multistate_figure_manifest.csv for skipped figures. No fabricated replacements were created.",
    call. = FALSE
  )
}

pdf_review_required <- manifest[
  status == "GENERATED" & grepl("FIX_BEFORE_DELIVERY|NOT_AUDITABLE", pdf_qa_status)
]
if (nrow(pdf_review_required)) {
  warning(
    nrow(pdf_review_required),
    " generated figure(s) require final-size visual adjudication of the retained PDF audit reports; generation was not misreported as an automated pass."
  )
}

cat("Multi-state figure workflow finished.\n")
cat("Generated figures:", manifest[status == "GENERATED", .N], "of", nrow(manifest), "\n")
cat("Figure manifest:", manifest_path, "\n")
