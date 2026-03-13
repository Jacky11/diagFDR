#' Write a human-readable README.md report
#'
#' Produces a lightweight narrative report (Markdown) mapping the export folder
#' to key diagnostics (scope/unit, stability, calibration) and linking plots/tables.
#'
#' @param diag Output of \code{\link{dfdr_run_all}}.
#' @param out_dir Output directory.
#' @param filename Markdown filename (default "README.md").
#' @return Invisibly returns the path to the README.
#' @export
dfdr_write_readme <- function(diag, out_dir, filename = "README.md") {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  path <- file.path(out_dir, filename)

  summ <- dfdr_summary_headline(diag)

  md_table <- function(df, digits = 4) {
    if (!is.data.frame(df) || nrow(df) == 0) return(c("(no rows)"))
    df_chr <- as.data.frame(lapply(df, function(z) {
      if (is.numeric(z)) format(z, digits = digits, scientific = TRUE) else as.character(z)
    }), stringsAsFactors = FALSE)
    header <- paste0("| ", paste(names(df_chr), collapse = " | "), " |")
    sep    <- paste0("| ", paste(rep("---", ncol(df_chr)), collapse = " | "), " |")
    rows   <- apply(df_chr, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
    c(header, sep, rows)
  }

  alpha_main <- (diag$params %||% list())$alpha_main %||% NA_real_
  lines <- c(
    "# diagFDR report",
    "",
    paste0("- Created: ", as.character(Sys.time())),
    paste0("- diagFDR version: ", as.character(utils::packageVersion("diagFDR"))),
    paste0("- Headline threshold (alpha_main): ", alpha_main),
    "",
    "## Headline summary (Table 6 essentials)",
    "",
    md_table(summ),
    "",
    "## Plots",
    "",
    "Plots are stored in `plots/` (PNG). Key plots include:",
    "- `plots/dalpha.png` (decoy support across alpha)",
    "- `plots/cv.png` (instability proxy across alpha)",
    "- `plots/dwin.png` (local boundary decoy support across alpha)",
    "- `plots/elasticity.png` (Jaccard stability across alpha)",
    "",
    "## Exported tables",
    "",
    "All diagnostic tables are exported as CSV in this folder. The recommended single-file summary is:",
    "- `summary_headline.csv`",
    "",
    "## Metadata",
    ""
  )

  # list metadata
  if (!is.null(diag$objects) && length(diag$objects)) {
    for (nm in names(diag$objects)) {
      meta <- attr(diag$objects[[nm]], "meta") %||% list()
      lines <- c(lines,
                 paste0("### ", nm),
                 paste0("- unit: ", meta$unit %||% "NA"),
                 paste0("- scope: ", meta$scope %||% "NA"),
                 paste0("- q_source: ", meta$q_source %||% "NA"),
                 paste0("- q_max_export: ", meta$q_max_export %||% "NA"),
                 if (!is.null(meta$p_source) && !is.na(meta$p_source)) paste0("- p_source: ", meta$p_source) else NULL,
                 ""
      )
    }
  } else {
    lines <- c(lines, "(No objects found in `diag$objects`.)", "")
  }

  warn <- diag$warnings %||% character()
  lines <- c(lines, "## Warnings / notes", "")
  lines <- c(lines, if (length(warn)) paste0("- ", warn) else "- (none)", "")

  readr::write_lines(lines, path)
  invisible(path)
}
