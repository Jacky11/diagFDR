#' Write a human-readable \code{README.md} report
#'
#' Produces a lightweight narrative report in Markdown that summarises the most
#' important diagnostics (headline summary, key plots, exported tables) and records
#' per-list metadata (unit, scope, q-source, etc.) from a \code{\link{dfdr_run_all}}
#' result. This is intended to accompany an export folder for auditing and sharing.
#'
#' @param diag A list as returned by \code{\link{dfdr_run_all}}.
#' @param out_dir Character scalar. Output directory (created if it does not exist).
#' @param filename Character scalar. Markdown filename (default \code{"README.md"}).
#'
#' @return
#' The path to the written README file, returned invisibly. The function is called
#' for its side effect of writing a file to disk.
#'
#' @examples
#' library(tibble)
#'
#' # Minimal dfdr_run_all-like object for demonstrating README writing
#' df <- tibble(
#'   id = c("1","2","3"),
#'   run = "run1",
#'   is_decoy = c(FALSE, TRUE, FALSE),
#'   score = c(10, 9, 8),
#'   q = c(0.01, 0.02, 0.03),
#'   pep = NA_real_
#' )
#' x <- as_dfdr_tbl(df, unit = "psm", scope = "global", q_source = "toy")
#'
#' diag <- list(
#'   objects = list(toy = x),
#'   params = list(alpha_main = 0.01),
#'   tables = list(headline = tibble(list = "toy", alpha = 0.01)),
#'   warnings = character()
#' )
#'
#' out_dir <- tempdir()
#' dfdr_write_readme(diag, out_dir = out_dir, filename = "README_example.md")
#'
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
