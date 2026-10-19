#' Render a human-readable HTML report from \code{dfdr_run_all} output
#'
#' Renders an HTML report that summarises key diagnostics and embeds plots
#' contained in \code{diag$plots}. The report is intended for interactive review
#' and auditable reporting.
#'
#' This function requires the \pkg{rmarkdown} package (suggested dependency) and
#' uses the built-in R Markdown template shipped with the package
#' (\code{inst/templates/dfdr_report.Rmd}).
#'
#' @param diag A list as returned by \code{\link{dfdr_run_all}}.
#' @param out_dir Character scalar. Output directory (created if it does not exist).
#' @param filename Character scalar. Output HTML filename (default
#'   \code{"diagFDR_report.html"}).
#' @param self_contained Logical. If \code{TRUE}, embed resources into a single HTML
#'   file (may be larger). Default \code{FALSE}.
#' @param open Logical. If \code{TRUE}, open the report in the default browser
#'   after rendering.
#'
#' @return
#' The path to the rendered HTML file, returned invisibly. The function is called
#' for its side effect of creating an HTML report on disk.
#'
#' @examples
#' # A minimal example that renders a report from a toy dataset.
#' # This example is conditional because rmarkdown is in Suggests.
#' if (requireNamespace("rmarkdown", quietly = TRUE)) {
#'   library(tibble)
#'   tmpdir <- tempdir()
#'
#'   set.seed(1)
#'   n <- 3000
#'   df <- tibble(
#'     id = as.character(seq_len(n)),
#'     run = "run1",
#'     is_decoy = sample(c(FALSE, TRUE), n, replace = TRUE, prob = c(0.95, 0.05)),
#'     score = rnorm(n),
#'     q = pmin(1, rank(-score) / n),
#'     pep = NA_real_
#'   )
#'   x <- as_dfdr_tbl(df, unit = "psm", scope = "global", q_source = "toy")
#'
#'   diag <- dfdr_run_all(
#'     xs = list(toy = x),
#'     alpha_main = 0.01,
#'     compute_pseudo_pvalues = FALSE
#'   )
#'
#'   # Render to a temporary directory (does not open a browser during checks)
#'   dfdr_render_report(diag, out_dir = tmpdir, open = FALSE)
#' }
#'
#' @export
dfdr_render_report <- function(diag, out_dir,
                               filename = "diagFDR_report.html",
                               self_contained = FALSE,
                               open = FALSE) {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    rlang::abort("HTML report rendering requires package 'rmarkdown' (in Suggests).")
  }
    if (!is.list(diag) || is.null(diag$tables) || is.null(diag$objects)) {
    rlang::abort("`diag` must be the output of dfdr_run_all().")
  }
  if (!is.character(out_dir) || length(out_dir) != 1 || !nzchar(out_dir)) {
    rlang::abort("`out_dir` must be a single directory path.")
  }
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  # rmarkdown resolves a relative output path against the template's folder,
  # so work with absolute paths throughout.
  out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)

  template <- system.file("templates", "dfdr_report.Rmd", package = "diagFDR")
  if (template == "") rlang::abort("Could not find template inst/templates/dfdr_report.Rmd in the installed package.")

  # Render from a private copy of the template: the installed library may be
  # read-only, and knitting must not write next to the installed file.
  work <- tempfile("diagFDR_report_")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  rmd <- file.path(work, "dfdr_report.Rmd")
  file.copy(template, rmd, overwrite = TRUE)

  out_file <- file.path(out_dir, filename)

  rmarkdown::render(
    input = rmd,
    output_file = out_file,
    output_format = rmarkdown::html_document(
      self_contained = self_contained, toc = TRUE, toc_depth = 3, theme = "readable",
      mathjax = NULL
    ),
    intermediates_dir = work,
    knit_root_dir = work,
    params = list(diag = diag, out_dir = out_dir),
    envir = new.env(parent = globalenv()),
    quiet = TRUE
  )


  if (isTRUE(open)) utils::browseURL(out_file)
  invisible(out_file)
}
