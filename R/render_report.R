#' Render a human-readable HTML report from dfdr_run_all output
#'
#' This renders an HTML report that summarizes key diagnostics
#' and embeds plots from \code{diag$plots}. It is intended for review and verifiable reporting.
#'
#' @param diag Output of \code{\link{dfdr_run_all}}.
#' @param out_dir Output directory.
#' @param filename Output HTML filename (default "diagFDR_report.html").
#' @param self_contained Logical. If TRUE, embed resources into a single HTML file
#'   (may be larger). Default FALSE.
#' @param open Logical. If TRUE, open the report in the default browser after rendering.
#' @return Path to the rendered HTML file (invisibly).
#' @export
dfdr_render_report <- function(diag, out_dir,
                               filename = "diagFDR_report.html",
                               self_contained = FALSE,
                               open = FALSE) {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    rlang::abort("HTML report rendering requires package 'rmarkdown' (in Suggests).")
  }
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  template <- system.file("templates", "dfdr_report.Rmd", package = "diagFDR")
  if (template == "") rlang::abort("Could not find template inst/templates/dfdr_report.Rmd in the installed package.")

  out_file <- file.path(out_dir, filename)

  rmarkdown::render(
    input = template,
    output_file = out_file,
    output_format = rmarkdown::html_document(
      self_contained = self_contained, toc = TRUE, toc_depth = 3, theme = "readable"
    ),
    params = list(diag = diag, out_dir = out_dir),
    envir = new.env(parent = globalenv()),
    quiet = TRUE
  )

  if (isTRUE(open)) utils::browseURL(out_file)
  invisible(out_file)
}
