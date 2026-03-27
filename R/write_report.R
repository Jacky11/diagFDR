#' Write diagnostic outputs to a folder
#'
#' Exports the output of \code{\link{dfdr_run_all}} to disk. Depending on selected
#' \code{formats}, this function can:
#' \itemize{
#'   \item write diagnostic tables as CSV;
#'   \item write plots as PNG;
#'   \item create a PowerPoint (\code{.pptx}) containing the plots (requires \pkg{officer} and \pkg{rvg} in Suggests);
#'   \item write a lightweight manifest (\code{\link{dfdr_write_manifest}});
#'   \item write a human-readable \code{README.md} (\code{\link{dfdr_write_readme}});
#'   \item write a compact headline summary table (\code{summary_headline.csv}).
#' }
#'
#' @param diag A list as returned by \code{\link{dfdr_run_all}}.
#' @param out_dir Character scalar. Output directory (created if it does not exist).
#' @param formats Character vector. Subset of
#'   \code{c("csv","png","pptx","manifest","readme","summary")}.
#' @param pptx_path Character scalar. Output path for the PPTX file (only used if
#'   \code{"pptx"} is included in \code{formats}).
#' @param width,height Numeric. Plot size (in inches) passed to \code{ggplot2::ggsave()}
#'   when writing PNGs.
#' @param dpi Numeric. DPI passed to \code{ggplot2::ggsave()}.
#'
#' @return
#' Returns \code{TRUE} invisibly. The function is called for its side effects of
#' writing files to \code{out_dir}.
#'
#' @examples
#' library(tibble)
#'
#' # Create a minimal diag object with one table and one plot
#' diag <- list(
#'   tables = list(headline = tibble(list = "toy", alpha = 0.01, D_alpha = 10)),
#'   plots = list(example_plot = ggplot2::ggplot(tibble(x = 1:3, y = 1:3),
#'                                              ggplot2::aes(x, y)) +
#'                 ggplot2::geom_point()),
#'   objects = list(),
#'   params = list(alpha_main = 0.01),
#'   warnings = character()
#' )
#'
#' out_dir <- tempdir()
#' dfdr_write_report(
#' diag, out_dir = out_dir,
#' formats = c("csv", "png", "summary", "manifest", "readme"))
#'
#' @export
dfdr_write_report <- function(diag,
                              out_dir,
                              formats = c("csv", "png", "manifest", "readme", "summary"),
                              pptx_path = file.path(out_dir, "diagFDR_report.pptx"),
                              width = 7.2, height = 4.4, dpi = 180) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  formats <- unique(formats)

  if ("manifest" %in% formats) dfdr_write_manifest(diag, out_dir)
  if ("readme" %in% formats) dfdr_write_readme(diag, out_dir)
  if ("summary" %in% formats) {
    summ <- dfdr_summary_headline(diag)
    readr::write_csv(summ, file.path(out_dir, "summary_headline.csv"))
  }

  # ---- CSV tables
  if ("csv" %in% formats) {
    tabs <- diag$tables
    for (nm in names(tabs)) {
      obj <- tabs[[nm]]
      nm2 <- dfdr_safe_filename(nm)

      if (is.data.frame(obj)) {
        readr::write_csv(obj, file.path(out_dir, paste0(nm2, ".csv")))
      } else if (is.list(obj)) {
        subdir <- file.path(out_dir, nm2)
        dir.create(subdir, showWarnings = FALSE, recursive = TRUE)
        for (k in names(obj)) {
          k2 <- dfdr_safe_filename(k)
          if (is.data.frame(obj[[k]])) {
            readr::write_csv(obj[[k]], file.path(subdir, paste0(k2, ".csv")))
          }
        }
      }
    }
  }

  # ---- PNG plots
  if ("png" %in% formats) {
    pdir <- file.path(out_dir, "plots")
    dir.create(pdir, showWarnings = FALSE, recursive = TRUE)
    for (nm in names(diag$plots)) {
      nm2 <- dfdr_safe_filename(nm)
      ggplot2::ggsave(file.path(pdir, paste0(nm2, ".png")),
                      diag$plots[[nm]], width = width, height = height, dpi = dpi)
    }
  }

  # ---- PPTX (optional)
  if ("pptx" %in% formats) {
    if (!requireNamespace("officer", quietly = TRUE) ||
        !requireNamespace("rvg", quietly = TRUE)) {
      rlang::abort("PPTX export requires packages 'officer' and 'rvg' in Suggests.")
    }
    doc <- officer::read_pptx()
    for (nm in names(diag$plots)) {
      doc <- officer::add_slide(doc, layout = "Title and Content", master = "Office Theme")
      doc <- officer::ph_with(doc, value = nm, location = officer::ph_location_type(type = "title"))
      doc <- officer::ph_with(doc, value = rvg::dml(ggobj = diag$plots[[nm]]),
                              location = officer::ph_location_type(type = "body"))
    }
    print(doc, target = pptx_path)
  }

  invisible(TRUE)
}

