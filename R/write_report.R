#' Write diagnostic outputs to a folder
#'
#' Writes tables as CSV and plots as PNG, optionally creates a PPTX (requires officer+rvg).
#' Also writes a lightweight manifest and README for human-readable reporting.
#'
#' @param diag Output of \code{\link{dfdr_run_all}}.
#' @param out_dir Output directory.
#' @param formats Character vector subset of \code{c("csv","png","pptx","manifest","readme")}.
#' @param pptx_path PPTX path (only used if \code{"pptx"} in formats).
#' @param width,height Plot size passed to \code{ggplot2::ggsave}.
#' @param dpi DPI passed to \code{ggplot2::ggsave}.
#'
#' @return Invisibly returns \code{TRUE}.
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

