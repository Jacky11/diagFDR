
#' Write diagnostic outputs to a folder
#'
#' Writes tables as CSV and plots as PNG, optionally creates a PPTX (requires officer+rvg).
#'
#' @param diag Output of \code{\link{dfdr_run_all}}.
#' @param out_dir Output directory.
#' @param formats Character vector subset of \code{c("csv","png","pptx")}.
#' @param pptx_path PPTX path (only used if \code{"pptx"} in formats).
#' @param width,height Plot size passed to \code{ggplot2::ggsave}.
#' @param dpi DPI passed to \code{ggplot2::ggsave}.
#'
#' @return Invisibly returns \code{TRUE}.
#' @export
dfdr_write_report <- function(diag,
                               out_dir,
                               formats = c("csv", "png"),
                               pptx_path = file.path(out_dir, "dfdr_plots.pptx"),
                               width = 7.2, height = 4.4, dpi = 180) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  formats <- unique(formats)

  # ---- CSV tables
  if ("csv" %in% formats) {
    tabs <- diag$tables
    # tables that are data frames
    for (nm in names(tabs)) {
      obj <- tabs[[nm]]
      if (is.data.frame(obj)) {
        readr::write_csv(obj, file.path(out_dir, paste0(nm, ".csv")))
      } else if (is.list(obj)) {
        # list of tables (e.g. sumpep per list)
        subdir <- file.path(out_dir, nm)
        dir.create(subdir, showWarnings = FALSE, recursive = TRUE)
        for (k in names(obj)) {
          if (is.data.frame(obj[[k]])) {
            readr::write_csv(obj[[k]], file.path(subdir, paste0(k, ".csv")))
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
      ggplot2::ggsave(file.path(pdir, paste0(nm, ".png")),
                      diag$plots[[nm]], width = width, height = height, dpi = dpi)
    }
  }

  # ---- PPTX (optional; requires officer+rvg)
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
