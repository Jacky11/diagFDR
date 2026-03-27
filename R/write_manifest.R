#' Write a manifest file describing parameters, list metadata, and warnings
#'
#' Writes a manifest alongside exported outputs to record key run parameters and
#' per-list metadata (unit, scope, q-source, export ceiling, etc.). This helps
#' tie diagnostic outputs to the exact inputs and settings used.
#'
#' If \pkg{jsonlite} is available (suggested dependency), the manifest is written
#' as JSON (\code{*.json}); otherwise a human-readable plain text file
#' (\code{*.txt}) is written.
#'
#' @param diag A list as returned by \code{\link{dfdr_run_all}}.
#' @param out_dir Character scalar. Output directory (created if it does not exist).
#' @param filename Character scalar. Base filename without extension (default \code{"manifest"}).
#'
#' @return
#' The path to the manifest file, returned invisibly. The function is called for
#' its side effect of writing a file to disk.
#'
#' @examples
#' library(tibble)
#'
#' # Create a minimal dfdr_run_all-like object to demonstrate manifest writing
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
#'   warnings = character()
#' )
#'
#' out_dir <- tempdir()
#' dfdr_write_manifest(diag, out_dir = out_dir, filename = "manifest_example")
#'
#' @export
dfdr_write_manifest <- function(diag, out_dir, filename = "manifest") {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # Collect list meta
  obj_meta <- lapply(diag$objects %||% list(), function(x) {
    attr(x, "meta") %||% list()
  })

  meta <- list(
    created = as.character(Sys.time()),
    package = "diagFDR",
    package_version = as.character(utils::packageVersion("diagFDR")),
    params = diag$params %||% list(),
    lists = obj_meta,
    warnings = diag$warnings %||% character()
  )

  if (requireNamespace("jsonlite", quietly = TRUE)) {
    path <- file.path(out_dir, paste0(filename, ".json"))
    jsonlite::write_json(meta, path = path, auto_unbox = TRUE, pretty = TRUE, null = "null")
  } else {
    # fall back to a readable plain text manifest
    path <- file.path(out_dir, paste0(filename, ".txt"))
    lines <- c(
      paste0("created: ", meta$created),
      paste0("package: ", meta$package),
      paste0("version: ", meta$package_version),
      "",
      "params:"
    )
    # params
    if (length(meta$params)) {
      for (nm in names(meta$params)) {
        lines <- c(lines, paste0("  ", nm, ": ", paste(meta$params[[nm]], collapse = ",")))
      }
    } else {
      lines <- c(lines, "  (none)")
    }

    lines <- c(lines, "", "lists:")
    if (length(meta$lists)) {
      for (nm in names(meta$lists)) {
        m <- meta$lists[[nm]]
        lines <- c(lines,
                   paste0("  - name: ", nm),
                   paste0("    unit: ", m$unit %||% NA),
                   paste0("    scope: ", m$scope %||% NA),
                   paste0("    q_source: ", m$q_source %||% NA),
                   paste0("    q_max_export: ", m$q_max_export %||% NA),
                   paste0("    p_source: ", m$p_source %||% NA)
        )
      }
    } else {
      lines <- c(lines, "  (none)")
    }

    lines <- c(lines, "", "warnings:")
    if (length(meta$warnings)) {
      lines <- c(lines, paste0("  - ", meta$warnings))
    } else {
      lines <- c(lines, "  (none)")
    }

    readr::write_lines(lines, path)
  }

  invisible(path)
}
