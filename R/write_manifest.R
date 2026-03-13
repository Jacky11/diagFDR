#' Write a manifest file describing parameters, list metadata, and warnings
#'
#' Writes \code{manifest.json} if \pkg{jsonlite} is available, otherwise writes
#' \code{manifest.txt}. This ties exported outputs to scope/unit/q-source settings.
#'
#' @param diag Output of \code{\link{dfdr_run_all}}.
#' @param out_dir Output directory.
#' @param filename Base filename (default "manifest").
#' @return Invisibly returns the path to the manifest file.
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
