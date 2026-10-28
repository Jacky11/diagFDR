#' Simulated competed target/decoy list (example data)
#'
#' A \code{dfdr_tbl} of 20,000 units produced by
#' \code{dfdr_simulate_tdc(n = 20000, pi0 = 0.667, delta = 4.65, seed = 20260425)}:
#' independent Gaussian nulls, correct specification. The attribute
#' \code{"truth"} records which units are incorrect, so realised FDPs can be
#' computed with \code{\link{dfdr_realised_fdp}}. Used by the examples and
#' the quick-start vignette so that they run without external data.
#'
#' @format A \code{dfdr_tbl} with columns \code{id}, \code{is_decoy},
#'   \code{q}, \code{pep} (true PEP of targets), \code{run}, \code{score}.
#' @source \code{\link{dfdr_simulate_tdc}}
"sim_tdc_example"
