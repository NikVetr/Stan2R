#' Stan2R: generative R code from Stan models
#'
#' Parse a Stan program with the Stan compiler ([stan_program()]) and translate
#' it into an R script for prior or posterior predictive simulation or pointwise
#' predictive densities ([stan2r()], [parse_Stan()]); rewrite multilevel models
#' with [flatten_model()].
#'
#' @importFrom stats setNames
#' @importFrom graphics par rect segments strheight strwidth text
#' @importFrom grDevices colorRampPalette
#' @keywords internal
"_PACKAGE"
