#' Extract the results of a Bayesian analysis of networks
#'
#' @param fit Fit object with a particular class that will dispatch to the respective package functions
#' @param ... Additional arguments to be passed onto the respective fitting functions



bgm_extract <- function(fit, ...){

  UseMethod("bgm_extract", fit)

}
