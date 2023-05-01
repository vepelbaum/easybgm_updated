#' Fit a Bayesian analysis of networks
#'
#' @param data An n x p matrix or dataframe containing the variables for n independent observations on p variables.
#' @param type What is the data type? Options: continuos, mixed, ordinal, binary
#' @param package The R-package that should be used for fitting the network model. Optional argument;
#'     default values are specified depending on the datatype.
#' @param not.cont If data-type is mixed, a vector of length p, specifying the not-continuous
#'     variables (1 = not continuous, 0 = continuous)
#' @param iter number of iterations for the sampler
#' @param save Logical. Should the posterior samples be obtained (default = FALSE)?
#' @param centrality Logical. Should the centrality measures be extracted (default = FALSE)? Note, that it will significantly increase the computation time.
#' @param progress Logical. Should a progress bar be shown (default = TRUE)?
#' @param edge.prior Default is 0.5. Single value or p x p matrix with one value per edge.
#' @param ... Additional arguments that are handed to the fitting functions of the packages, e.g., informed prior specifications.
#'
#' @return
#' @export
#' @import bgms
#' @import BDgraph
#' @import BGGM
#'
#' @examples
easybgm <- function(data, type, package = NULL, not.cont = NULL, iter = 1e4,
                    save = FALSE, centrality = FALSE, progress = TRUE,
                    edge.prior = 0.5, ...){


  if(type == "mixed" & is.null(not.cont)){
    stop("Please provide a binary vector of length p specifying the not continuous variables
         (1 = not continuous, 0 = continuous).",
         call. = FALSE)
  }


  # Set default values for fitting if package is unspecified
  if(is.null(package)){
    if(type == "continuous") package <- "BDgraph"
    if(type == "mixed") package <- "BDgraph"
    if(type == "ordinal") package <- "bgms"
    if(type == "binary") package <- "bgms"
  }

  if((package == "bgms") & (type %in% c("continuous", "mixed"))){
    warning("bgms can only fit ordinal or binary datatypes. For continuous or mixed data,
           choose either the BDgraph or BGGM package. By default we have changed the package to BDgraph",
            call. = FALSE)
    package <- "BDgraph"
  }


  fit <- list()
  class(fit) <- c(package, "easybgm")

  # Fit the model
  tryCatch(
    {fit <- bgm_fit(fit, data = data, type = type, not.cont = not.cont, iter = iter,
                    save = save, centrality = centrality, progress = progress, ...)
    },
    error = function(e){
      # If an error occurs, stop running the code
      stop("Error: ", e$message)
    })

  # Extract the results
  res <- bgm_extract(fit, model = fit$model,
                     save = save,
                     data = data, centrality = centrality, ...)

  # Output results
  class(res) <- c(package, "easybgm")
  return(res)
}
