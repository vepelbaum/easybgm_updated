
#' Extract results from a Bayesian analysis of networks
#'
#' @param fit fit object of the respective package used
#' @param method type of model estimated, e.g., ggm, gcgm, dgm-binary, Ising
#' @param package package used to obtain the fit object
#' @param posterior_samples binary indicating whether the posterior samples should be extracted
#' @param not.cont only if method = "gcgm" vector indicating the not-continuous variables
#' @param data if posterior_samples = T, provide the raw data used to estimate the network
#'
#' @export
#'
extract_results <- function(fit, method, package = "BDgraph", posterior_samples=F, not.cont=NULL, data=NULL){

  # ----------------------------------------
  # Obtain the output
  # ----------------------------------------
  if(package=="BDgraph"){
    if(method %in% c("ggm")){
      bdgraph_res <- list()
      #Bayesian model-averaged estimates
      bdgraph_res$estimates_bma <- pr2pc(fit$K_hat)
      diag(bdgraph_res$estimates_bma) <- 0
      bdgraph_res$inc_probs <- as.matrix(BDgraph::plinks(fit))
      bdgraph_res$inc_probs  <- bdgraph_res$inc_probs + t(bdgraph_res$inc_probs)
      bdgraph_res$BF <- bdgraph_res$inc_probs / (1 - bdgraph_res$inc_probs)
      bdgraph_res$structure_bma <- 1*(bdgraph_res$inc_probs > 0.5)
      bdgraph_res$structure_probabilities <- fit$graph_weights/sum(fit$graph_weights)
      bdgraph_res$graph_weights <- fit$graph_weights
      bdgraph_res$sample_graph <- fit$sample_graphs
      bdgraph_res$package <- "BDgraph"
      bdgraph_res$model <- "ggm"

      if(posterior_samples == TRUE){
        if(is.null(data)){
          stop("Provide the raw data with the \"data\" argument",
               call. = FALSE)
        }
        # Extract posterior samples
        data<-as.matrix(data)
        bdgraph_res$samples_posterior <- extract_posterior(fit, data=data, method = method, not.cont)[[1]]
        # Centrality indices
        bdgraph_res$centrality_strength <- centrality_strength(bdgraph_res)
        bdgraph_res$centrality <- centrality(bdgraph_res)
      }

      output <- bdgraph_res
    }

    if(method %in% c("gcgm")){
      bdgraph_res <- list()
      #Bayesian model-averaged estimates
      bdgraph_res$estimates_bma <- pr2pc(fit$K_hat)
      diag(bdgraph_res$estimates_bma) <- 0
      bdgraph_res$inc_probs <- as.matrix(BDgraph::plinks(fit))
      bdgraph_res$inc_probs  <- bdgraph_res$inc_probs + t(bdgraph_res$inc_probs)
      bdgraph_res$BF <- bdgraph_res$inc_probs / (1 - bdgraph_res$inc_probs)
      bdgraph_res$structure_bma <- 1*(bdgraph_res$inc_probs > 0.5)
      bdgraph_res$structure_probabilities <- fit$graph_weights/sum(fit$graph_weights)
      bdgraph_res$graph_weights <- fit$graph_weights
      bdgraph_res$sample_graph <- fit$sample_graphs
      bdgraph_res$package <- "BDgraph"
      bdgraph_res$model <- "gcgm"

      if(posterior_samples == TRUE){
        if(is.null(not.cont)){
          stop("Specify a vector indicating variables are continuos with the not.cont argument (1 indicates not continuous)",
               call. = FALSE)
        }
        if(is.null(data)){
          stop("Provide the raw data with the \"data\" argument",
               call. = FALSE)
        }
        data<-as.matrix(data)
        # Extract posterior samples
        bdgraph_res$samples_posterior <- extract_posterior(fit, data, method = method, not.cont = not.cont)[[1]]
        # Centrality indices
        bdgraph_res$centrality_strength <- centrality_strength(bdgraph_res)
        bdgraph_res$centrality <- centrality(bdgraph_res)
      }
      output <- bdgraph_res
    }
    if(method %in% c("dgm-binary")){
      bdgraph_res <- list()
      #Bayesian model-averaged estimates
      bdgraph_res$inc_probs <- as.matrix(BDgraph::plinks(fit))
      bdgraph_res$inc_probs  <- bdgraph_res$inc_probs + t(bdgraph_res$inc_probs)
      bdgraph_res$BF <- bdgraph_res$inc_probs / (1 - bdgraph_res$inc_probs)
      bdgraph_res$package <- "BDgraph"
      bdgraph_res$model <- "dgm-binary"

      if(posterior_samples == TRUE){
        stop("Posterior samples cannot be obtained for \"dgm-binary\".",
             call. = FALSE)
      }
      output <- bdgraph_res
    }
  }
  if(package == "rbinnet" & method == "Ising"){
    rbinnet_res <- list()
    rbinnet_res$estimates_bma <- vector2matrix(fit$parameters$sigma_eap, fit$nodes)
    diag(rbinnet_res$estimates_bma) <- 0
    rbinnet_res$inc_probs <- vector2matrix(fit$parameters$inclusion_probabilities, fit$nodes)
    rbinnet_res$BF <- vector2matrix(fit$parameters$inc_BF, fit$nodes)
    rbinnet_res$structure_bma <- 1*(rbinnet_res$inc_probs > 0.5)
    rbinnet_res$samples_posterior <- fit$parameters$sigma_samples
    rbinnet_res$package <- "rbinnet"

    # Centrality indices
    rbinnet_res$centrality_strength <- centrality_strength(rbinnet_res)
    if(centrality_samples == TRUE){
      rbinnet_res$centrality <- centrality(rbinnet_res)
    }
    output <- rbinnet_res
  }

  return(output)
}
