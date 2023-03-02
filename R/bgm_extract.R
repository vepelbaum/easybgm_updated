#' Extract results from a Bayesian analysis of networks
#'
#' @param fit fit object of the respective package used. Note for objects from the package 'BGGM', the package requires the input from explore(data)
#' @param method type of model estimated, e.g., ggm, gcgm, dgm-binary, Ising
#' @param package package used to obtain the fit object
#' @param posterior_samples if TRUE, the posterior samples will be extracted. Note will significantly increase the computation time for 'BDgraph'.
#' @param not.cont only if method = "gcgm" vector indicating the not-continuous variables
#' @param data if posterior_samples = T, provide the raw data used to estimate the network
#' @param centrality if TRUE, the centrality samples will be extracted. Note will significantly increase the computation time.
#'
#' @export
#' @import BDgraph dplyr
#'
bgm_extract <- function(fit, method, edge.prior = 0.5, package = "BDgraph", posterior_samples=F, not.cont=NULL, data=NULL, centrality =F){

  # ----------------------------------------
  # Obtain the output
  # ----------------------------------------
  if(package=="BDgraph"){
    if(method %in% c("ggm")){
      bdgraph_res <- list()
      #Bayesian model-averaged estimates
      bdgraph_res$sigma <- pr2pc(fit$K_hat)
      diag(bdgraph_res$sigma) <- 0
      bdgraph_res$inc_probs <- as.matrix(BDgraph::plinks(fit))
      bdgraph_res$inc_probs  <- bdgraph_res$inc_probs + t(bdgraph_res$inc_probs)
      bdgraph_res$BF <- (bdgraph_res$inc_probs / (1 - bdgraph_res$inc_probs))/(edge.prior /(1-edge.prior))
      bdgraph_res$structure <- 1*(bdgraph_res$inc_probs > 0.5)
      bdgraph_res$structure_probabilities <- fit$graph_weights/sum(fit$graph_weights)
      bdgraph_res$graph_weights <- fit$graph_weights
      bdgraph_res$sample_graph <- fit$sample_graphs
      bdgraph_res$package <- "BDgraph"
      bdgraph_res$model <- "ggm"

      if(centrality == TRUE){
        posterior_samples <- TRUE
      }

      if(posterior_samples == TRUE){
        if(is.null(data)){
          stop("Provide the raw data with the \"data\" argument",
               call. = FALSE)
        }
        # Extract posterior samples
        data<-as.matrix(data)
        bdgraph_res$samples_posterior <- extract_posterior(fit, data=data, method = method, not.cont)[[1]]

        if(centrality == TRUE){
          # Centrality indices
          #bdgraph_res$centrality_strength <- centrality_strength(bdgraph_res)
          bdgraph_res$centrality <- centrality(bdgraph_res)
        }

      }

      output <- bdgraph_res
    }

    if(method %in% c("gcgm")){
      bdgraph_res <- list()
      #Bayesian model-averaged estimates
      bdgraph_res$sigma <- pr2pc(fit$K_hat)
      diag(bdgraph_res$sigma) <- 0
      bdgraph_res$inc_probs <- as.matrix(BDgraph::plinks(fit))
      bdgraph_res$inc_probs  <- bdgraph_res$inc_probs + t(bdgraph_res$inc_probs)
      bdgraph_res$BF <- (bdgraph_res$inc_probs / (1 - bdgraph_res$inc_probs))/(edge.prior/(1-edge.prior))
      bdgraph_res$structure <- 1*(bdgraph_res$inc_probs > 0.5)
      bdgraph_res$structure_probabilities <- fit$graph_weights/sum(fit$graph_weights)
      bdgraph_res$graph_weights <- fit$graph_weights
      bdgraph_res$sample_graph <- fit$sample_graphs
      bdgraph_res$package <- "BDgraph"
      bdgraph_res$model <- "gcgm"

      if(centrality == TRUE){
        posterior_samples <- TRUE
      }
      if(posterior_samples == TRUE){
          stop("Posterior samples cannot be extracted for GCGMs at the moment.")

        # if(is.null(not.cont)){
        #   stop("Specify a vector indicating variables are continuos with the not.cont argument (1 indicates not continuous)",
        #        call. = FALSE)
        # }
        # if(is.null(data)){
        #   stop("Provide the raw data with the \"data\" argument",
        #        call. = FALSE)
        # }
        # data<-as.matrix(data)
        # # Extract posterior samples
        # bdgraph_res$samples_posterior <- extract_posterior(fit, data, method = method, not.cont = not.cont)[[1]]
        #
        # if(centrality == TRUE){
        # # Centrality indices
        # # bdgraph_res$centrality_strength <- centrality_strength(bdgraph_res)
        # bdgraph_res$centrality <- centrality(bdgraph_res)
        # }
      }
      output <- bdgraph_res
    }
    if(method %in% c("dgm-binary")){
      bdgraph_res <- list()
      #Bayesian model-averaged estimates
      bdgraph_res$inc_probs <- as.matrix(BDgraph::plinks(fit))
      bdgraph_res$BF <- (bdgraph_res$inc_probs / (1 - bdgraph_res$inc_probs))/(edge.prior/(1-edge.prior))
      bdgraph_res$package <- "BDgraph"
      bdgraph_res$model <- "dgm-binary"

      if(posterior_samples == TRUE){
        stop("Posterior samples cannot be obtained for \"dgm-binary\".",
             call. = FALSE)
      }
      output <- bdgraph_res
    }
  }
  if(package == "BGGM") {
    out_select <- BGGM::select(fit)
    bggm_res <- list()
    bggm_res$sigma <- out_select$pcor_mat
    colnames(bggm_res$sigma) <- colnames(fit$Y)
    bggm_res$BF <- out_select$BF_10
    bggm_res$inc_probs <- out_select$BF_10/(out_select$BF_10 + 1)
    bggm_res$structure <- out_select$Adj_10

    if(centrality == TRUE){
      posterior_samples <- TRUE
    }
    if(posterior_samples == TRUE){
      p <- ncol(bggm_res$sigma)
      samples <- matrix(0, ncol = p*(p-1)/2, nrow = fit$iter)
      for(i in 1:fit$iter){
        sample <- fit$post_samp$pcors[, , i]
        samples[i, ] <- as.vector(sample[upper.tri(sample)])
      }
      bggm_res$samples_posterior <- samples

      if(centrality == TRUE){
        # bggm_res$centrality_strength <- centrality_strength(bggm_res)
        bggm_res$centrality <- centrality(bggm_res)
      }
    }
    bggm_res$package <- "BGGM"
    if(fit$type == "continuous"){
      bggm_res$model <- "ggm"
    } else {
      bggm_res$model <- "gcgm"
    }

    output <- bggm_res
  }
  if(package == "rbinnet" & method == "Ising"){
    rbinnet_res <- list()
    rbinnet_res$sigma <- vector2matrix(fit$parameters$sigma_eap, fit$nodes)
    diag(rbinnet_res$sigma) <- 0
    rbinnet_res$inc_probs <- vector2matrix(fit$parameters$inclusion_probabilities, fit$nodes)
    rbinnet_res$BF <- vector2matrix(fit$parameters$inc_BF, fit$nodes)
    rbinnet_res$structure <- 1*(rbinnet_res$inc_probs > 0.5)
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
