bgm_fit <- function(data, type, package = NULL, not.cont = NULL, iter = 1e4, save = FALSE, centrality = FALSE, progress = TRUE, ...){

  if(is.null(type)){
    stop("Please specify the data type (i.e., continuous, mixed, binary, ordinal) with the 'type' argument.",
         call. = FALSE)
  }
  if(type == "mixed" & is.null(not.cont)){
    stop("Please provide a binary vector of length p specifying the not continuous variables
         (1 = not continuous, 0 = continuous)",
         call. = FALSE)
  }


  # Set default values for fitting if package is unspecified
  if(is.null(package)){
    if(type == "continuos") package <- "BDgraph"
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

  #----------------------------------------------------------------------------------
  # fitting with BDgraph

  if(package == "BDgraph"){
    if(type == "continuous"){
      bdgraph_fit <- BDgraph::bdgraph(data=data,               #(M) n*p matrix of responses
                                      method="ggm",           #(M) type of data
                                      algorithm="rjmcmc",      #(O) type of sampling algorithm
                                      iter=iter,           #(O) no. iterations sampler
                                      save=save,               #(O) Should samples be stored
                                      burnin=iter/10,      #(O) no. burnin iterations sampler
                                      g.start = "empty",       #(O) starting point of graph
                                      df.prior = 3,            #(M) prior: degree of freedom for G-Wishart distribution
                                      g.prior = 0.5)           #(M) prior: inclusion probability for edges

      # extracting results
      fit <- bgm_extract(bdgraph_fit, method = "ggm",
                         package = "BDgraph", posterior_samples = save,
                         data = data, centrality = centrality)
    }
    if(type %in% c("mixed", "ordinal")){
      if(type == "ordinal") not.cont <- rep(1, ncol(data))
      # fitting the model
      bdgraph_fit <- BDgraph::bdgraph(data=data,               #(M) n*p matrix of responses
                                      method="gcgm",           #(M) type of data
                                      algorithm="rjmcmc",      #(O) type of sampling algorithm
                                      iter=iter,           #(O) no. iterations sampler
                                      save= save,               #(O) Should samples be stored
                                      burnin=iter/10,      #(O) no. burnin iterations sampler
                                      g.start = "empty",       #(O) starting point of graph
                                      not.cont = c(rep(0, 3), rep(1, 10)), #(O) Specifies not continuous variables if method is gcgm
                                      df.prior = 3,            #(M) prior: degree of freedom for G-Wishart distribution
                                      g.prior = 0.5)           #(M) prior: inclusion probability for edges

      # extracting results
      fit <- bgm_extract(bdgraph_fit, method = "gcgm",
                         package = "BDgraph", posterior_samples = save,
                         not.cont = not.cont,
                         data = data, centrality = centrality)

    }
    if(type == "binary"){
      bdgraph_fit <- bdgraph.mpl(data,
                                 method = "dgm-binary",
                                 iter = iter,
                                 save = save)
    }
    fit <- bgm_extract(bdgraph_fit, method = "dgm-binary",
                       package = "BDgraph", posterior_samples = save
                       data = data, centrality = centrality)
  }

  #----------------------------------------------------------------------------------
  # fitting with BGGM
  if(package == "BGGM"){
    # Fit the model
    bggm_fit <- BGGM::explore(data,                        #(M) n*p matrix of responses
                              type = type,                 #(O) type of data
                              mixed_type = not.cont,       #(O) which data should be treated as ranks
                              prior_sd = 0.25,             #(M) prior distribution standard deviation
                              iter = iter,             #(O) no. iterations sampler
                              progress = progress,            #(O) Should a progress bar be plotted?
                              impute = FALSE,              #(O) Should missings be imputed?
                              seed = NULL)                    #(O) Integer for random seed

    # Extracting results
    fit <- bgm_extract(bggm_fit, package = "BGGM",
                       posterior_samples = save, centrality = centrality)
  }


  #----------------------------------------------------------------------------------
  # fitting with bgms
  if(package == "bgms"){


    bgms_fit <- bgm(x = data,    #(M) n * p matrix of binary responses
                    iter = iter,        #(O) no. iterations Gibbs sampler
                    burnin = iter/10,           #(O) no. burnin iterations Gibbs sampler
                    interaction_prior = c("UnitInfo"), # (O) Type of prior distribution on the interactions
                    cauchy_scale = 2.5,     #(O) if interaction_prior is "Cauchy", the scale of the cauchy distribution
                    threshold_alpha = 1,    #(O) prior threshold parameter for the shape of the Beta-distribution
                    threshold_beta = 1,     #(O) prior threshold parameter for the shape of the Beta-distribution
                    save = save,            #(O) if TRUE, outputs posterior draws
                    display_progress = progress
    )

    fit <- bgm_extract(bgms_fit, package = "bgms", posterior_samples = save, centrality = centrality)
  }

  class(fit) <- "easybgm"
  return(fit)
}
