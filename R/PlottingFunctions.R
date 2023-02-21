

#' Plot Posterior Structure Probabilities
#'
#' @param output Output object from the bgm_extract function
#' @param as.BF if TRUE plots the y-axis as Bayes factor instead of posterior structure probability
#'
#' @export
#' @import ggplot2
#'

plot_structure_probability <- function(output, as.BF = FALSE) {

  sorted_structure_prob <- as.data.frame(sort(output$structure_probabilities, decreasing=T))
  colnames(sorted_structure_prob) <- "posterior_prob"
  if(as.BF){
    BF1s <- sorted_structure_prob$posterior_prob[1] / sorted_structure_prob$posterior_prob # BF best structure vs. others
    data <- data.frame(structures = 1:length(BF1s), BayesFactor = BF1s)
    ggplot2::ggplot(data, aes(x = structures, y = BayesFactor)) +
      geom_point(size = 4, shape = 1) +
      scale_y_continuous(trans = "log10") +
      theme_classic()+
      labs(x = "Structures",
           y = expression(log(BF[1][s])))+
      geom_hline(yintercept = 1/10, linetype = "dashed", size = 1.5)  +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", linewidth = 1.1),
            axis.text = element_text(size = 14), axis.title = element_text(size = 16),
            axis.ticks.length = unit(.25, "cm"))
  } else {
    data <- data.frame(structures = 1:nrow(sorted_structure_prob), Probs = sorted_structure_prob)
    ggplot2::ggplot(data, aes(x = structures, y = posterior_prob)) +
      geom_point(size = 4, shape = 1) +
      theme_classic()+
      labs(x = "Structures",
           y = "Posterior Structure Probability")+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", linewidth = 1.1),
            axis.text = element_text(size = 14), axis.title = element_text(size = 16),
            axis.ticks.length = unit(.25, "cm"))
  }
}

# ---------------------------------------------------------------------------------------------------------------


#' Plot posterior structure complexity
#'
#' @param output Output object from the bgm_extract function
#'
#' @export
#' @import ggplot2
#'

plot_posteriorcomplexity <- function(output) {
  if (output$package == "rbinnet") {
    stop("Plot not implemented for rbinnet",
         call. = FALSE)
  }
  complexity <- c()
  for(i in 1:length(output$sample_graph)){
    complexity[i] <- sum(as.numeric(unlist(strsplit(res$sample_graph[i], ""))))
  }

  data_complexity <- tibble(complexity, weights = res$graph_weights)  %>%
    group_by(complexity) %>%
    summarise(complexity_weight = sum(weights)) %>%
    mutate(complexity_weight = complexity_weight/sum(complexity_weight))

  ggplot(data_complexity, aes(x = complexity, y = complexity_weight))+
    geom_point() +
    ylab("Posterior Complexity Probability") +
    xlab("Complexity")  +
    theme_minimal()+
    theme(legend.position = c(.85, 0.25), axis.text=element_text(size=20),
          legend.background = element_rect(fill = NULL), panel.border = element_blank(),
          axis.line = element_line(colour = "black", size = 1.1), axis.ticks.length=unit(.2, "cm"),
          axis.ticks = element_line(size= .8), legend.text = element_text(size=14),
          axis.title.x = element_text(size=18,face="bold"),
          axis.title.y = element_text(size=18,face="bold"),
          text=element_text(  family="Times New Roman"),
          panel.grid.major = element_blank()
    )
}

# ---------------------------------------------------------------------------------------------------------------

#' Edge evidence plot
#'
#' @param output Output object from the bgm_extract function
#' @param evidence_thresh BF which will be considered sufficient evidence for in-/exclusion
#' @param layout Layout of the network; qgraph argument
#' @param edge.width Layout of the network; qgraph argument
#' @param split if TRUE, plot is split in included and excluded edges
#' @param ... Additional `qgraph` arguments
#'
#' @export
#' @import qgraph
#'
plot_edgeevidence <- function(output, evidence_thresh = 10, split = F, ...) {

  graph <- output$BF
  diag(graph) <- 1

  # assign a color to each edge (inclusion - blue, exclusion - red, no conclusion - grey)
  graph_color <- graph
  graph_color <-  ifelse(graph < evidence_thresh & graph > 1/evidence_thresh, graph_color <- "#bfbfbf", graph_color <- "#36648b")
  graph_color[graph < (1/evidence_thresh)] <- "#990000"

  if(split == F){
    graph[output$inc_probs <= 1] <- 1
    diag(graph) <- 1
    colnames(graph) <- colnames(output$sigma)
    qgraph::qgraph(graph,
                   edge.color = graph_color, # specifies the color of the edges
                   ...
    )
  }
  if(split==T){
    par(mfrow=c(1, 2))
    graph_inc <- graph_exc <- graph
    # plot included graph
    graph_inc[output$inc_probs >= .5] <- 1
    graph_inc[output$inc_probs < .5] <- 0
    diag(graph_inc) <- 1
    colnames(graph_inc) <- colnames(output$sigma)
    qgraph::qgraph(graph_inc,
                   edge.color = graph_color, # specifies the color of the edges
                   ...
    )
    # Plot excluded graph
    graph_exc[output$inc_probs >= .5] <- 0
    graph_exc[output$inc_probs < .5] <- 1
    diag(graph_exc) <- 1
    colnames(graph_exc) <- colnames(output$sigma)
    qgraph::qgraph(graph_exc,
                   edge.color = graph_color, # specifies the color of the edges
                   ...
    )
    par(mfrow=c(1, 1))
  }

}

# ---------------------------------------------------------------------------------------------------------------

#' Network plot
#'
#' @param output Output object from the bgm_extract function
#' @param exc_prob threshold for excluding edges; all edges with a lower inclusion probability will not be shown
#' @param ... Additional `qgraph` arguments
#'
#' @export
#' @import qgraph

plot_network <- function(output, exc_prob = .5, ...) {


  graph <- output$sigma

  # Exclude edges with a inclusion probability lower .5
  inc_probs_m <- output$inc_probs
  graph[inc_probs_m < exc_prob] <- 0
  diag(graph) <- 1

  # Plot
  qgraph::qgraph(graph, ...)

}

# ---------------------------------------------------------------------------------------------------------------

#' Plot of interaction parameters and their 95% highest density intervals
#'
#' @param output Output object from the extract_results function
#'
#' @export
#' @import ggplot2 HDInterval
#'
plot_parameterHDI <- function(output) {

  package <- output$package
  if(is.null(output$samples_posterior)){
    stop("Samples of the posterior distribution required. When extracting the results, set \"posterior_samples = TRUE\".")
  }

  hdi_intervals <- as.data.frame(apply(output$samples_posterior, MARGIN = 2, FUN = hdi))
  posterior_medians <- apply(output$samples_posterior, MARGIN = 2, FUN = median)

  names <- colnames(output$sigma)
  combinations <- combn(nrow(output$sigma), 2)
  index <- vector(length = ncol(combinations))
  for(i in 1:ncol(combinations)) index[i] <- paste0(names[combinations[1, i]],"-",names[combinations[2, i]])

  posterior <- cbind(data.frame(posterior_medians, row.names = NULL),
                     data.frame(t(hdi_intervals), row.names = NULL), index)
  colnames(posterior) <- c("posterior_medians", "lower", "upper", "names")
  posterior <- posterior[order(posterior$posterior_medians, decreasing = F),]
  posterior$names <- factor(posterior$names, levels = posterior$names)


  ggplot2::ggplot(data = posterior, aes(x = names, y = posterior_medians, ymin = lower,
                                        ymax = upper)) +
    geom_pointrange(position=position_dodge(width=c(0.3)), size = .5) +
    theme_bw() +
    coord_flip() +
    ylab("Highest Density Interval of Parameter")+
    xlab("") +
    geom_hline(yintercept = 0, linetype = "dashed", size = 1.3) +
    theme(axis.text=element_text(size=8), panel.border = element_blank(),
          axis.line = element_line(colour = "black", size = 1.1), axis.ticks.length=unit(.2, "cm"),
          axis.ticks = element_line(size= .8),
          axis.title.x = element_text(size=16,face="bold"), plot.title = element_text(size = 18, face = "bold"))
}

