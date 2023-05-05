# easybgm: Easy Bayesian Graphical Modeling

DISCLAIMER: The package is still undergoing rapid development and changes to functions or arguments may occur. 

The `R` package `easybgm` provides a user-friendly package for performing a Bayesian analysis of psychometric networks. In particular, it helps to fit, extract, and visualize the results of a Bayesian graphical model commonly used in the social-behavioral sciences. The package is a wrapper around existing packages. So far, the package supports fitting and extracting results of cross-sectional network models using `BDgraph` (Mohammadi \& Wit, 2015), `BGGM`(Williams \& Mulder, 2019), `bgms` (Marsman \& Haslbeck, 2023). As output, the package extracts the parameter estimates, the posterior inclusion probability, the inclusion Bayes factor, and optionally posterior samples of the parameters and the centrality. The package furthermore provides an extensive suite of visualization functions. 

## Installation

To install this package from Github use

```r
install.packages("remotes")
remotes::install_github("KarolineHuth/easybgm")
```

To rather install the most up-to-date developer version, use 

```r
install.packages("remotes")
remotes::install_github("KarolineHuth/easybgm", ref = "developer")
```

## Overview

### Estimation

For the estimation, our package builds wrapper functions around existing R packages (i.e., `BDgraph`, `bgms`, and `BGGM`). To initiate estimation, researchers specify the data set and the data type (i.e., continuous, mixed, ordinal, or binary). Based on the data type specification, `easybgm` estimates the network using the appropriate R package (i.e., `BDgraph` for continuous and mixed data, and `bgms` for ordinal and binary data). Users can deviate from the default package choice by specifying their preferred R package with the `package` argument. Based on the different estimates, `easybgm` outputs the results, including the parameter estimates, the posterior inclusion probability, the inclusion Bayes factor, and optionally the posterior samples of the parameters as well as strength centrality samples.

### Visualization

The package comes with an extensive suite to visualize the results of the Bayesian analysis of networks. We provide more information on each of the plots below. The visualization function use `qgraph` (Epskamp et al., 2012) and `ggplot2` (Wickham, 2016) as the backbone for the plots. 

#### Edge Evidence Plot 

In the edge evidence plot, edges represent the inclusion Bayes factor $\text{BF}_{10}$. The edge evidence plot aids researchers in deciding which edges provide robust inferential conclusions: red edges indicate evidence for edge absence (i.e., conditional independence), grey edges indicate the absence of evidence, and blue edges indicate evidence for edge presence (i.e., conditional dependence).

#### Network Plot

In the network plot in which edges indicate the strength of partial association between two nodes. The network plot shows all edges with an inclusion Bayes factor larger than $1$. Edge thickness and saturation represent the strength of the association; the thicker the edge, the stronger the association. Red edges indicate negative relations and blue edges indicate positive associations.

#### Structure Plots

!!! EDIT!!!! 

Posterior probabilities of the visited structures, sorted from the most to
the least probable. Each dot represents one structure. The grey dashed line denotes the
prior probability for each structure. The right plot indicates the posterior probability for
different structure complexities, where complexity comprises the network density. 

#### 95 \% Highest Density Intervals of the Posterior Parameter Distribution

The 95 \% highest density interval (HDI) of the parameters is visualized with a parameter forest plot. In the plot, dots represent the median of the posterior samples and the lines indicate the shortest interval that covers 95\% of the posterior distribution. The narrower an interval, the more stable a parameter.  

#### Centrality Plot

!!! EDIT!!!! 

Researchers often use centrality measures to obtain aggregated information for each node, for example, the connectedness quantified in the strength centrality. Credible intervals for strength centrality can be obtained by calculating the centrality measure for each sample of the posterior distribution. The higher the centrality, the more strongly connected the node; error bars represent the 95\% highest density interval (HDI). 

## Example

We want to illustrate the package use with an example. In particular, we use the Wenchuan data which can be loaded with 
the package `bgms`. We fit the model and extract its results with the function `easybgm`. We specify the data and the data type, which in this case is `ordinal`. 

```r
library(easybgm)
library(bgms)

data <- Wenchuan
res <- easybgm(data, type = "ordinal")
```

Having fitted the model, we can now visualize its results. In a first step, we assess the edge evidence plot in which edges represent the inclusion Bayes factor $\text{BF}_{10}$. Especially in a large network like ours, it is useful to split the edge evidence plot in two parts by setting the `split` argument to `TRUE`. As such, the left plot shows edges with some evidence for inclusion (i.e., $\text{BF}_{10} > 1$), where blue edges represent evidence for inclusion ($\text{BF}_{10} > 10$) and grey edges absence of evidence ($1 < \text{BF}_{10} < 10$). The right edge evidence plot shows edges with some evidence for exclusion (i.e., $\text{BF}_{10} < 1$) with evidence for exclusion shown as red ($\text{BF}_{01} > 10$) and inconclusive evidence as grey ($0.1 < \text{BF}_{10} < 1$).  

```r
plot_edgeevidence(res, edge.width = 2, split = T)
```

Furthermore, we can look at the network plot in which edges represent the partial associations.

```r
plot_network(res, layout = "spring", 
             layoutScale = c(.8,1), palette = "R",
             theme = "TeamFortress", vsize = 6)
```

In addition we can obtain posterior samples from the posterior distribution by setting `save = TRUE` in the `easybgm` function and thereby open up new possibilities of assessing the model. We can extract the posterior density of the parameters with a parameter forest plot. 

```r
res <- easybgm(data, type = "ordinal", save = TRUE, centrality = TRUE)
plot_parameterHDI(res)
```

Lastly, we can assess the structure specifically with three plots. Note that this only works, if we save the posterior samples and use the model fit of `BDgraph` or `bgms`. 

```r
plot_posteriorstructure(res, as.BF = F)
plot_posteriorcomplexity(res, as.BF = F)
plot_structure(res, layoutScale = c(.8,1), palette = "R",
               theme = "TeamFortress", vsize = 6, edge.width = .3, layout = "spring")
```

Furthermore, researcher can wish to aggregate the findings of the network model, commonly done with centrality measures. Due to the discussion around the meaningfulness of centrality measures in psychometric network models, we recommend users to stick to the strength centrality. To obtain the centrality measures, users need to set `save = TRUE` and `centrality = TRUE`, when estimating the network model with `easybgm`. The centrality measures can be inspected with the centrality plot. 

```r
plot_centrality(res, measures = "Strength")
```

## Background Information

For more information on the Bayesian background, its application to networks and the respective plots, check out: 

Huth, K., de Ron, J., Luigjes, J., Goudriaan, A., Mohammadi, R., van Holst, R., Wagenmakers, E.J., \& Marsman, M. (2023). Bayesian Analysis of Cross-sectional Networks: A Tutorial in R and JASP. PsyArXiv https://doi.org/10.31234/osf.io/ub5tc.

## Bug Reports, Feature Request, or Contributing

If you encounter any bugs or have ideas for new features, you can submit them by creating an issue on Github. Additionally, if you want to contribute to the development of `easybgm`, you can initiate a branch with a pull request; we can review and discuss the proposed changes.

## References

Epskamp, S., Cramer, A. O. J., Waldorp, L. J., Schmittmann, V. D., & Borsboom, D. (2012). qgraph: Network Visualizations of Relationships in Psychometric Data. Journal of Statistical Software, 48 . doi: 10.18637/jss.v048.i04

Huth, K., de Ron, J., Luigjes, J., Goudriaan, A., Mohammadi, R., van Holst, R., Wagenmakers, E.J., \& Marsman, M. (2023). Bayesian Analysis of Cross-sectional Networks: A Tutorial in R and JASP. PsyArXiv https://doi.org/10.31234/osf.io/ub5tc.

Marsman, M., Haslbeck, J. M. B. (2023). Bayesian Analysis of the Ordinal Markov Random Field. PsyArXiv https://doi.org/10.31234/osf.io/ukwrf. 

Mohammadi, Reza, and Ernst C Wit. (2015). “BDgraph: An R Package for Bayesian Structure Learning in Graphical Models.” Journal of Statistical Software 89 (3). https://doi.org/10.18637/jss.v089.i03.

Wickham, H. (2016). ggplot2: Elegant graphics for data analysis. Springer-Verlag New York. Retrieved from https://ggplot2.tidyverse.org

Williams, Donald R, and Joris Mulder. (2019). “Bayesian Hypothesis Testing for Gaussian Graphical Models: Conditional Independence and Order Constraints.” PsyArXiv. https://doi.org/10.31234/osf.io/ypxd8.


