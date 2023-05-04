# easybgm: Easy Bayesian Graphical Modeling

DISCLAIMER: The package is still undergoing rapid changes and changes to function or arguments names may occur. 

The `R` package `easybgm` provides a user-friendly package for performing a Bayesian analysis of psychometric networks. In particular, it helps to fit, extract, and visualize the results of a Bayesian graphical model commonly used in the social-behavioral sciences. The package is a wrapper around existing packages. So far, the package supports fitting and extracting results of cross-sectional network models using `BDgraph` (Mohammadi \& Wit, 2015), `BGGM`(Williams \& Mulder, 2019), `bgms` (Marsman \& Haslbeck, 2023). As output, the package extracts the parameter estimates, the posterior inclusion probability, the inclusion Bayes factor, and optionally posterior samples of the parameters and the centrality. The package furthermore provides an extensive suite of visualization functions. 

## Installation

To install this package from Github use

```r
install.packages("remotes")
remotes::install_github("KarolineHuth/easybgm")
```

To rather install the most up-to-date version in development use 

```r
install.packages("remotes")
remotes::install_github("KarolineHuth/easybgm", ref = "developer")
```

## Overview

### Estimation

Note that for the data type ordinal and binary, easybgm by default uses `bgms`, for mixed and continuous data it uses by default `BDgraph`. Users can specify the preferred package, with the `package` argument. Until now, easybgm supports `BDgraph`, `BGGM`, and `bgms`. 

### Visualization


#### Edge Evidence Plot 

In the edge evidence plot, edges represent the inclusion Bayes factor $\text{BF}_{10}$. The edge evidence plot aids researchers in deciding which edges provide robust inferential conclusions: red edges indicate evidence for edge absence (i.e., conditional independence), grey edges indicate the absence of evidence, and blue edges indicate evidence for edge presence (i.e., conditional dependence).

#### Network Plot

#### Structure Plots

#### 95 \% Highest Density Intervals of the Posterior Parameter Distribution

#### Centrality Plot

## Example

We want to illustrate the package use with an example. In particular, we use the Wenchuan data which can be loaded with 
the package `bgms`. We fit the model and extract its results with the function `easybgm`. We specify the data and the data type, which in this case is `ordinal`. 

```r
library(easybgm)
library(bgms)

data <- Wenchuan
res <- easybgm(data, type = "ordinal")
```

Having fitted the model, we can now visualize its results. In a first step, we assess the edge evidence plot in which edges represent the inclusion Bayes factor $\text{BF}_{10}$. The edge evidence plot aids researchers in deciding which edges provide robust inferential conclusions: red edges indicate evidence for edge absence (i.e., conditional independence), grey edges indicate the absence of evidence, and blue edges indicate evidence for edge presence (i.e., conditional dependence). Especially in a large network, it is recommended to split the edge evidence plot in two parts by setting the `split` argument to `TRUE`. As such, the left plot shows present edges (i.e., $\text{BF}_{10} > 1$), where blue edges represent evidence for inclusion ($\text{BF}_{10} > 10$) and grey edges absence of evidence ($1 < \text{BF}_{10} < 10$). The right edge evidence plot depicting absent edges (i.e., $\text{BF}_{10} < 1$) with evidence for exclusion shown as red ($\text{BF}_{01} > 10$) and inconclusive evidence as grey ($0.1 < \text{BF}_{10} < 1$).  

```r
plot_edgeevidence(res, edge.width = 2, split = T)
```

Furthermore, we can look at the network plot in which edges indicate the strength of partial association between two nodes. The network plot shows all edges with an inclusion Bayes factor larger than $1$. Edge thickness and saturation represent the strength of the association; the thicker the edge, the stronger the association. Red edges indicate negative relations and blue edges indicate positive associations.

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

Huth, K., de Ron, J., Luigjes, J., Goudriaan, A., Mohammadi, R., van Holst, R., Wagenmakers, E.J., \& Marsman, M. (2023). Bayesian Analysis of Cross-sectional Networks: A Tutorial in R and JASP. PsyArXiv https://doi.org/10.31234/osf.io/ub5tc.

Marsman, M., Haslbeck, J. M. B. (2023). Bayesian Analysis of the Ordinal Markov Random Field. PsyArXiv https://doi.org/10.31234/osf.io/ukwrf. 

Mohammadi, Reza, and Ernst C Wit. (2015). “BDgraph: An R Package for Bayesian Structure Learning in Graphical Models.” Journal of Statistical Software 89 (3). https://doi.org/10.18637/jss.v089.i03.

Williams, Donald R, and Joris Mulder. (2019). “Bayesian Hypothesis Testing for Gaussian Graphical Models: Conditional Independence and Order Constraints.” PsyArXiv. https://doi.org/10.31234/osf.io/ypxd8.


