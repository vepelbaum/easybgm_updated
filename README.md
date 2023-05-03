# easybgm: Easy Bayesian Graphical Modeling

The `R` package `easybgm` provides a user-friendly package to conduct a Bayesian analysis of psychometric networks. In particular, it helps fitting, extracting and visualizing the results of a Bayesian graphical model commonly found in social-behavioral science. The package is a wrapper around the implemented Until now, the package supports fitting and extracting results of cross-sectional network models using `BDgraph`, `BGGM`, `bgms`. The package extracts the parameter estimates, posterior inclusion probability, inclusion Bayes factor and the posterior density of the parameters.

## Installation

To install this package from Github use

```r
install.packages("remotes")
library("remotes")
install_github("KarolineHuth/easybgm")
```

## Overview

### Estimation

### Visualization

## Example

We want to illustrate the package use with an example. In particular, we use the Wenchuan data which can be loaded with 
the package `bgms`. We fit the model and extract its results with `bgm_fit`

```r
library(easybgm)
library(bgms)

data <- Wenchuan
res <- easybgm(data, type = "ordinal")
```

Having fitted the model, we can now visualize its results. In a first step, we assess the edge evidence plot. 

```r
plot_edgeevidence(res, edge.width = 2, split = T)
```

Furthermore, we can look at the network plot. 

```r
plot_network(res, layout = "spring", 
             layoutScale = c(.8,1), palette = "R",
             theme = "TeamFortress", vsize = 6)
```

In addition we can obtain posterior samples from the posterior distribution by setting `save = TRUE` and thereby open up new possibilities of assessing the model. We can extract the posterior density of the parameters with a parameter forest plot. 

```r
plot_parameterHDI(fit)
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
plot_centrality(fit, measures = "Strength")
```

## Background Information

For more information on the Bayesian background, its application to networks and the respective plots, check out: 

Huth, K., de Ron, J., Luigjes, J., Goudriaan, A., Mohammadi, R., van Holst, R., Wagenmakers, E.J., \& Marsman, M. (2023). Bayesian Analysis of Cross-sectional Networks: A Tutorial in R and JASP. PsyArXiv https://doi.org/10.31234/osf.io/ub5tc.

## Bug Reports, Feature Request, or Contributing



## References

Huth, K., de Ron, J., Luigjes, J., Goudriaan, A., Mohammadi, R., van Holst, R., Wagenmakers, E.J., \& Marsman, M. (2023). Bayesian Analysis of Cross-sectional Networks: A Tutorial in R and JASP. PsyArXiv https://doi.org/10.31234/osf.io/ub5tc.

Mohammadi, Reza, and Ernst C Wit. 2015. “BDgraph: An R Package for Bayesian Structure Learning in Graphical Models.” Journal of Statistical Software 89 (3). https://doi.org/10.18637/jss.v089.i03.

Williams, Donald R, and Joris Mulder. 2019. “Bayesian Hypothesis Testing for Gaussian Graphical Models: Conditional Independence and Order Constraints.” PsyArXiv. https://doi.org/10.31234/osf.io/ypxd8.


