# easybgm: Easy Bayesian Graphical Modeling

The `R` package `easybgm` provides a user-friendly pacakge to conduct a Bayesian analysis of networks. In particular, it helps fitting, extracting and visualizing the results of a Bayesian graphical model commonly found in psychology. Until now, the package supports fitting and extracting results of cross-sectional network models using `BDgraph`, `BGGM`, `bgms`. The package extracts the parameter estimates, posterior inclusion probability, inclusion Bayes factor and the posterior density of the parameters.

## Installation

To download this package use

```r
install.packages("remotes")
library("remotes")
install_github("KarolineHuth/easybgm")
```

## Overview

### Model estimation

### Result extraction

### Visualization

## Example

We want to illustrate the package use with an example. In particular, we use the Wenchuan data which can be loaded with 
the package `bgms`. We fit the model and extract its results with `bgm_fit`

```r
library(easybgm)
library(bgms)

data <- Wenchuan
res <- bgm_fit(data, type = "ordinal")
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

Lastly, we can assess the structure specifically with three plots. Note that this only works, if we save the posterior samples and use the model fit of `BDgraph` or `bgms`. 

```r
plot_posteriorstructure(res, as.BF = F)
plot_posteriorcomplexity(res, as.BF = F)
plot_structure(res, layoutScale = c(.8,1), palette = "R",
               theme = "TeamFortress", vsize = 6, edge.width = .3, layout = "spring")
```


## Background Information

For more information on the Bayesian background, its application to networks and the respective plots, check out: Huth, K., de Ron, J., Luigjes, J., Goudriaan, A., Mohammadi, R., van Holst, R., Wagenmakers, E.J., \& Marsman, M. (2023). Bayesian Analysis of Cross-sectional Networks: A Tutorial in R and JASP. PsyArXiv https://doi.org/10.31234/osf.io/ub5tc.
