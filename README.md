# easybgm

The package helps extracting and visualizing the results of a Bayesian analysis of networks commonly found in psychology. Until now it supports the extraction of cross-sectional network models fitted using `BDgraph` and `BGGM`. The package allows to extract the parameter estimates, posterior inclusion probability, inclusion Bayes factor and the posterior density of the parameters. In addition, for `BDgraph` it allows to assess and visualize the posterior probability of the structure space. 

For more information on the Bayesian background, its application to networks and the respective plots, check out: Huth, K., de Ron, J., Luigjes, J., Goudriaan, A., Mohammadi, R., van Holst, R., Wagenmakers, E.J., \& Marsman, M. (2023). Bayesian Analysis of Cross-sectional Networks: A Tutorial in R and JASP. PsyArXiv https://doi.org/10.31234/osf.io/ub5tc.

To download this package use

```r
library("remotes")
install_github("KarolineHuth/easybgm")
```
