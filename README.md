# gps (version 1.1)

**General P-Splines**

General P-splines are non-uniform B-splines penalized by a general difference penalty, proposed by [Li and Cao (2022)](https://arxiv.org/abs/2201.06808). Constructible on arbitrary knots, they extend the standard P-splines of [Eilers and Marx (1996)](https://doi.org/10.1214/ss/1038425655). They are also related to the O-splines of [O'Sullivan (1986)](https://doi.org/10.1214/ss/1177013525) via a sandwich formula that links a general difference penalty to a derivative penalty. The package includes routines for setting up and handling difference and derivative penalties. It also fits P-splines and O-splines to (x, y) data (optionally weighted) for a grid of smoothing parameter values in the automatic search intervals of [Li and Cao (2023)](https://doi.org/10.1007/s11222-022-10178-z). It aims to facilitate other packages to implement P-splines or O-splines as a smoothing tool in their model estimation framework.

# Installation from CRAN

Version 1.1 will appear [on CRAN](https://CRAN.R-project.org/package=gps) sometime in May, 2023.

```r
install.packages("gps")
```

# Installation from GitHub

```r
## you may need to first install package 'remotes' from CRAN
remotes::install_github("ZheyuanLi/gps")
```

# Vignette

R code for [Li and Cao (2022)](https://arxiv.org/abs/2201.06808) is at: https://github.com/ZheyuanLi/gps-vignettes/blob/main/gps1.pdf

R code for [Li and Cao (2023)](https://doi.org/10.1007/s11222-022-10178-z) is at: https://github.com/ZheyuanLi/gps-vignettes/blob/main/gps2.pdf
