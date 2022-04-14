# gps (version 1.1, under development)

**General P-Splines**

The package contain routines to construct general P-splines for non-uniform B-splines on arbitrary knots, proposed by [Li and Cao (2022)](https://arxiv.org/abs/2201.06808). This P-spline variant extends the standard P-splines of [Eilers and Marx (1996)](https://doi.org/10.1214/ss/1038425655) that are tailored for uniform B-splines on equidistant knots. The package includes `SparseD` for computing general difference matrices, `SparseS` for computing derivative penalty matrix or its sparse root using a sandwich formula, as well as several demos on B-splines and P-splines. It aims to facilitate other packages to implement general P-splines as a smoothing tool in their model estimation framework. See for example, [**gps.mgcv**](https://github.com/ZheyuanLi/gps.mgcv).

# Installation from CRAN

Version 1.0 is [on CRAN](https://CRAN.R-project.org/package=gps). However, **gps.mgcv** requires version 1.1, so I recommend installation from this GitHub repository instead.

Version 1.1 is under active development, and we will release both **gps_1.1** and **gps.mgcv_1.0** on CRAN when they are ready.

# Installation from GitHub

```r
## you may need to first install package 'devtools' from CRAN
devtools::install_github("ZheyuanLi/gps")
```

# Vignette

R code for [Li and Cao (2022)](https://arxiv.org/abs/2201.06808) is at: https://github.com/ZheyuanLi/gps-vignettes/blob/main/gps1.pdf
