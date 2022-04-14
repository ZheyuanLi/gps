\name{prior}

\alias{prior}
\alias{PriorCoef}
\alias{MPinverse}

\title{
Sample B-spline coefficients from their prior distribution
}

\description{
In the Bayesian view, the \eqn{L_2} penalty \eqn{\|\bm{D\beta}\|^2} for B-spline coefficients \eqn{\bm{\beta}} is an improper Gaussian prior: \eqn{\bm{\beta} \sim \textrm{N}(\bm{0},\ (\bm{D'D})^-)}, where the inverse denotes the Moore-Penrose generalized inverse.
}

\usage{
PriorCoef(n, D)

MPinverse(D, diag.only = FALSE)
}

\arguments{
\item{n}{number of samples to draw.}
\item{D}{the matrix}
\item{diag.only}{if TRUE, only diagonal elements}
}

\value{
A vector of coefficients when \code{n = 1}, or a coefficient matrix with \code{n} columns when \code{n > 1}.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
require(gps)

## 11 domain knots at equal quantiles of Beta(3, 3) distribution
xd <- qbeta(seq.int(0, 1, by = 0.1), 3, 3)
## full knots (with clamped boundary knots) for constructing cubic B-splines
xt <- c(0, 0, 0, xd, 1, 1, 1)

## 2nd order general difference matrix
D <- SparseD(xt, d = 4)[[2]]
## sparse "root" of the 2nd order derivative penalty matrix
K <- SparseS(xt, d = 4, m = 2, root = TRUE)

## sample B-spline coefficients from their prior distribution
b1 <- PriorCoef(n = 5, D)
b2 <- PriorCoef(n = 5, K)
op <- par(mfrow = c(1, 2), mar = c(2, 2, 1.5, 0.5), oma = c(0, 0, 1, 0))
## prior B-spline coefficients for general difference penalty
matplot(b1, type = "l", lty = 1, ann = FALSE)
title("general difference penalty")
## prior B-spline coefficients for derivative penalty
matplot(b2, type = "l", lty = 1, ann = FALSE)
title("derivative penalty")
title("random B-spline coefficients from their prior", outer = TRUE)
par(op)

## plot the corresponding cubic splines for these B-spline coefficients
x <- MakeGrid(xd, n = 11)
B <- splines::splineDesign(xt, x, ord = 4, sparse = TRUE)
y1 <- B \%*\% b1
y2 <- B \%*\% b2
op <- par(mfrow = c(1, 2), mar = c(2, 2, 1.5, 0.5), oma = c(0, 0, 1, 0))
matplot(x, y1, type = "l", lty = 1, ann = FALSE)
title("general difference penalty")
matplot(x, y2, type = "l", lty = 1, ann = FALSE)
title("derivative penalty")
title("random cubic splines with prior B-spline coefficients", outer = TRUE)
par(op)
}
