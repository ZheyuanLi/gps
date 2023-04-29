\name{gps2GS}

\alias{gps2GS}
\alias{DemoRhoLim}

\title{
Penalized B-splines estimation with automatic grid search of their smoothing parameter
}

\description{
Fit penalized B-splines (including standard or general P-splines and O-splines) to \code{(x, y, w)} for a grid of smoothing parameter values in the automatic search intervals of Li and Cao (2023). The GCV score and effective degree of freedom of each fit are also returned.
}

\usage{
gps2GS(x, y, w = NULL, xt, d = 4, m = 2, gps = TRUE, periodic = FALSE,
       ng = 20, scalePen = TRUE)
       
DemoRhoLim(fit, plot = TRUE)
}

\arguments{
\item{x, y, w}{a vector of \eqn{x}-values, \eqn{y}-values and weights.}
\item{xt}{full knot sequence for ordinary B-splines (\code{length(xt) >= 2 * d}).}
\item{d}{B-spline order (\eqn{d \ge 2}).}
\item{m}{penalty order (\eqn{1 \le m \le d - 1}).}
\item{gps}{if TRUE, use a difference penalty; if FALSE, use a derivative penalty.}
\item{periodic}{if TRUE, periodic boundary conditions are applied to B-splines and their penalty, so that periodic P-splines are estimated.}
\item{ng}{number of grid points in the grid search of \eqn{\rho}; can be set to 0 to set up the grid search only, without actual P-splines estimation.}
\item{scalePen}{if TRUE, scale the penalty matrix \eqn{\bm{S}} (as \pkg{mgcv} does).}
\item{fit}{fitted P-splines returned by \code{gps2GS}.}
\item{plot}{if TRUE, produce summary plots.}
}

\details{
We smooth \eqn{y_i} using \eqn{f(x_i) = \bm{B_i\beta}}, where \eqn{\bm{B_i}} is \eqn{i}-th row of the B-spline design matrix \eqn{\bm{B}} and \eqn{\bm{\beta}} is a vector of B-spline coefficients. These coefficients are estimated by minimizing: \deqn{\|\bm{y} - \bm{B\beta}\|^2 + \exp(\rho)\cdot\|\bm{D\beta}\|^2,} where the \eqn{L_2} penalty \eqn{\|\bm{D\beta}\|^2} is some wiggliness measure for \eqn{f(x)} and \eqn{\rho \in (-\infty, +\infty)} is a smoothing parameter.

}

\value{
\code{gps2GS} returns a large list with the following components:

\itemize{
\item \code{eqn}
\item \code{eigen}
\item \code{rho.lim}
\item \code{E}
\item \code{pwls}
}
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\references{
Zheyuan Li and Jiguo Cao (2023). Automatic search intervals for the smoothing parameter in penalized splines, \emph{Statistics and Computing}, \doi{10.1007/s11222-022-10178-z}
}

\examples{
require(gps)

x <- rnorm(100)
xt <- PlaceKnots(x, d = 4, k = 10)

## set ng = 0 to set up grid search only
## here the y-values does not matter; we simply use the x-values
setup <- gps2GS(x, x, xt = xt, d = 4, m = 2, ng = 0)

## compute exact eigenvalues
DemoResult <- DemoRhoLim(setup)

## simulate 100 (x, y) data from g(x) = sin(2 * pi * x) on [0, 1]
## x-values are not equidistant but at quantiles of Beta(2, 2)
## note that g(x) is a periodic function
x <- qbeta(seq.int(0, 1, length.out = 100), 2, 2)
gx <- sin(2 * pi * x)
y <- rnorm(length(x), gx, sd = 0.1)

## place quantile knots with clamped boundary knots
xt <- PlaceKnots(x, d = 4, k = 10)

## fit a general P-spline with different boundary constraints
ordinary <- gps2GS(x, y, xt = xt, d = 4, m = 2)
periodic <- gps2GS(x, y, xt = xt, d = 4, m = 2, periodic = TRUE)

## identify the optimal fit minimizing GCV score
opt.ordinary <- which.min(ordinary$pwls$gcv)
opt.periodic <- which.min(periodic$pwls$gcv)

## inspect grid search result
## column 1: ordinary cubic spline
## column 2: periodic cubic spline
op <- par(mfcol = c(2, 2), mar = c(2, 2, 1.5, 0.5))
## ordinary spline
with(ordinary$pwls, plot(rho, edf, ann = FALSE))
title("edf v.s. log(lambda)")
with(ordinary$pwls, plot(rho, gcv, ann = FALSE))
with(ordinary$pwls, points(rho[opt.ordinary], gcv[opt.ordinary], pch = 19))
title("GCV v.s. log(lambda)")
## periodic spline
with(periodic$pwls, plot(rho, edf, ann = FALSE))
title("edf v.s. log(lambda)")
with(periodic$pwls, plot(rho, gcv, ann = FALSE))
with(periodic$pwls, points(rho[opt.periodic], gcv[opt.periodic], pch = 19))
title("GCV v.s. log(lambda)")
par(op)

## inspect fitted splines
yhat.ordinary <- with(ordinary, eqn$B \%*\% pwls$beta)
yhat.periodic <- with(periodic, eqn$B \%*\% pwls$beta)
op <- par(mfrow = c(1, 2), mar = c(2, 2, 1.5, 0.5))
## ordinary spline
matplot(x, yhat.ordinary, type = "l", lty = 1, ann = FALSE)
title("ordinary")
## periodic spline
matplot(x, yhat.periodic, type = "l", lty = 1, ann = FALSE)
title("periodic")
par(op)

## pick and plot the optimal fit minimizing GCV score
best.ordinary <- yhat.ordinary[, opt.ordinary]
best.periodic <- yhat.periodic[, opt.periodic]
op <- par(mfrow = c(1, 2), mar = c(2, 2, 1.5, 0.5))
## ordinary spline
plot(x, y, ann = FALSE)
lines(x, gx, lwd = 2, col = 2)
lines(x, best.ordinary, lwd = 2)
title("ordinary")
## periodic spline
plot(x, y, ann = FALSE)
lines(x, gx, lwd = 2, col = 2)
lines(x, best.periodic, lwd = 2)
title("periodic")
par(op)
}
