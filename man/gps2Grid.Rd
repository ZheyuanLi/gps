\name{gps2Grid}

\alias{gps2Grid}

\title{
General P-spline estimation via grid search
}

\description{
Fit a general P-spline to (x, y, w) for an automatically generated grid of log smoothing parameter values.
}

\usage{
gps2Grid(x, y, w = NULL, xt, d = 4, m = 2, deriv.pen = FALSE,
         boundary = "none", ng = 20)
}

\arguments{
\item{x, y, w}{a vector of \eqn{x}-values, \eqn{y}-values and weights.}
\item{xt}{full knot sequence for ordinary B-splines (\code{length(xt) >= 2 * d} required).}
\item{d}{B-spline order (\eqn{d \ge 2} required).}
\item{m}{penalty order (\eqn{1 \le m \le d - 1} required).}
\item{deriv.pen}{TRUE to use derivative penalty; FALSE to use general difference penalty.}
\item{boundary}{type of boundary constraint: "none", "natural" (applicable only to B-splines of even order) or "periodic".}
\item{ng}{number of grid points; can be set to 0 to }
}

\value{
to do
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
require(gps)

## simulate 100 (x, y) data from g(x) = sin(2 * pi * x) on [0, 1]
## x-values are not equidistant but at quantiles of Beta(2, 2)
## note that g(x) is a periodic function
x <- qbeta(seq.int(0, 1, length.out = 100), 2, 2)
gx <- sin(2 * pi * x)
y <- rnorm(length(x), gx, sd = 0.1)

## place quantile knots with clamped boundary knots
xt <- PlaceKnots(x, d = 4, k = 10)

## fit a general P-spline with different boundary constraints
ordinary <- gps2Grid(x, y, xt = xt, d = 4, m = 2)
natural <- gps2Grid(x, y, xt = xt, d = 4, m = 2, boundary = "natural")
periodic <- gps2Grid(x, y, xt = xt, d = 4, m = 2, boundary = "periodic")

## identify the optimal fit minimizing GCV score
opt.ordinary <- which.min(ordinary$pls$gcv)
opt.natural <- which.min(natural$pls$gcv)
opt.periodic <- which.min(periodic$pls$gcv)

## inspect grid search result
## column 1: ordinary cubic spline
## column 2: natural cubic spline
## column 3: periodic cubic spline
op <- par(mfcol = c(2, 3), mar = c(2, 2, 1.5, 0.5))
## ordinary spline
with(ordinary$pls, plot(rho, edf, ann = FALSE))
title("edf v.s. log(lambda)")
with(ordinary$pls, plot(rho, gcv, ann = FALSE))
with(ordinary$pls, points(rho[opt.ordinary], gcv[opt.ordinary], pch = 19))
title("GCV v.s. log(lambda)")
## natural spline
with(natural$pls, plot(rho, edf, ann = FALSE))
title("edf v.s. log(lambda)")
with(natural$pls, plot(rho, gcv, ann = FALSE))
with(natural$pls, points(rho[opt.natural], gcv[opt.natural], pch = 19))
title("GCV v.s. log(lambda)")
## periodic spline
with(periodic$pls, plot(rho, edf, ann = FALSE))
title("edf v.s. log(lambda)")
with(periodic$pls, plot(rho, gcv, ann = FALSE))
with(periodic$pls, points(rho[opt.periodic], gcv[opt.periodic], pch = 19))
title("GCV v.s. log(lambda)")
par(op)

## inspect fitted splines
yhat.ordinary <- with(ordinary, eqn$B \%*\% pls$beta)
yhat.natural <- with(natural, eqn$B \%*\% pls$beta)
yhat.periodic <- with(periodic, eqn$B \%*\% pls$beta)
op <- par(mfcol = c(1, 3), mar = c(2, 2, 1.5, 0.5))
## ordinary spline
matplot(x, yhat.ordinary, type = "l", lty = 1, ann = FALSE)
title("ordinary")
## natural spline
matplot(x, yhat.natural, type = "l", lty = 1, ann = FALSE)
title("natural")
## periodic spline
matplot(x, yhat.periodic, type = "l", lty = 1, ann = FALSE)
title("periodic")
par(op)

## pick and plot the optimal fit minimizing GCV score
best.ordinary <- yhat.ordinary[, opt.ordinary]
best.natural <- yhat.natural[, opt.natural]
best.periodic <- yhat.periodic[, opt.periodic]
op <- par(mfcol = c(1, 3), mar = c(2, 2, 1.5, 0.5))
## ordinary spline
plot(x, y, ann = FALSE)
lines(x, gx, lwd = 2, col = 2)
lines(x, best.ordinary, lwd = 2)
title("ordinary")
## natural spline
plot(x, y, ann = FALSE)
lines(x, gx, lwd = 2, col = 2)
lines(x, best.natural, lwd = 2)
title("natural")
## periodic spline
plot(x, y, ann = FALSE)
lines(x, gx, lwd = 2, col = 2)
lines(x, best.periodic, lwd = 2)
title("periodic")
par(op)
}
