\name{DReigen}

\alias{DReigen}
\alias{Rho2EDF}

\title{
Compute exact Demmler-Reinsch eigenvalues
}

\description{
\code{Rho2EDF} maps smoothing parameter value to effective degree of freedom using Demmler-Reinsch eigenvalues
}

\usage{
DReigen(E)

Rho2EDF(values, rho)
}

\arguments{
\item{E}{in the returned list of \code{\link{gps2Grid}}.}
\item{values}{a sequence of Demmler-Reinsch eigenvalues.}
\item{rho}{a grid of log smoothing parameters.}
}

\value{
A sequence of Demmler-Reinsch eigenvalues.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
require(gps)

x <- rnorm(100)
xt <- PlaceKnots(x, d = 4, k = 10)

## set ng = 0 to skip solving penalized least squares
## y = x
setup <- gps2Grid(x, x, xt = xt, d = 4, m = 2, ng = 0)

## compute exact eigenvalues
ei <- DReigen(setup$E)

## compare approximated eigenvalues with exact eigenvalues on the log scale
log.exact <- log(ei$values)
log.approx <- log(setup$eigen$approx.values)
op <- par(mar = c(2, 2, 1.5, 0.5))
plot(log.exact, pch = 19, col = 2, ann = FALSE)
points(log.approx, pch = 19)
legend("topright", legend = c("approx", "exact"), col = 1:2, pch = 19)
title("log(eigenvalues)")
par(op)
}
