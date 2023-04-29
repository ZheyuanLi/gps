\name{PlaceKnots}

\alias{PlaceKnots}

\title{
Automatically place knots according to data
}

\description{
Place knots for ordinary B-splines or periodic B-splines using automatic strategies.
}

\usage{
PlaceKnots(x, d, k, domain = NULL, uniform = FALSE, periodic = FALSE)
}

\arguments{
\item{x}{observed \eqn{x}-values.}
\item{d}{B-spline order.}
\item{k}{number of interior knots.}
\item{domain}{a vector of two values giving domain interval \eqn{[a, b]}. Will use \code{min(x)} and \code{max(x)} if not specified.}
\item{uniform}{if TRUE, place equidistant knots; if FALSE, place quantile knots with clamped boundary knots.}
\item{periodic}{if TRUE, return the domain knot sequence that is sufficient for constructing periodic B-splines (see \code{\link{pbsDesign}}); if FALSE, return the full knot sequence that is required for constructing ordinary B-splines.}
}

\value{
A vector of \eqn{K = k + 2d} knots for ordinary B-splines, or \eqn{k + 2} knots for periodic B-splines.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
require(gps)

x <- rnorm(50)

## uniform knots for uniform cubic B-splines
xt1 <- PlaceKnots(x, d = 4, k = 5, uniform = TRUE)
B1 <- splines::splineDesign(xt1, x, ord = 4)

## clamped quantile knots for clamped non-uniform cubic B-splines
xt2 <- PlaceKnots(x, d = 4, k = 5, uniform = FALSE)
B2 <- splines::splineDesign(xt2, x, ord = 4)

## uniform knots for uniform periodic cubic B-splines
xd1 <- PlaceKnots(x, d = 4, k = 5, uniform = TRUE, periodic = TRUE)
PB1 <- pbsDesign(x, xd1, d = 4)

## quantile knots for non-uniform periodic cubic B-splines
xd2 <- PlaceKnots(x, d = 4, k = 5, uniform = FALSE, periodic = TRUE)
PB2 <- pbsDesign(x, xd2, d = 4)
}
