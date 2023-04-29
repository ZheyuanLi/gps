\name{DemoNull}

\alias{DemoNull}

\title{
Demonstrate the null space of P-splines
}

\description{
Cubic P-splines set up with non-uniform B-splines and a 2nd order standard or general difference penalty are fitted to observations simulated from \eqn{y = x}. Should the resulting standard or general P-splines have the correct null space, the limiting fit at \eqn{\lambda = +\infty} will be a straight line regardless of knot locations. In this demo, non-uniform knots from different distributions (primarily Beta distributions with varying shape parameters) are attempted. Results show that standard P-splines have an incorrect and unpredictable limiting behavior that is sensitive to knot locations, whereas general P-splines have a correct and consistent limiting behavior.
}

\usage{
DemoNull(n, k, gps = FALSE)
}

\arguments{
\item{n}{number of simulated observations from \eqn{y = x}.}
\item{k}{number of interior knots to place.}
\item{gps}{if TRUE, fit general P-splines; if FALSE, fit standard P-splines.}
}

\value{
This function has no returned values.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
require(gps)

## standard P-splines
DemoNull(n = 100, k = 10, gps = FALSE)

## general P-splines
DemoNull(n = 100, k = 10, gps = TRUE)
}
