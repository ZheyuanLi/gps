\name{rspl}

\alias{rspl}
\alias{random.splines}

\title{
Simulate random cubic splines
}

\description{
Simulate random cubic splines.
}

\usage{
rspl(x, domain = NULL, n = 1)
}

\arguments{
\item{x}{\eqn{x}-values where simulated cubic splines are evaluated.}
\item{domain}{a vector of two values giving domain interval \eqn{[a, b]}. Will use \code{min(x)} and \code{max(x)} if not specified.}
\item{n}{number of replicates to simulate.}
}

\value{
A list of four components:
\itemize{
\item \code{y} is a vector of random cubic spline values evaluated at \code{x} when \code{n = 1}, or a matrix of \code{n} columns when \code{n > 1}, where each column is an independent replicate of random cubic splines;
\item \code{b} is a vector of random B-spline coefficients when \code{n = 1}, or a matrix of \code{n} columns when \code{n > 1}, where each column is an independent replicate of random B-spline coefficients;
\item \code{xt} is the full knot sequence for B-splines;
\item \code{domain} gives the domain of the simulated spline(s).
}
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
require(gps)

x <- seq.int(0, 1, 0.01)

## a random cubic spline
y <- rspl(x, n = 1)$y
op <- par(mar = c(2, 2, 1.5, 0.5))
plot(x, y, type = "l", ann = FALSE)
title("a random cubic spline")
par(op)

## 5 random cubic splines
Y <- rspl(x, n = 5)$y
op <- par(mar = c(2, 2, 1.5, 0.5))
matplot(x, Y, type = "l", lty = 1, ylab = "y")
title("5 random cubic splines")
par(op)
}
