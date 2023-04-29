\name{GramBS}

\alias{GramBS}

\title{
Gram matrix of B-splines
}

\description{
Compute the Gram matrix \eqn{\bm{G}}, i.e., the matrix of inner products between B-splines \eqn{b_1(x), b_2(x), \ldots}. Precisely, its element is \eqn{G_{uv} = \int b_u(x)b_v(x)\textrm{d}x}. Such matrix is useful for estimating functional linear models.

The Gram matrix of differentiated B-splines gives the derivative penalty matrix \eqn{\bm{S}} for O-splines. Precisely, its element is \eqn{S_{uv} = \int b_u^{(m)}(x)b_v^{(m)}(x)\textrm{d}x}. Such matrix is straightforward to compute using the results of \code{\link{SparseD}}; see Examples.
}

\usage{
GramBS(xt, d)
}

\arguments{
\item{xt}{full knot sequence for ordinary B-splines (\code{length(xt) >= 2 * d}).}
\item{d}{B-spline order (\eqn{d \ge 2}).}
}

\value{
A sparse matrix of "dsCMatrix" class.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
require(gps)
require(Matrix)

## 11 domain knots at equal quantiles of Beta(3, 3) distribution
xd <- qbeta(seq.int(0, 1, by = 0.1), 3, 3)
## full knots (with clamped boundary knots) for constructing cubic B-splines
xt <- c(0, 0, 0, xd, 1, 1, 1)
## compute Gram matrix of B-splines
G <- GramBS(xt, d = 4)
round(G, digits = 3)

## Gram matrix of differentiated B-splines, i.e., a derivative penalty matrix
## compute derivative penalty matrices of all orders (m = NULL in SparseD)
D <- SparseD(xt, d = 4, gps = FALSE)
S <- lapply(D, crossprod)
lapply(S, round, digits = 1)
}
