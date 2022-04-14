\name{Ospline}

\alias{Ospline}
\alias{SandBar}
\alias{SparseS}
\alias{btSb}

\title{
Build derivative penalty matrix or its sparse "root"
}

\description{
Compute order-\eqn{m} derivative penalty matrix \eqn{\bm{S}} or its sparse "root" \eqn{\bm{K}} such that \eqn{\bm{K'K} = \bm{S}}.

Evaluate derivative penalty 

Evaluate derivative penalty \eqn{\int_a^b f^{(m)}(x)^2\mathrm{d}x = \bm{\beta'S\beta}}, where \eqn{\bm{S}} is the order-\eqn{m} derivative penalty matrix and \eqn{\bm{\beta}} is the vector of B-spline coefficients for \eqn{f(x)}.
}

\usage{
SandBar(xt, d, m)

SparseS(xt, d, m, root = FALSE)

btSb(b, xt, d, m)
}

\arguments{
\item{xt}{full knot sequence for ordinary B-splines (\code{length(xt) >= 2 * d} required).}
\item{d}{B-spline order (\eqn{d \ge 2} required).}
\item{m}{order of the derivative penalty (\eqn{0 \le m \le d - 1} required).}
\item{root}{if TRUE, return the sparse "root" of the derivative penalty matrix.}
\item{b}{B-spline coefficients (\code{length(b) == length(xt) - d} required).}
}

\details{
Build \eqn{\bm{\bar{S}}} in the sandwich formula \eqn{\bm{S} = \bm{D'\bar{S}D}}.

Given that \eqn{\bm{S}} can be obtained from \code{\link{SparseS}}, it is straightforward to compute the quadratic form \eqn{\bm{\beta'S\beta}}. However, the function takes an alternative approach.

Denote the domain knot sequence by \eqn{s_0, s_1, s_2, \ldots, s_k, s_{k + 1}}, where \eqn{(s_j)_1^k} are interior knots and \eqn{s_0 = a}, \eqn{s_{k + 1} = b} are domain endpoints. The derivative penalty aggregates roughness penalties on all knot spans: \deqn{\int_a^b f^{(m)}(x)^2\mathrm{d}x = \sum_{j = 0}^k\int_{s_j}^{s_{j + 1}} f^{(m)}(x)^2\mathrm{d}x.} The function calculates and returns those additive components in a vector.
}

\value{
\code{SandBar} returns a sparse matrice of "dsCMatrix" (when \eqn{m < d - 1}) or "ddiMatrix" (when \eqn{m = d - 1}) class.

\code{SparseS} returns a sparse matrice giving the order-\eqn{m} derivative penalty matrix (of "dsCMatrix" class) or its sparse "root" (of "dgCMatrix" class).

\code{btSb} returns a vector of evaluated derivative penalties at all knot spans.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
require(gps)

## example with ordinary B-splines
xt <- c(0, 0, 0, 0, 1, 3, 4, 4, 4, 4)

SandBar(xt, d = 4, m = 2)

SparseS(xt, d = 4, m = 2)

b <- rnorm(6)
btSb(b, xt, d = 4, m = 2)
}
