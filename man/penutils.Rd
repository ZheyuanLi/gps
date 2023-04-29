\name{penutils}

\alias{penutils}
\alias{penalty.utils}
\alias{penalty.utilities}
\alias{DiffCoef}
\alias{btSb}

\title{
Utility functions for working with wiggliness penalties
}

\description{
Evaluate \eqn{\|\bm{D\beta}\|^2} without using \eqn{\bm{D}}.
}

\usage{
DiffCoef(b, xt, d, m)

btSb(b, xt, d, m)
}

\arguments{
\item{b}{a vector of B-spline coefficients (\code{length(b) == length(xt) - d}).}
\item{xt}{full knot sequence for ordinary B-splines (\code{length(xt) >= 2 * d}).}
\item{d}{B-spline order (\eqn{d \ge 2}).}
\item{m}{penalty order (\eqn{1 \le m \le d - 1}).}
}

\details{

\subsection{Implicit Evaluation of the Penalty}{

Sometimes we want to evaluate the penalty \eqn{\|\bm{D\beta}\|^2} for some \eqn{\bm{\beta}}. The obvious way is to do the matrix-vector multiplication \eqn{\bm{D\beta}} then compute its \eqn{L_2} norm, however, implicit evaluation without using \eqn{\bm{D}} is possible. For general P-splines, we can calculate \eqn{\bm{D}_{\textrm{gps}}\bm{\beta}} by taking order-\eqn{m} general differences between elements of \eqn{\bm{\beta}}, and function \code{DiffCoef} does this. For O-splines, the evaluation can be more refined. Denote domain knots by \eqn{s_0,\ s_1,\ s_2,\ \ldots,\ s_k,\ s_{k + 1}}, where \eqn{(s_j)_1^k} are interior knots and \eqn{s_0 = a}, \eqn{s_{k + 1} = b} are domain endpoints. The derivative penalty adds up local wiggliness measure on each interval: \eqn{\int_a^b f^{(m)}(x)^2\mathrm{d}x = \sum_{j = 0}^k\int_{s_j}^{s_{j + 1}} f^{(m)}(x)^2\mathrm{d}x}. Function \code{btSb} calculates each integral in the summation and returns those additive components in a vector.
}
}

\value{
\code{DiffCoef} (for general P-splines only) returns \eqn{\bm{D}_{\textrm{gps}}\bm{\beta}} as a vector.

\code{btSb} (for O-splines only) returns a vector with element \eqn{\int_{s_j}^{s_{j + 1}} f^{(m)}(x)^2\mathrm{d}x}.
}

\examples{
require(Matrix)
require(gps)

## 11 domain knots at equal quantiles of Beta(3, 3) distribution
xd <- qbeta(seq.int(0, 1, by = 0.1), 3, 3)
## full knots (with clamped boundary knots) for constructing cubic B-splines
xt <- c(0, 0, 0, xd, 1, 1, 1)

## compute 2nd order D matrix for O-splines
D.os <- SparseD(xt, d = 4, m = 2, gps = FALSE)
D2.os <- D.os$order.2

## get 2nd order D matrix for general P-splines
## we can of course compute it with D.gps <- SparseD(xt, d = 4, m = 2, gps = TRUE)
## but it is readily stored in the "sandwich" attribute of 'D.os'
D.gps <- attr(D.os, "sandwich")$D
D2.gps <- D.gps$order.2

## random B-spline coefficients
b <- rnorm(ncol(D2.gps))

## two ways to evaluate a difference penalty
diff.b1 <- DiffCoef(b, xt, d = 4, m = 2)  ## implicit
diff.b2 <- as.numeric(D2.gps \%*\% b)       ## explicit
range(diff.b1 - diff.b2) / max(abs(diff.b1))

## several ways to evaluate a derivative penalty
sum(btSb(b, xt, d = 4, m = 2))  ## recommended
sum(as.numeric(D2.os \%*\% b) ^ 2)
S2.os <- crossprod(D2.os); sum(b * as.numeric(S2.os \%*\% b))
}
