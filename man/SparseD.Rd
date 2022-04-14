\name{SparseD}

\alias{SparseD}
\alias{DiffCoef}
\alias{NullD}


\title{
General difference matrices for ordinary B-splines
}

\description{
Compute general difference matrices of order 1 to \eqn{d - 1} for ordinary B-splines of order \eqn{d}. Compute general differences of \eqn{f(x)}'s B-spline coefficients, which are \eqn{f^{(m)}(x)}'s coefficients. Compute general differences of B-spline coefficients. Find null space of a general difference matrix. Compute a basis matrix for the null space of a general difference matrix.
}

\usage{
SparseD(xt, d)

DiffCoef(b, xt, d, m)

NullD(xt, d, m)
}

\arguments{
\item{xt}{full knot sequence for ordinary B-splines (\code{length(xt) >= 2 * d} required).}
\item{d}{B-spline order (\eqn{d \ge 2} required).}
\item{m}{order of the general difference matrix (\eqn{1 \le m \le d - 1} required).}
\item{b}{B-spline coefficients (\code{length(b) == length(xt) - d} required).}
}

\value{
\code{SparseD} returns a list of sparse matrices (of "dgCMatrix" class), giving general difference matrices of order 1 to \eqn{d - 1}.

\code{DiffCoef} returns a vector of coefficients.

\code{NullD} returns a matrix of \eqn{m} columns.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
require(gps)

## example with ordinary B-splines
xt <- c(0, 0, 0, 0, 1, 3, 4, 4, 4, 4)

SparseD(xt, d = 4)

b <- rnorm(6)
DiffCoef(b, xt, d = 4, m = 2)

NullD(xt, d = 4, m = 2)
}
