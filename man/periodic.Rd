\name{periodic}

\alias{periodic}
\alias{pbsDesign}
\alias{SparsePD}

\title{
Design matrix and general difference matrices for periodic B-splines
}

\description{
For order-\eqn{d} periodic B-splines, \code{pbsDesign} evaluates B-splines or their derivatives at given \eqn{x}-values, and \code{SparsePD} computes general difference matrices of order 1 to \eqn{d - 1}.
}

\usage{
pbsDesign(x, xd, d, nDeriv = 0, sparse = FALSE, wrap = TRUE)

SparsePD(xd, d, wrap = TRUE)
}

\arguments{
\item{x}{\eqn{x}-values where periodic B-splines are to be evaluated.}
\item{xd}{domain knot sequence for periodic B-splines (\code{length(xd) >= d + 1}).}
\item{d}{B-spline order (\eqn{d \ge 2}).}
\item{nDeriv}{derivative order.}
\item{sparse}{if TRUE, create a sparse design matrix of "dgCMatrix" class.}
\item{wrap}{if TRUE, the knots wrapping strategy is used; if FALSE, the linear constraint strategy is used.}
}

\details{
These functions perform type-2 construction, by transforming design matrix and general difference matrices for ordinary B-splines to satisfy periodic boundary constraints (see Details). By contrast, \code{pbsDesign} and \code{SparsePD} in \pkg{gps} perform type-1 construction by basis wrapping.

A spline \eqn{f(x)} on domain \eqn{[a, b]} can be constructed to satisfy periodic boundary constraints, that is, \eqn{f^{(q)}(a) = f^{(q)}(b)}, \eqn{q} = 0, 1, ..., degree - 1. These are actually linear equality constraints

Unlike ordinary B-splines, period B-splines do not require explicit auxiliary boundary knots for their construction. The magic is that auxiliary boundary knots will be automatically positioned by periodic extension of interior knots.

Denote the domain knot sequence by \eqn{s_0, s_1, s_2, \ldots, s_k, s_{k + 1}}, where \eqn{(s_j)_1^k} are interior knots and \eqn{s_0 = a}, \eqn{s_{k + 1} = b} are domain endpoints. For order-\eqn{d} B-splines, we replicate the first \eqn{d - 1} interior knots (after adding \eqn{b - a}) to the right of \eqn{[a, b]} for an augmented set of \eqn{K = k + d + 1} knots, which spawns \eqn{p = K - d = k + 1} ordinary B-splines. It turns out that periodic B-splines can be obtained by wrapping segments of those ordinary B-splines that stretch beyond \eqn{[a, b]} to the start of the domain (a demo is offered by \code{\link{DemoPBS}}).

Note that we must have at least \eqn{d - 1} interior knots to do such periodic extension. This means that \eqn{d + 1} domain knots are required at a minimum for construction of periodic B-splines.
}

\value{
\code{pbsDesign} returns a design matrix with \code{length(x)} rows and \code{length(xd) - 1} columns. \code{SparsePD} returns a list of sparse matrices (of "dgCMatrix" class), giving general difference matrices of order 1 to \eqn{d - 1}.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
require(gps)

## 5 domain knots: three interior knots 0.5, 1.5 and 1.8 in domain [0, 3]
xd <- c(0, 0.5, 1.5, 1.8, 3)

## make a grid
x <- MakeGrid(xd, n = 10)

## construct periodic cubic B-splines
PB1 <- pbsDesign(x, xd, d = 4, wrap = TRUE)
PB2 <- pbsDesign(x, xd, d = 4, wrap = FALSE)

## construct general difference matrices of order 1 to 3
SparsePD(xd, d = 4, wrap = TRUE)
SparsePD(xd, d = 4, wrap = FALSE)
}
