\name{penalty}

\alias{penalty}
\alias{penalties}
\alias{SparseD}
\alias{PriorCoef}
\alias{MPinv}

\title{
Wiggliness penalties for penalized B-splines
}

\description{
For penalized B-splines (including standard or general P-splines and O-splines), (1) construct matrix \eqn{\bm{D}} in the wiggliness penalty \eqn{\|\bm{D\beta}\|^2}; (2) sample B-spline coefficients from their prior distribution \eqn{\textrm{N}(\bm{0},\ (\bm{D'D})^-)}; (3) compute the Moore-Penrose generalized inverse matrix \eqn{(\bm{D'D})^-}.
}

\usage{
SparseD(xt, d, m = NULL, gps = TRUE)

PriorCoef(n, D)

MPinv(D, only.diag = FALSE)
}

\arguments{
\item{xt}{full knot sequence for ordinary B-splines (\code{length(xt) >= 2 * d}).}
\item{d}{B-spline order (\eqn{d \ge 2}).}
\item{m}{penalty order (\eqn{1 \le m \le d - 1}). Can be a vector of multiple values for \code{SparseD}.}
\item{gps}{if TRUE, return \eqn{\bm{D}_{\textrm{gps}}}; if FALSE, return \eqn{\bm{D}_{\textrm{os}}}.}
\item{n}{number of samples to draw from the prior distribution.}
\item{D}{matrix \eqn{\bm{D}_{\textrm{gps}}} or \eqn{\bm{D}_{\textrm{os}}}.}
\item{only.diag}{if TURE, only diagonal elements are computed.}
}

\details{

\subsection{General Difference Penalty for General P-Splines}{

A general P-spline is characterized by an order-\eqn{m} general difference matrix \eqn{\bm{D}_{\textrm{gps}}}, which can be computed by \code{SparseD(\dots, gps = TRUE)}. For interpretation, the differenced coefficients \eqn{\bm{D}_{\textrm{gps}}\bm{\beta}} are in fact \eqn{f^{(m)}(x)}'s B-spline coefficients, so the penalty is their squared \eqn{L_2} norm.
}

\subsection{Derivative Penalty for O-Splines}{

An O-spline is characterized by \eqn{\bm{D}_{\textrm{os}}} such that \eqn{\|\bm{D}_{\textrm{os}}\bm{\beta}\|^2 = \int_a^b f^{(m)}(x)^2\textrm{d}x}. Since \eqn{f^{(m)}(x)} has B-spline coefficients \eqn{\bm{D}_{\textrm{gps}}\bm{\beta}}, the integral can be shown to be \eqn{\bm{\beta'}\bm{D}_{\textrm{gps}}'\bm{\bar{S}}\bm{D}_{\textrm{gps}}\bm{\beta}}, where \eqn{\bm{\bar{S}}} is the Gram matrix of those B-splines representing \eqn{f^{(m)}(x)}. Following the Cholesky factorization \eqn{\bm{\bar{S}} = \bm{U'U}}, the quadratic form becomes \eqn{\|\bm{U}\bm{D}_{\textrm{gps}}\bm{\beta}\|^2}, so that \eqn{\bm{D}_{\textrm{os}} = \bm{U}\bm{D}_{\textrm{gps}}}. This matrix can be computed by \code{SparseD(\dots, gps = FALSE)}, with \eqn{\bm{\bar{S}}} and \eqn{\bm{D}_{\textrm{gps}}} also returned in a "sandwich" attribute.
}

\subsection{Penalty Matrix}{

We can express the \eqn{L_2} penalty \eqn{\|\bm{D\beta}\|^2} as quadratic form \eqn{\bm{\beta'S\beta}}, where \eqn{\bm{S} = \bm{D'D}} is called a penalty matrix. It is trivial to compute \eqn{\bm{S}} (using function \code{crossprod}) once \eqn{\bm{D}} is available, so we don't feel the need to provide a function for this. Note that the link between \eqn{\bm{D}_{\textrm{os}}} and \eqn{\bm{D}_{\textrm{gps}}} implies a sandwich formula \eqn{\bm{S}_{\textrm{os}} = \bm{D}_{\textrm{gps}}'\bm{\bar{S}}\bm{D}_{\textrm{gps}}}, wherease \eqn{\bm{S}_{\textrm{gps}} = \bm{D}_{\textrm{gps}}'\bm{D}_{\textrm{gps}}}.
}

\subsection{The Bayesian View}{

In the Bayesian view, the penalty \eqn{\bm{\beta'S\beta}} is a Gaussian prior for B-spline coefficients \eqn{\bm{\beta}}. But it is an improper one because \eqn{\bm{S}} has a null space where an unpenalized order-\eqn{m} polynomial lies. Let's decompose \eqn{\bm{\beta} = \bm{\xi} + \bm{\theta}}, where \eqn{\bm{\xi}} (the projection of \eqn{\bm{\beta}} on this null space) is the coefficients of this order-\eqn{m} polynomial, and \eqn{\bm{\theta}} (orthogonal to \eqn{\bm{\xi}}) is the component that can be shrunk to zero by the penalty. As a result, \eqn{\bm{\xi} \propto \bm{1}} is not proper, but \eqn{\bm{\theta} \sim \textrm{N}(\bm{0},\ \bm{S}^-)} is. Function \code{PriorCoef} samples this distribution, and the resulting B-spline coefficients can be used to create random spline curves. The algorithm behind \code{PriorCoef} bypasses the Moore-Penrose generalized inverse and is very efficient. We don't recommend forming this inverse matrix because it, being completely dense, is expensive to compute and store. But if we need it anyway, it can be computed using function \code{MPinv}.
}
}

\value{
\code{SparseD} returns a list of sparse matrices (of "dgCMatrix" class), giving \eqn{\bm{D}_{\textrm{gps}}} or \eqn{\bm{D}_{\textrm{os}}} of order \code{m[1]}, \code{m[2]}, ..., \code{m[length(m)]}. In the latter case, \eqn{\bm{\bar{S}}} (sparse matrices of "dsCMatrix" or "ddiMatrix" class) and \eqn{\bm{D}_{\textrm{gps}}} for computing \eqn{\bm{D}_{\textrm{os}}} are also returned in a "sandwich" attribute.

\code{PriorCoef} returns a list of two components:
\itemize{
\item \code{coef} gives a vector of B-spline coefficients when \code{n = 1}, or a matrix of \code{n} columns when \code{n > 1}, where each column is an independent sample;
\item \code{sigma} is a vector, giving the marginal standard deviation for each B-spline coefficient.
}

\code{MPinv} returns the dense Moore-Penrose generalized inverse matrix \eqn{\bm{(D'D})^-} if \code{only.diag = FALSE}, and the diagonal entries of this matrix if \code{only.diag = TRUE}.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\references{
Zheyuan Li and Jiguo Cao (2022). General P-splines for non-uniform splines, \doi{10.48550/arXiv.2201.06808}
}

\examples{
require(Matrix)
require(gps)

## 11 domain knots at equal quantiles of Beta(3, 3) distribution
xd <- qbeta(seq.int(0, 1, by = 0.1), 3, 3)
## full knots (with clamped boundary knots) for constructing cubic B-splines
xt <- c(0, 0, 0, xd, 1, 1, 1)

## compute D matrices of order 1 to 3 for O-splines
D.os <- SparseD(xt, d = 4, gps = FALSE)
D1.os <- D.os[[1]]; D2.os <- D.os[[2]]; D3.os <- D.os[[3]]

## get D matrices of order 1 to 3 for general P-splines
## we can of course compute them with D.gps <- SparseD(xt, d = 4, gps = TRUE)
## but they are readily stored in the "sandwich" attribute of 'D.os'
D.gps <- attr(D.os, "sandwich")$D
D1.gps <- D.gps[[1]]; D2.gps <- D.gps[[2]]; D3.gps <- D.gps[[3]]

## we can compute the penalty matrix S = D'D
S.gps <- lapply(D.gps, crossprod)
S1.gps <- S.gps[[1]]; S2.gps <- S.gps[[2]]; S3.gps <- S.gps[[3]]
S.os <- lapply(D.os, crossprod)
S1.os <- S.os[[1]]; S2.os <- S.os[[2]]; S3.os <- S.os[[3]]

## if we want to verify the sandwich formula for O-splines
## extract 'Sbar' matrices stored in the "sandwich" attribute
## and compute the relative error between S and t(D) \%*\% Sbar \%*\% D
Sbar <- attr(D.os, "sandwich")$Sbar
Sbar1 <- Sbar[[1]]; Sbar2 <- Sbar[[2]]; Sbar3 <- Sbar[[3]]
range(S1.os - t(D1.gps) \%*\% Sbar1 \%*\% D1.gps) / max(abs(S1.os))
range(S2.os - t(D2.gps) \%*\% Sbar2 \%*\% D2.gps) / max(abs(S2.os))
range(S3.os - t(D3.gps) \%*\% Sbar3 \%*\% D3.gps) / max(abs(S3.os))

## sample B-spline coefficients from their prior distribution
b.gps <- PriorCoef(n = 5, D2.gps)$coef
b.os <- PriorCoef(n = 5, D2.os)$coef
op <- par(mfrow = c(1, 2), mar = c(2, 2, 1.5, 0.5), oma = c(0, 0, 1, 0))
## prior B-spline coefficients with a general difference penalty
matplot(b.gps, type = "l", lty = 1, ann = FALSE)
title("general difference penalty")
## prior B-spline coefficients with a derivative penalty
matplot(b.os, type = "l", lty = 1, ann = FALSE)
title("derivative penalty")
title("random B-spline coefficients from their prior", outer = TRUE)
par(op)

## plot the corresponding cubic splines with these B-spline coefficients
x <- MakeGrid(xd, n = 11)
B <- splines::splineDesign(xt, x, ord = 4, sparse = TRUE)
y.gps <- B \%*\% b.gps
y.os <- B \%*\% b.os
op <- par(mfrow = c(1, 2), mar = c(2, 2, 1.5, 0.5), oma = c(0, 0, 1, 0))
matplot(x, y.gps, type = "l", lty = 1, ann = FALSE)
title("general difference penalty")
matplot(x, y.os, type = "l", lty = 1, ann = FALSE)
title("derivative penalty")
title("random cubic splines with prior B-spline coefficients", outer = TRUE)
par(op)
}
