\name{DemoBS}

\alias{DemoBS}

\title{
Demonstrate the construction of ordinary B-splines
}

\description{
Demonstrate the construction of 4 ordinary cubic B-splines on 8 knots.
}

\usage{
DemoBS(uniform = TRUE, clamped = FALSE)
}

\arguments{
\item{uniform}{if TRUE, place uniform knots; if FALSE, place non-uniform knots.}
\item{clamped}{if TRUE, place clamped boundary knots when \code{uniform = FALSE}. For aesthetic reason, only boundary knots on the left end are clamped. This parameter is ignored when \code{uniform = TRUE}.}
}

\value{
This function has no returned values.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
require(gps)

## uniform B-splines
DemoBS(uniform = TRUE)

## non-uniform B-splines
DemoBS(uniform = FALSE, clamped = FALSE)

## non-uniform & clamped B-splines
DemoBS(uniform = FALSE, clamped = TRUE)
}
