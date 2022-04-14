\name{DemoBS}

\alias{DemoBS}

\title{
Demonstrate construction of ordinary B-splines
}

\description{
Demonstrate construction of 4 ordinary cubic B-splines on 8 knots.
}

\usage{
DemoBS(uniform = TRUE, clamped = FALSE)
}

\arguments{
\item{uniform}{TRUE for uniform knots and FALSE for non-uniform knots.}
\item{clamped}{TRUE for clamped boundary knots (for aesthetic reason only boundary knots on the left end are clamped).}
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
