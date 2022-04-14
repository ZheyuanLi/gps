\name{DemoPBS}

\alias{DemoPBS}

\title{
Demonstrate construction of periodic B-splines
}

\description{
Demonstrate construction of 6 periodic cubic B-splines on 7 domain knots.
}

\usage{
DemoPBS(uniform = TRUE)
}

\arguments{
\item{uniform}{TRUE for uniform knots and FALSE for non-uniform knots.}
}

\value{
This function has no returned values.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
require(gps)

## uniform periodic cubic B-splines
DemoPBS(uniform = TRUE)

## non-uniform periodic cubic B-splines
DemoPBS(uniform = FALSE)
}
