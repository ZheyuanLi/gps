\name{DemoSpl}

\alias{DemoSpl}

\title{
Demonstrate a polynomial spline and its B-spline representation
}

\description{
Demonstrate a cubic spline and its B-spline representation.
}

\usage{
DemoSpl(uniform = TRUE)
}

\arguments{
\item{uniform}{if TRUE, place uniform knots; if FALSE, place non-uniform knots.}
}

\value{
A list giving domain knots, B-spline coefficients and piecewise polynomial coefficients of the illustrated cubic spline.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
require(gps)

## a cubic spline with uniform knots
DemoSpl(uniform = TRUE)

## a cubic spline with non-uniform knots
DemoSpl(uniform = FALSE)
}
