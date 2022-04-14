\name{DemoKnots}

\alias{DemoKnots}

\title{
Demonstrate ordinary cubic B-splines on three types of knots
}

\description{
Demonstrate ordinary cubic B-splines on three types of knots: (a) uniform knots; (b) non-uniform knots; (c) non-uniform knots with clamped boundary knots.
}

\usage{
DemoKnots(aligned = TRUE)
}

\arguments{
\item{aligned}{if TRUE, case (b) and (c) are aligned for nice display, because they have identical interior knots.}
}

\value{
This function has no returned values.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
require(gps)

DemoKnots(aligned = TRUE)
}
