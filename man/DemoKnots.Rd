\name{DemoKnots}

\alias{DemoKnots}

\title{
Demonstrate ordinary cubic B-splines on three types of knots
}

\description{
Demonstrate ordinary cubic B-splines on three types of knots: (a) uniform knots; (b) non-uniform knots; (c) non-uniform knots with clamped boundary knots. The same interior knots are positioned in cases (b) and (c).
}

\usage{
DemoKnots(aligned = TRUE)
}

\arguments{
\item{aligned}{if TRUE, interior knots in cases (b) and (c) are aligned for a better display.}
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
