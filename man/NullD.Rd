\name{NullD}

\alias{NullD}

\title{
Find null space of a general difference matrix
}

\description{
Compute a basis matrix for the null space of a general difference matrix.
}

\usage{
NullD(xt, d, m)
}

\arguments{
\item{xt}{full knot sequence for ordinary B-splines (\code{length(xt) >= 2 * d} required).}
\item{d}{B-spline order (\eqn{d \ge 2} required).}
\item{m}{order of the general difference matrix (\eqn{1 \le m \le d - 1} required).}
}

\value{
A matrix of \eqn{m} columns.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
require(gps)

xt <- c(0, 0, 0, 0, 1, 3, 4, 4, 4, 4)
NullD(xt, d = 4, m = 2)
}
