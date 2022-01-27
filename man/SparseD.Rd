\name{SparseD}

\alias{SparseD}
\alias{SparsePD}

\title{
Build general difference matrices
}

\description{
Compute general difference matrices of order 1 to \eqn{d - 1}, for ordinary/periodic B-splines of order \eqn{d}.
}

\usage{
SparseD(xt, d)

SparsePD(xd, d)
}

\arguments{
\item{xt}{full knot sequence for ordinary B-splines (\code{length(xt) >= 2 * d} required).}
\item{xd}{domain knot sequence for periodic B-splines (\code{length(xd) >= d + 1} required).}
\item{d}{B-spline order (\eqn{d \ge 2} required).}
}

\value{
A list of sparse matrices (of "dgCMatrix" class), giving general difference matrices of order 1 to \eqn{d - 1}.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
require(gps)

## example with ordinary B-splines
xt <- c(0, 0, 0, 0, 1, 3, 4, 4, 4, 4)
SparseD(xt, d = 4)

## example with periodic B-splines
xd <- c(0, 0.5, 1.5, 1.8, 3)
SparsePD(xd, d = 4)
}
