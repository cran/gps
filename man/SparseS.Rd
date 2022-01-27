\name{SparseS}

\alias{SparseS}

\title{
Build derivative penalty matrix or its sparse "root"
}

\description{
Compute order-\eqn{m} derivative penalty matrix \eqn{\bm{S}} or its sparse "root" \eqn{\bm{K}} such that \eqn{\bm{K'K} = \bm{S}}.
}

\usage{
SparseS(xt, d, m, root = FALSE)
}

\arguments{
\item{xt}{full knot sequence for ordinary B-splines (\code{length(xt) >= 2 * d} required).}
\item{d}{B-spline order (\eqn{d \ge 2} required).}
\item{m}{order of the derivative penalty (\eqn{0 \le m \le d - 1} required).}
\item{root}{if TRUE, return the sparse "root" of the derivative penalty matrix.}
}

\value{
A sparse matrice giving the order-\eqn{m} derivative penalty matrix (of "dsCMatrix" class) or its sparse "root" (of "dtCMatrix" or "dgCMatrix" class).
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
require(gps)

## example with ordinary B-splines
xt <- c(0, 0, 0, 0, 1, 3, 4, 4, 4, 4)
SparseS(xt, d = 4, m = 2)
}
