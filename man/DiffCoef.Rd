\name{DiffCoef}

\alias{DiffCoef}

\title{
Compute general differences of B-spline coefficients
}

\description{
Compute general differences of B-spline coefficients in \eqn{f(x)}, which are coefficients of order-\eqn{(d - m)} B-splines in \eqn{f^{(m)}(x)}.
}

\usage{
DiffCoef(b, xt, d, m)
}

\arguments{
\item{b}{B-spline coefficients (\code{length(b) == length(xt) - d} required).}
\item{xt}{full knot sequence for ordinary B-splines (\code{length(xt) >= 2 * d} required).}
\item{d}{B-spline order (\eqn{d \ge 2} required).}
\item{m}{order of the general difference matrix (\eqn{1 \le m \le d - 1} required).}
}

\value{
A vector of coefficients.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
require(gps)

## example with ordinary B-splines
xt <- c(0, 0, 0, 0, 1, 3, 4, 4, 4, 4)
b <- rnorm(6)
DiffCoef(b, xt, d = 4, m = 2)
}
