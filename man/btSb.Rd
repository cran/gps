\name{btSb}

\alias{btSb}

\title{
Evaluate derivative penalty 
}

\description{
Evaluate derivative penalty \eqn{\int_a^b f^{(m)}(x)^2\mathrm{d}x = \bm{\beta'S\beta}}, where \eqn{\bm{S}} is the order-\eqn{m} derivative penalty matrix and \eqn{\bm{\beta}} is the vector of B-spline coefficients in \eqn{f(x)}.
}

\usage{
btSb(b, xt, d, m)
}

\arguments{
\item{b}{B-spline coefficients (\code{length(b) == length(xt) - d} required).}
\item{xt}{full knot sequence for ordinary B-splines (\code{length(xt) >= 2 * d} required).}
\item{d}{B-spline order (\eqn{d \ge 2} required).}
\item{m}{order of the derivative penalty (\eqn{0 \le m \le d - 1} required).}
}

\details{
Given that \eqn{\bm{S}} can be obtained from \code{\link{SparseS}}, it is straightforward to compute the quadratic form \eqn{\bm{\beta'S\beta}}. However, the function takes an alternative approach.

Denote the domain knot sequence by \eqn{s_0, s_1, s_2, \ldots, s_k, s_{k + 1}}, where \eqn{(s_j)_1^k} are interior knots and \eqn{s_0 = a}, \eqn{s_{k + 1} = b} are domain boundaries. The derivative penalty aggregates roughness penalties on all knot spans: \deqn{\int_a^b f^{(m)}(x)^2\mathrm{d}x = \sum_{j = 0}^k\int_{s_j}^{s_{j + 1}} f^{(m)}(x)^2\mathrm{d}x.} The function calculates and returns those additive components in a vector.
}

\value{
A vector of evaluated derivative penalties at all knot spans.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
require(gps)

## example with ordinary B-splines
xt <- c(0, 0, 0, 0, 1, 3, 4, 4, 4, 4)
b <- rnorm(6)
btSb(b, xt, d = 4, m = 2)
}
