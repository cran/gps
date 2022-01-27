\name{pbsDesign}

\alias{pbsDesign}

\title{
Build design matrix for periodic B-splines
}

\description{
Evaluate periodic B-splines or their derivatives at given x-values.
}

\usage{
pbsDesign(x, xd, d, nDeriv = 0, sparse = FALSE)
}

\arguments{
\item{x}{x-values where periodic B-splines are to be evaluated.}
\item{xd}{domain knot sequence for periodic B-splines (\code{length(xd) >= d + 1} required).}
\item{d}{B-spline order (\eqn{d \ge 2} required).}
\item{nDeriv}{number of derivatives to take.}
\item{sparse}{if TRUE, create a sparse design matrix of "dgCMatrix" class.}
}

\details{
Unlike ordinary B-splines, period B-splines do not require explicit auxiliary boundary knots for their construction. The magic is that auxiliary boundary knots will be automatically positioned by periodic extension of interior knots.

Denote the domain knot sequence by \eqn{s_0, s_1, s_2, \ldots, s_k, s_{k + 1}}, where \eqn{(s_j)_1^k} are interior knots and \eqn{s_0 = a}, \eqn{s_{k + 1} = b} are domain boundaries. For order-\eqn{d} B-splines, we replicate the first \eqn{d - 1} interior knots (after adding \eqn{b - a}) to the right of \eqn{[a, b]} for an augmented set of \eqn{K = k + d + 1} knots, which spawns \eqn{p = K - d = k + 1} ordinary B-splines. It turns out that periodic B-splines can be obtained by wrapping segments of those ordinary B-splines that stretch beyond \eqn{[a, b]} to the start of the domain (a demo is offered by \code{\link{DemoPBS}}).

Note that we must have at least \eqn{d - 1} interior knots to do such periodic extension. This means that \eqn{d + 1} domain knots are required at a minimum for construction of periodic B-splines.
}

\value{
A \eqn{n \times (k + 1)} design matrix, where \eqn{n} is the number of x-values (i.e., \code{length(x)}) and \eqn{k} is the number of interior knots (i.e., \code{length(xd) - 2}).
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
require(gps)

## 5 domain knots: three interior knots 0.5, 1.5 and 1.8 in domain [0, 3]
xd <- c(0, 0.5, 1.5, 1.8, 3)

## make a grid
x <- MakeGrid(xd, n = 10)

## construct periodic cubic B-splines
PB <- pbsDesign(x, xd, d = 4)
}
