\name{DemoNull}

\alias{DemoNull}

\title{
Demonstrate null space of a P-spline
}

\description{
Demonstrate the limiting behaviour of a P-spline at \eqn{\lambda = +\infty}. A cubic P-spline with a 2nd order difference penalty is fitted to observations simulated from \eqn{g(x) = x}. If it has the correct null space, it will be able to produce a linear line fit with \eqn{\lambda = +\infty}.
}

\usage{
DemoNull(n, k, type = 1)
}

\arguments{
\item{n}{number of simulated observations from \eqn{g(x) = x}.}
\item{k}{number of interior knots to place.}
\item{type}{type of illustrations (1 or 2); see Details.}
}

\details{
Three P-spline variants are experimented:
\itemize{
\item standard P-spline (uniform B-splines with standard difference penalty);
\item naive P-spline (non-uniform B-splines with standard difference penalty);
\item general P-spline (non-uniform B-splines with general difference penalty).
}
The orders of B-splines and difference penalty are chosen to be 4 and 2, respectively.

The function offers two types of illustrations.

In type-1 illustration, \eqn{(x_i)_1^n} is unevenly spaced on \eqn{[0, 1]} (following Beta(3, 3) distribution), and we use automatic knot placement strategy to position \eqn{k} interior knots at equal distances or quantiles. This figure shows that a naive P-spline, with a wrong null space, is unable to produce a linear line fit when \eqn{\lambda = +\infty}.

In type-2 illustration, \eqn{(x_i)_1^n} is equidistant on \eqn{[0, 1]}, and we experiment non-uniform knot sequences that follow different distributions (primarily Beta distributions with varying shape parameters). This figure conveys that a naive P-spline has an undesired, arbitrary limiting behavior because knots can be arbitrary in general.
}

\value{
This function has no returned values.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
require(gps)

## type-1 illustration
DemoNull(n = 400, k = 10, type = 1)

## type-2 illustration
DemoNull(n = 400, k = 50, type = 2)
}
