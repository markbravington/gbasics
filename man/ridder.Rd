\name{ridder}
\alias{ridder}
\title{Parallel root-finder}
\description{
Ridder's root-finder for many univariate functions in parallel. Faster than looping over each function in turn, IME.

Ridder's method is a very good 1D root-finding algorithm; see \code{https://doi.org/10.1109\%2FTCS.1979.1084580} for original, or the section in "Numerical Recipes" (probably chapter 9, for the 2007 edition) by Press, Teukolsky, Vetterling, Flannery.

One weakness of this implementation is that all components of \code{FUN} get evaluated in each iteration, until the very last component has converged. It would be faster if I allowed \code{FUN} to have a logical-index argument, saying which components to actually evaluate each time.
}
\usage{ridder(FUN, lo = NULL, hi = NULL, tol = 0.000001, skip_bounds = FALSE)
}
\arguments{
\item{  FUN}{function with one \emph{vector} argument. If your underlying function has other args, you need to wrap it first so that you can pass a one-argument function to \code{ridder}.}
\item{  lo, hi}{bounds for search that should bracket the roots (can be vectors)}
\item{  tol}{for roots (absolute tolerance- a scalar, duhhh)}
\item{  skip_bounds}{if TRUE, start the searches exactly at the bounds; quicker, but a bad idea if infinite values result. If FALSE, the search will start inside the bounds and automatically do some bracketing, being sure to avoid the bounds.}
}
\value{The roots (a vector).
}
\keyword{misc}
