\name{diploido}
\alias{diploido}
\title{Class for one diploid locus}
\description{The ancient and probably obsolete class \code{diploido} lets a diploid locus be treated as a 1D object (a vector), even though it's really a 2D matrix; no, I can't remember exactly why. The usual generics are available. Contrast with \code{\link{loc.ar}}, which works with entire data.frames most of which consists of a multilocus genotype.
}
\usage{diploido(cop1, cop2)
}
\arguments{
\item{  cop1}{integer vector}
\item{  cop2}{integer vector}
}
\seealso{loc.ar
}
\keyword{misc}
