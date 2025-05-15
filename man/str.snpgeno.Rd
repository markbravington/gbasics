\name{str.snpgeno}
\alias{str.snpgeno}
\alias{str.loc.ar}
\alias{str.NGS_count_ar}
\alias{str}
\title{Summaries for various genotype classes}
\description{Default \code{str} causes horrible crashes on these objects: too long, or something. So, use these instead.
}
\usage{
\method{str}{snpgeno}( object, loci = TRUE, keys = TRUE, ...) # S3 method for snpgeno
\method{str}{loc.ar}( object, loci = TRUE, keys = TRUE, ...) # S3 method for loc.ar
\method{str}{NGS_count_ar}( object, loci = TRUE, keys = TRUE, ...) # S3 method for NGS_count_ar
}
\arguments{
\item{  object}{of whatever class}
\item{  loci}{TRUE or FALSE to show (some) locus names}
\item{  keys}{TRUE or FALSE to show (most or all) sample-specific fields}
\item{  ...}{ignored AFAIK}
}
\keyword{misc}
