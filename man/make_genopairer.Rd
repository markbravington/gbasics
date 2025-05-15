\name{make_genopairer}
\alias{make_genopairer}
\title{Mapping between pairs of possible genotypes and compressed form}
\description{
You almost certainly do not need, nor want, to call this. But it needs to be exported so that \code{kinference} and other packages can find it.

With diploid genotypes such as "AA" and "AB" at a locus, the pairwise probs under given co-inheritance (0, 1, or 2 copies) are symmetric. It is daft to store the full symm matrix, and causes some boring technical problems too. So, \code{make_genopairer} creates a mapping between a vector of pair-names such as "AA/AB" and the underlying pairs. It allows eg a 2000x6x6 array to be stored as 2000x15, and makes passing it to Rcpp faaaaar easier.
}
\usage{make_genopairer(genotypes)
}
\arguments{\item{  genotypes}{character vector, produced by e.g. \code{define_genotypes}, which you probably shouldn't be calling either.}
}
\value{Matrix with rows & columns being the individual genotypes, and entries a lookup from that pair into a character vector of possible \emph{ordered} pairs, which is stored as the "what" attribute. See \bold{Examples}.
}
\examples{
make_genopairer( c( 'AA', 'AB', 'BB'))
#   AA AB BB
#AA  1  2  3
#AB  2  4  5
#BB  3  5  6
#attr(,"what")
#[1] "AA/AA" "AB/AA" "BB/AA" "AB/AB" "BB/AB" "BB/BB"
}
\keyword{misc}
