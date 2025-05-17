\name{read_snpgds2snpgeno}
\alias{read_snpgds2snpgeno}
\title{Create a 'snpgeno' from 'snpgds' file or object}
\description{\code{read_snpgds2snpgeno} is a wrapper allowing one-line conversion from a \code{snpgds}-format \emph{file} to a \code{\link{snpgeno}}. Needs the \pkg{SNPRelate} package.
}
\usage{
read_snpgds2snpgeno(
    filename,
    locusID,
    sampleID,
    infoFrame = NULL,
    infoFields = NULL,
    locinfoFrame = NULL,
    locinfoFields = NULL,
    plateField = NULL)
}
\arguments{
\item{  filename}{string giving the snpgds file name, contents formatted as per \code{snpgdsExampleFileName} in package \pkg{SNPRelate}.}
\item{  locusID}{string namgiving the name for locus IDs. Cannot be blank}
\item{  sampleID}{string naming the column for sample IDs. Cannot be blank}
\item{  infoFrame}{optional string naming the sample metadata dataframe}
\item{  infoFields}{optional character vector naming the sample metadata variables}
\item{  locinfoFrame}{optional string naming the locus metadata dataframe}
\item{  locinfoFields}{optional character vector naming the locus metadata variables}
\item{  plateField}{optional string naming the column for sample-specific plate ID}
}
\details{Locus-specific metadata and sample-specific metadata may each be supplied as either a single dataframe (using a non-NULL \code{infoFrame} and/or \code{locinfoFrame}), or as a vector of field names (using \code{infoFields} and/or \code{locinfoFields}). If all are NULL, no metadata other than the locus ID and sample ID are read in.
}
\examples{
# From SMB:
#if (!requireNamespace("BiocManager", quietly=TRUE))
#    install.packages("BiocManager")
# BiocManager::install("SNPRelate")
if( requireNamespace( 'SNPRelate'))\{
  # Simplest possible case (no locus or sample metadata other than IDs;
  # will add an all-one plate field):
  sg <- read_snpgds2snpgeno(filename = SNPRelate::snpgdsExampleFileName(),
      locusID = "snp.id", sampleID = "sample.id")
  # More likely: either sample or locus metadata is a single frame,
  # and the other is a bunch of separate fields. In this case, sample
  # metadata is a single frame and locus metadata is separate fields.
  sg <- read_snpgds2snpgeno(filename = SNPRelate::snpgdsExampleFileName(),
     locusID = "snp.id", sampleID = "sample.id",
     infoFrame = "sample.annot",
     locinfoFields = c("snp.rs.id", "snp.position", "snp.chromosome", "snp.allele"))
\}
}
\keyword{misc}
