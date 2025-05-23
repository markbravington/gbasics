\name{NGS_count_ar}
\alias{NGS_count_ar}
\alias{NGS_count_ar.matrix}
\alias{print.NGS_count_ar}
\alias{print}
\title{Next-Gen sequencer count data}
\description{
\code{NGS_count_ar} is a class for raw counts by sequence-variant within locus, per sample. The number of alleles (variants) can differ between loci. It's not used by any functions in the \pkg{kinference} package, which expects \code{\link{snpgeno}} objects throughout. However, it's used extensively in CSIRO's private \code{genocalldart} package. It's an early stage in processing DartSeq/DartCap/DartTag data--- long before genotypes are called--- and can also be used to hold similar data from other formats such as VCF, via \code{read_vcf2NGS_count_ar()}. It's also a generic S3 method for creating such an object, with one low-level method for matrices (which you probably shouldn't call yourself) and another method for converting from \code{\link{loc.ar}} objects, which will presumably have come from \code{genocalldart::read_count_dart} or such. The default method (ie if \code{NGS_count_ar} is called on any old rubbish) will \code{stop()}, by design. See \bold{Details} for internal structure, and what you can add/expect.

There are methods for basic stuff: try \code{methods(class="NGS_count_ar")}. \code{dim(x)} returns \code{#samples*#loci}; for the number of alleles, use \code{nrow( x$seqinfo)}. The \code{print} method tries to make things clearly readable, and has a few extra arguments to help control that- see \bold{Arguments}. Subsetting can have up to 3 indices \code{i}, \code{j}, and \code{k}; \code{j} can be numeric, logical, or character (i.e. locus names, which are matched against \code{x$locinfo$Locus}). If \code{k} is omitted, another \code{NGS_count_ar} object is returned, as controlled by \code{i} and/or \code{j}. If \code{length(k)==1}, then the result is a 2D array of sample-by-locus. If \code{length(k)>1}, the result is a 3D array of sample-by-locus-by-allele, with NA counts added when \code{k} refers to an allele "number" that doesn't exist for that locus; obviously, all loci then have the same number of "alleles" in the new structure. \code{dimnames} are set for the 2nd (locus) dimension only. These "pure array" results don't preserve the detailed sample or locus information (except via the \code{dimnames}) but can be useful for subsequent manipulation.
\subsection{Obscure note}{A similar class is \code{\link{loc.ar}}, which preserves counts but allows max 3 non-null alleles per locus. It's not clear which of \code{\link{loc.ar}} or \code{NGS_count_ar} is "more advanced"- depends on the context. However, most later stages in \code{genocalldart} pipeline (eg for bait selection) require \code{\link{loc.ar}}. To go back the other way, you probably need to call \code{pick_ref_alt} (qv); see the pipeline examples.
}
}
\usage{
NGS_count_ar( x, ...) # generic
\method{NGS_count_ar}{matrix}( x,  # S3 method for "matrix'
    sampinfo, locinfo, seqinfo,
    strip_numerics = TRUE,  rename = TRUE, ...) # S3 method for matrix
\method{print}{NGS_count_ar}( x,  # S3 method for "NGS_count_ar'
    trailing_dot = getOption("trailing_dot_NGS_count_ar", FALSE),
    dot_for_0 = getOption("dot_for_zero_NGS_count_ar", FALSE),
    center_dot=centre_dot,
    centre_dot= getOption( 'centre_dot', getOption( 'center_dot', FALSE)),
    ...) # S3 method for NGS_count_ar
}
\arguments{
\item{  x}{for the constructor, an integer matrix with dimensions (allele,sample), i.e. the \emph{transpose} of final storage. Or a \code{loc.ar}. For \code{print}, an \code{NGS_count_ar}}
\item{  sampinfo}{dataframe of sample info}
\item{  locinfo}{dataframe of locus info}
\item{  seqinfo}{dataframe of sequence info (including which locus the sequence belongs to)}
\item{  strip_numerics}{if TRUE, remove any non-integer columns from \code{locinfo}- typically, these are summary data in a CSV that you can live without}
\item{  rename}{certain names are expected by fun/xctions that handle \code{NGS_count_ar} objects; this patches up "alternative" names that are generated by \code{read_cluster_dart3}}
\item{  trailing_dot}{set TRUE if you want to show that count data is non-integer (eg after norming by sample-total-reads) so that all counts end with a period. The post-decimal-place digits are not normally important, but if you really need to see them, you can do so by setting a \code{k} subscript. Or set the global option, as per \bold{Usage}.}
\item{  dot_for_0}{set TRUE if you want zero-counts replaced by a dot. Using the global option is often helpful. Or set the global option, as per \bold{Usage}.}
\item{  center_dot, centre_dot}{iff \code{dot_for_0} is TRUE, whether to use a central dot (Latin-1 and Unicode 0xb7) or a period (if FALSE) to replace leading zeros.}
\item{  ...}{passed to other methods}
}
\details{
The allele counts are stored as a (sample*sequence) matrix. There are also three key attributes, all accessible using the \code{$} operator, and each a dataframe:

\describe{
\item{info}{Metadata on each sample, including "Our_sample" and, usually for Dart data at least, "TargetID"}
\item{locinfo}{Metadata on each locus; see below}
\item{seqinfo}{Metadata on each allele, including "Locus", "FullAlleleSeq", and probably others}
}

The key fields in \code{x$locinfo} (there might be others) are

\describe{
\item{Locus}{Unique strings}
\item{n_alleles}{Number of sequence variants (presumably 2 or more)}
\item{consensus}{String of the overall sequence, with dots at SNP sites}
\item{var_pos}{String showing the positions of SNP sites, eg "17,41"}
\item{end_col}{of the sequences found at that locus (needed for lookup into the main matrix)}
}

You can also add other attributes that apply to the whole object. One that's used in \code{genocalldart} is \code{mean_fish_tot} (used in norming count data from new samples), but it's rather dicey for general use because it does depend on the set of loci; so if you subset by loci, things (should) change...
}
\seealso{\code{\link{loc.ar}}, \code{\link{snpgeno}}
}
\keyword{misc}
