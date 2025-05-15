\name{loc.ar}
\alias{loc.ar}
\alias{unloc.ar}
\title{Create & manipulate locus-arrays}
\description{
A \code{loc.ar} object is meant for "numerical diploid genotype data". The class is \emph{not} used by the \pkg{kinference} package, is not particularly well thought out, and may not be properly documented here. We no longer use \code{loc.ar}s much at CSIRO except as intermediate steps in constructing a \code{\link{snpgeno}}, but in the past we have used them for: (i) microsat genotypes, with each sample-locus being a pair of numbers for allele-length; and (ii) biallelic SNP read-depth pairs, with the genotypes (if decided) stored as a hidden extra. The latter purpose is generalized to "microhaplotypes" (multiple alleles per locus) by \code{\link{NGS_count_ar}} (qv), which is probably more "future-proofed". \code{loc.ar} probably won't go away from \code{\link{gbasics}}, but it is not guaranteed to be actively maintained (whereas \code{\link{snpgeno}} and \code{\link{NGS_count_ar}} probably will be); then again, perhaps it will be improved in future. Who knows?

There are methods for the usual things like subsetting (see \bold{Details}), \code{rbind}, \code{cbind}, \code{\link{print}}, \code{head}, \code{tail}, \code{\link{str}}. There is also a constructor \code{loc.ar}- an S3 generic- but currently there are \emph{no} methods for it (at least not in package \pkg{gbasics}; a dodgy one has been moved to \code{genocalldart} for now). Any \code{loc.ar} object has to be built manually, as done by a few functions in package \pkg{genocalldart}, such as \code{read_count_dart}.

At heart, a \code{loc.ar} is a 3D numeric (or integer) array with dimension \code{(n_samples,n_loci,X)} where X is a fixed maximum number of non-null alleles for each locus (always either 2 or 3 in our applications to date), but with attributes that contains auxiliary data: sample-specific, locus-specific, optionally sample-by-locus, plus whatever other fixed attributes you want to give it. It is broadly similar to \code{\link{snpgeno}}, which is better documented, except for the details of its contents (3D numeric, vs 2D \code{raw}-mode interpreted as actual called genotypes). Currently, there are methods for printing, subsetting, rbind/cbind, and for transforming to/from dataframes with two cols per locus. Subsetting can handle the auxiliary data automatically, as described below.

Auxiliary data is stored in attributes which can be accessed via the \code{$} operator. Sample-specific auxiliary info (which must exist) can be accessed via \code{x$info}; it should be a dataframe containing a column "Our_sample". Locus-specific auxiliary info is optional (but you'd be brave not to have it...); it should live in a dataframe containing one column "Locus" (the name), and is stored in \code{x$locinfo}. Loci are \emph{rows} here, whereas loci are \emph{columns} in \code{x} itself. Other than subsetting, there are currently no special operations on \code{locinfo}; for example, \code{print( x)} does not show it (whereas, for \code{\link{snpgeno}} objects, the field "Locus" is shown vertically as the column name).

Further auxiliaries, such as SNP genotypes (sample-by-locus), can also be tacked onto the attributes as you please. All attributes are copied by subsetting. By default, the whole attribute is copied, except for (i) \code{info} and \code{locinfo}, which subset automatically in the appropriate way, and (ii) any auxiliaries named in the attribute \code{x$subset_like_both} (a character vector that you can set manually), which should be arrays/matrices whose first two dimensions pertain to sample and to locus respectively (just like the main data in a \code{loc.ar}); they can have more dimensions, too.

I sometimes add class \code{specialprint} from package \pkg{mvbutils} to \code{loc.ar} objects, to \emph{suppress} printing of some auxiliaries. See EXAMPLES.
\subsection{Internal note}{I used to use class \code{\link{diploido}} to hold microsat genotypes inside dataframes, but that was less flexible. The \code{genocalldart::SBTlike_loc.ar} creator function can take input from "diploido"-class objects in dataframes, rather than just dataframes with e.g. "D225_A.1" and "D225_A.2" columns. If any column is of class "diploido", all are assumed to be, and "normal" loci are ignored.
}
}
\usage{
loc.ar( x, ...) # generic
unloc.ar( la)
}
\arguments{
\item{  x}{in theory, something to be converted into a \code{loc.ar} (but there are no methods yet!)}
\item{  la}{a \code{loc.ar} with exactly 2 alleles per locus (thus an \code{S*L*2} array) to be turned into a dataframe with the \code{info} columns, plus two columns for each locus. Why, you may ask?}
\item{  ...}{just in case}
}
\details{
Sample info and locus info are stored as dataframes in \code{attr(<x>,"info")} and \code{attr(<x>,"locinfo")} respectively. Internally, the code uses the \pkg{atease} package so it can just write e.g. \code{x@info} instead, but you should not need to load \code{atease} because you can always use the \code{$} and \code{$<-} accessors instead.

For subsetting, note that \code{k}, or \code{j} \emph{and} \code{k}, can be omitted, in which cases \code{i} and \code{j}, or \code{i} alone, will be used as the \emph{first} subscript; this is \emph{unlike} dataframes, where a single subscript gets used as the \emph{second} subscript. If the resulting array contains only one locus, it will by default be collapsed to an n-by-2 matrix, unless \code{drop=FALSE} explicitly. Otherwise, the result will stay as a \code{loc.ar} provided the 3rd dimension stays at length 2. The first dimension (samples) is never dropped.

\code{unloc.ar} should probably acquire a "diploido.use" argument, to optionally return loci as \code{\link{diploido}} objects rather than pairs of columns.
}
\value{
\item{  loc.ar, [, rbind}{a \code{loc.ar}}
\item{  unloc.ar}{a dataframe with two cols per locus}
}
\seealso{NGS_count_ar, snpgeno
}
\keyword{misc}
