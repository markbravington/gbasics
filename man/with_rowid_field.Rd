\name{with_rowid_field}
\alias{with_rowid_field}
\alias{rowid_field}
\title{Facilitate sample-based subscripting of sngpenos}
\description{
\code{with_rowid_field} can be applied to a \code{\link{snpgeno}} to specify which "sample ID field" to use when subsetting a \code{\link{snpgeno}} by row (ie sample), using character "sample ID" of your choice instead of numeric or logical index.

You can already subscript a \code{\link{snpgeno}} by locus with a character index, which is looked up in \code{<snpgeno>@locinfo$Locus}. The choice of "Locus" for field name is pretty obvious and uncontroversial! But for samples, there's no corresponding uniquely obvious field name so you can use \code{with_rowid_field} to specify the one you want.

Specifically, that field in \code{<snpgeno>$info} will be used to look up character-mode \code{i} in \code{<snpgeno>[i,...]} subsetting. Then you can subset to samples using their "identifiers", rather than having to use logical or numeric lookup. Also, if \code{with_rowid_field} has been called, then subsequent calls to \code{find_HSPs} (qv) and friends will return the rowid_fields in \code{i} and \code{j} rather than row-numbers, which is much less error-prone for downstream use.

You certainly don't have to call \code{with_rowid_field} on your \code{\link{snpgeno}} (it didn't exist until 2023...), but unless you do, you'll be stuck with logical/character subsetting by sample. Of course, you can still use those afterwards, as well.

\code{rowid_field} is a convenience lookup function to remind yourself which field you specified...
\subsection{Details}{
\code{with_rowid_field} is actually a generic which can mark two (currently) classes of object: \code{\link{snpgeno}} and \code{data.frame} (or descendents of the latter). If you call it on a \code{\link{snpgeno}}, it actually gets applied to \code{<snpgeno>$info} which is a \code{data.frame}; this means that you can subset not just \code{<snpgeno>} itself, but also \code{<snpgeno>$info}, with character sample info. Note that this might theoretically \emph{contradict} the default behaviour of a \code{data.frame}, which is to look up character row indices in the highly-unreliable \code{rownames(<data.frame>)}; I force rownames \emph{off} (ie NULL) in a \code{<snpgeno>$info} so this shouldn't be a problem.

\code{with_rowid_field} augments the S3 class of a \code{data.frame}, but not of a \code{\link{snpgeno}} (the functionality of using rowid is built into the latter, and the rowid label is attached to \code{<snpgeno>$info} not to \code{<snpgeno>} itself). I know it does say "DETAILS" here but there are limits, so I'm not going to explain why.

There are methods for \code{[}, \code{[<-}, \code{rbind}, and \code{cbind}, which attempt to do the right thing...
}
}
\usage{
with_rowid_field( x, rowid_field) # generic
rowid_field( x) # generic
}
\arguments{
\item{  x}{thing that you'll want to do character-based row subsets on}
\item{  rowid_field}{which column name to look up.}
}
\value{\code{with_rowid_field} returns the original object marked (somewhere in its internals) with a \code{rowid_field} attribute. \code{rowid_field} returns a string, or NULL.
}
\examples{
\dontrun{
sn <- snpgeno(...) # with "sampID" as a field in 'info'
try( sn[ 'samp_X12',]) # crash
sn <- with_rowid_field( sn, 'sampID')
sn[ 'samp_X12',] # goodo
rowid_field( sn) # [1] "sampID"
rowid_field( sn@info) # ditto
sn@info[ 'samp_X12', ] # whatever
}
}
\keyword{misc}
