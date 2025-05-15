\name{get_genotype_encoding}
\alias{get_genotype_encoding}
\alias{define_genotypes}
\title{Diploid genotype encodings}
\description{
Many of the functions in \code{kinference} etc can accept several different encodings for diploid genotypes, depending on how the genotyeps are distinguished. For example, one encoding specifies that single-nulls are called separately from homozygotes (perhaps with error), whereas another says that single-nulls and homozygotes are not distinguished; some encodings allow a 3rd or 4th etc option for the allele; etc.

Encoding is specified by the \code{\link{diplos}} attribute in the \code{\link{snpgeno}} object. The hope is that it should normally be handled for you automatically, but you can set the encoding manually to one of a few pre-specified options that \code{kinference} will understand; the choices are shown by \code{get_genotype_encoding()}, though only a few will actually work with the \pkg{kinference} package. An encoding is stored as a character vector of recognized genotype-categories; some are described next.
\subsection{Some popular encodings}{
This is probably the \emph{wrong} place to describe the various encodings currently used by \code{kinference} (a few others are or were used by CSIRO's own genotyping pipeline before the results are \code{kinference}-ready). Anyway, the commonest ones are:

\itemize{
\item 4-way (\code{genotypes4_ambig}) allows the values "OO", "AB", "AAO", "BBO". This means nulls are allowed, double-nulls ("OO") are taken notice of, but single-nulls and homozygotes are not differentiated, so that "AAO" means "A" was present but "B" was not, and conversely for "BBO". You can use this encoding even if your data is guaranteed null-free; there just won't be any "OO", and the null-allele frequency can be set to 0, which means \code{kinference} will always interpret "AAO" as a homozygote.
\item 3-way (\code{genotypes3}) comprises "AB", "AAO", and "BBOO" where the latter covers "BB" homozygotes and "BO" single-nulls \emph{and} "OO" double-nulls. The motivation is for loci where nulls are rare but not non-existent, and where double-nulls (which \emph{might} in practice be artefacts, although we assume they aren't) get wildly over-interpreted in a kinship setting. Merging genotype-categories never causes bias, but it does sacrifice a small amount of statistical information to gain (much more important) robustness. Nulls can actually be informative, so "blanket 3-way" is mostly not the default for CSIRO's datasets; problems only come when nulls are rare. So CSIRO has sometimes used \code{genotypes4_ambig} for all loci, but setting the \code{locinfo$useN} field to 3 for loci where double-nulls are too risky. Various functions in \code{kinference} look at \code{locinfo$useN} to decide whether an individual locus is treated \emph{during analysis} as 3-way, 4-way, or perhaps 6-way (next) rather than following the general encoding of the whole dataset.
\item 6-way (\code{genotypes6}) with possibilities "AA", "AB", "AO", "BB", "BO", and "OO". Here, the genotyping process has tried to distinguish between single-nulls and true homozygotes, based on read-depth. The discrimination will never be 100\% accurate, even for "good" loci, so explicit allowance must be made for getting it wrong sometimes (if you see references to \code{snerr}, that's why). For some loci it's not even worth trying, so for them \code{locinfo$useN} is set to 4 or 3 (and it is set to 6 for loci where 6-way is really being used). The pipeline for 6-way genotyping uses in-house CSIRO code and is quite fiddly--- so it's likely to \emph{stay} in-house! CSIRO has used 6-way genotyping for \emph{large} ongoing studies of SBTuna and school sharks, so this encoding will continue to be supported, but for new projects we would recommend 4-way instead, using more loci to compensate.
}

You will see that 4-way is a coarse-grained version of 6-way, and 3-way is a coarsened version of 4-way. The \code{locinfo$useN} field (a column of numbers) can be used to selectively coarsen certain loci, but of course it can't "uncoarsen" any part of a 3-way-encoded dataset if the overall encoding was too coarse.
}
\subsection{Internal representation}{
Internally (though you \emph{shouldn't} need to know this), the genotypes in a \code{\link{snpgeno}} are stored in a \code{raw}-mode matrix, where each element takes 1 byte. The elements of the encoding correspond in order to the numbers 1,2,3,.... With \code{genotypes4_ambig}, for example, an element with raw value 1 means "OO", with value 2 means "AB", etc. Value 0 should never be seen in real datasets, but is sometimes useful as a temporary placeholder to show that "something needs fixing"! It will display as "?".

Generally speaking, \code{sngpeno} objects use raw value \code{ff} (255) to indicate "missing" (as distinct from double-null). However, the \pkg{kinference} package deliberately does not handle missing data; \emph{all} genotypes have to get called, even though some of those calls might be in error, and genotyping errors are allowed for (or ignored) in the subsequent statistical analyses.

Special S3 methods exist for printing, comparing, and assigning genotype encodings, so it should at least \emph{look} clear to you, even if you have no idea what all this is about. By design, you can't directly assign integer values into genotypes; see \bold{Examples}.
}
\subsection{Changing encodings}{Mostly, encoding should be handled automatically during the creation of your \code{\link{snpgeno}}. However, you do sometimes need do access to it, and \emph{occasionally} you might want to manually recode your data. That takes some care because the internal representation of the encodings won't generally match up; for example, "AB" is encoded as 1 for \code{genotypes3} but as 2 in \code{genotypes4_ambig}. You can use \code{get_genotype_encoding} to change the encoding, then manually change old genotypes to new ones. Here's an example that changes all genotypes from 3-way to 4-way (which you also could do less painfully via \code{useN}--- hopefully that is documented somewhere).
}
\subsection{Expandability}{
You can in fact use any character vector you like; there is no need to stick to those provided by \code{define_genotypes}, and the interpretation of the character-strings is entirely at the discretion of downstream software. For example, you could specify your own encoding scheme for microhaplotypes; from memory, I think up to 7 variants could be accommodated within the 8-bit \code{raw} storage. However, most of the \code{kinference} functions can only handle one or two of the pre-specified encodings, so you'd need new software.

.

\code{define_genotypes} is meant only for use inside other functions, and to help explain the flexible encoding system for genotype storage used by \code{\link{snpgeno}}. Once you've called it inside a function, you can refer directly to a genotype \code{"AB"} as just \code{AB} etc (ie no quotes), and you can also refer to sets of genotype-encoding that are commonly used in \code{\link{snpgeno}} objects; those sets are just character vectors. It creates objects such as \code{AB} \bold{in its caller}, so a lot of kinference functions run \code{define_genotypes} internally as a first step to simplify their subsequent code- and that is the intended use.

If \emph{you} run \code{define_genotypes()} from the R{} prompt, those things will be created in \code{.GlobalEnv}. So, most people shouldn't. There is usually no specific need to run \code{define_genotypes} in your own code either, unless you're aiming to extend \code{kinference}; the point of documenting it here, is to explain and illustrate how the genotype encodings work. For All Practical Purposes, \code{get_genotype_encoding} should be enough.
}
}
\usage{
# Really meant to be used early inside another function
define_genotypes(nlocal = sys.parent())
get_genotype_encoding()
}
\arguments{\item{  nlocal}{the frame number, or environment, to create things in- see \code{\link{mlocal}}. Leave this alone unless you \emph{really} know what you're doing...}
}
\value{
\code{get_genotype_encoding} returns a list of character vectors. Their names aren't "sacred"; the only thing that matters to \code{kinference} is the character-vector of strings in the encoding.
\code{define_genotypes} returns nothing, but various objects are created; see the code.
}
\seealso{\code{\link{snpgeno}}; \code{cq} and \code{mlocal} in package \pkg{mvbutils}
}
\examples{
# Use of define_genotypes() inside a function:
library( mvbutils) # just for '\%is.a\%'
count_hetz <- function( x)\{
    stopifnot( x \%is.a\% 'snpgeno')
    define_genotypes()
    # Now eg genotypes4_ambig is directly available...
    # and so is AB instead of "AB", etc. Cor blimey.
  return( sum( x==AB)) # thus saving a pair of quotes...
  \}
data( snpgarbage)
count_hetz( snpgarbage)
sum( snpgarbage=='AB') # had to type two quotes...
# Recode snpgarbage as 3-way (no good reason)
get_genotype_encoding() # for inspection; returns them all
diplos( snpgarbage) # aha! Looks like 'genotypes4_ambig'
# Make a copy; *don't* try to modify-in-place!
sg3 <- snpgarbage
sg3[] <- as.raw( 0) # will flag an error later if we forget something
# sg3[] <- 0 won't work; guard against user boo-boo
# Set the desired encoding...
diplos( sg3) <- get_genotype_encoding()$genotypes3
# ... and re-map _all_ values from the old encoding
# Note that AB and AAO may not use the same raw code in both
# encodings, so you gotta explicitly set those ones too
sg3[ snpgarbage=='OO'] <- 'BBOO'
sg3[ snpgarbage=='AAO'] <- 'AAO'
sg3[ snpgarbage=='AB'] <- 'AB'
sg3[ snpgarbage=='BBO'] <- 'BBOO'
snpgarbage # out with the old...
sg3 # ... in with the new
}
\keyword{misc}
