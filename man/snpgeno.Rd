\name{snpgeno}
\alias{snpgeno}
\alias{snpgeno.default}
\alias{diplos}
\alias{diplos<-}
\alias{set_rowid_field}
\title{Class for SNP genotypes}
\description{
\code{snpgeno} is an S3 class for storing \emph{already-called} genotypes from multiple samples and loci, as well as associated information about those samples and loci. Pretty much all the functions in the \pkg{kinference} package expect a \code{snpgeno} input, in many cases requiring extra information about the loci, such as allele frequency estimates, which has usually been added by other functions in \code{kinference}. You can extend the per-sample and per-locus information by adding fields in the usual R{} manner, as well as adding extra information that is per-sample-and-locus (such as number-of-reads); there is generally no need to make an "inherited class" for such ad-lib extensions. Printing looks IMO nice, and is succinct. See how elegantly the locus-names are shown (vertically), to save space!

There are a couple of ways to create a \code{snpgeno}. One is to call \code{\link{read_vcf2snpgeno}} (qv) on a VCF file that already contains called genotypes, or \code{\link{read_snpgds2snpgeno}} (qv) for a SNPDGS file. Another is via the constructor \code{snpgeno}, which is an S3 generic with a default method (plus a method for the obscure class \code{\link{loc.ar}}; see \bold{Note}). The normal way to call the default, would be to give it a character matrix of genotypes, as well as sample-specific and locus-specific data. You can also use the default to give you an "empty" \code{snpgeno} of the right size but without genocalls. That might be useful if and only if you are planning to call your own genotypes and want to create the \code{snpgeno} manually, in which case the code of \code{snpgeno.default} and \code{snpgeno.loc.ar} may be informative.

A \code{snpgeno} can be subsetted by samples and/or loci, using numeric, logical, or character indices; see \bold{Subsetting}.

\code{snpgeno} specific methods currently exist for: \code{rbind}, \code{cbind}, dollar and dollar-assign (direct access to attributes such as sample-specific covariates and locus-specific information), subset (see \bold{Subsetting}), (in)equality (\code{==} and \code{!=}), \code{as.character}, \code{\link{print}}, \code{\link{str}}.
\subsection{Genotype encoding and storage}{
(See also \code{\link{get_genotype_encoding}} for the same information, written slightly differently.)

The word "genotype" can be a bit ambiguous, so we here use "genocall" to describe an already-called genotype (which may be ambiguous, null, etc) specifically for \emph{one} sample and \emph{one} locus. Genocalls are stored internally in a matrix of mode \code{raw} (1 byte per entry), for efficiency; if you do \code{unclass(x)[1:5, 1:3]} you can see the guts. As the next paragraph explains, the raw values are automatically converted to characters for printing, assignment, and equality-testing (about the only test you can apply to a genocall). Thus, things like \code{x[1,1]=="AA"} or \code{x[2,3] <- "BBO" or }x[4,5]==x[17,5]' should work fine. To do anything more sophisticated, first call \code{as.character(x)} to convert the raw values into a character matrix; don't try to handle the raw elements directly.

Normally all the genocall-encoding stuff (the mapping between raws and characters) will be set up for you automatically, through the constructor call or via \code{\link{read_vcf2snpgeno}} etc. You should almost never need to worry about the raw values themselves- and it's dangerous to mess with them. But for the record: the encoding, which applies to \emph{all} loci, is via the \code{diplos} attribute, a character vector which you can obtain by \code{diplos(x)}. It is indexed by the raw values, starting at \code{01}. So, if your \code{diplos} attribute was \code{c("AB","AA","BB")}, then a raw value of \code{02} would denote an "AA" homozygote (and would print as such). This allows for flexible encodings according to the nature of the data, including null alleles (via deliberately ambiguous genocalls such as"AAO"- either single-null, or reference homozygote) and multallelic loci with up to 6 alleles (though multiallelics are not used by package \pkg{kinference} version 1.x).

The interpretation of the encoding is determined entirely by whatever code processes it; there's no intrinsic meaning attached to the character representations. Most \code{kinference} functions assume some specific encoding, amongst one of those pre-defined by \code{\link{define_genotypes}} (qv), and they will check. You can convert manually between encodings, if you are feeling brave; see \code{get_genotype_encodings} for an example.

Missing genocalls should be encoded as raw value \code{ff}. They (along with any raw values outside the range \code{1:length(diplos(x)))}, such as \code{as.raw(0)}) will print as \code{"?"} and are converted to NA by \code{as.character(x)}, but there's otherwise no guarantee they're properly supported. Missing genocalls are specifically \emph{not} allowed in our \code{kinference} applications, so I haven't made much effort.

You \emph{can} do some dodgy things, such as assigning non-character items directly to elements of \code{x}, but you almost certainly \emph{shouldn't}. I'm not going to list them. Just don't.
}
\subsection{Other data}{
Data about individual samples ("metadata" for geneticists; "covariates" for statisticians) is stored in the \code{info} attribute, a dataframe which you can access via e.g. \code{my_sngeno$info} (thanks to S3 magic, because \code{snpgeno} is internally \emph{not} a list; it's just a matrix, with extra attributes). There must be as many rows in \code{info} as there are in the main genocall matrix. You can add whatever sample-specific fields you like to \code{x$info}. As SUBSETTING explains, you can later nominate one field as "sample ID".

The "locinfo" attribute, which must have as many rows as there are columns in \code{x}, can be accessed via \code{x$locinfo}. One column \emph{must} be called "Locus". Again, you can add whatever locus-specific fields you like (including matrices) to \code{x$locinfo}.

\code{snpgeno} methods do their damnedest to ensure that \code{rownames} are NULL for \code{info}, \code{locinfo}, and any other data named in \code{<snpgeno>$subset_like_both} (see \bold{User defined extras}). This might not always work... nevertheless, \code{rownames} on dataframes are pretty disastrous IMO (and I suspect they're an early design decision much regretted by guRus) and if you \emph{try} to use them with \code{snpgeno} objects, you \emph{will} hit some kind of trouble. But, see next subsection for a solution...
}
\subsection{Subsetting}{
You can subscript via \code{x[i,]} or \code{x[,j]} or \code{x[i,j]}- but not as \code{x[single_subscript]}, nor via matrix-subscripting. As usual in R{}, the subscripts can be integers, logicals, missing, or character. If the \code{j} subscript (for loci) is character, then it should match elements of \code{x$locinfo$Locus}. To use characters as \code{i} subscripts (i.e. by sample), see \code{\link{with_rowid_field}}. However, note that many functions in \code{kinference} require an \code{Our_sample} field (which should certainly be unique) and will not pay attention to \code{\link{rowid_field}} even if it nominates a different field.

\preformatted{
  my_snpg <- with_rowid_field( mysnpg, 'UniqueSampleID')
  my_snpg[ c( 'Abner12', 'Zoe9'),] # assuming...
  # ... that those names appear in 'my_snpg$info$UniqueSampleID'
  \%\%#
}

 Subsetting gets applied automatically to the \code{x$info} and \code{x$locinfo} attributes, and also to any attributes whose names are given in \code{x$subset_like_both}. Subset-replacement also works as usual.
\subsection{User defined extras}{
Aside from things that belong in \code{x$info} or \code{x$locinfo}, you can add arbitrary extra attributes via the "$<-" operator, e.g. \code{x$extra <- stuff}. These will not be printed by \code{print(x)}, but of course you can extract and manipulate them yourself if you want.

Occasionally, you might want to add something that is per-genocall, so also has dimensions of sample * locus (and possibly extra dimensions afterwards, too): per-allele counts would be an example. In that case, just add its name to the character-vector attribute \code{subset_like_both} (which will be non-existent if there are none such). For example:

\preformatted{
  x$qualitee <- array( 0, c( nrow( x), ncol( x), 3))
  x$subset_like_both <- c( x$subset_like_both, 'qualitee')
}

Thereafter, \code{$qualitee} should behave sensibly when \code{x} is subsetted (and allegedly also with \code{rbind} and \code{cbind}, though I'm skeptical). Note that it's up to you to make sure the something really has the its first two dimensions correct; this is S3, there is no safety net. There are no \code{dimnames} built into these attributes, but you can add them yourself; if so, presumably the first two should always equal \code{x$info$Our_sample} and \code{x$locinfo$Locus}.

For extra attributes that are not part of \code{x$subset_like_both}: \code{cbind(x,y,...)} and \code{rbind(x,y,...)} should preserve their values from \code{x} ie the first argument, but will discard all other copies (in \code{y} and so on).
}
}
}
\usage{
snpgeno( x, ...) # generic
\method{snpgeno}{default}( x, diplos, n_samples, n_loci, info, locinfo,  allow_nonchar, ...) # S3 method for default
diplos( x)
diplos( x) <- value
}
\arguments{
\item{  x}{thing to convert. For the \emph{default} constructor, this would normally be a character matrix of genocalls, but NULL is also allowed for an empty result that you can fill in manually later; iff you set \code{allow_nonchar=TRUE}, you can also pass in \code{raw} or \code{integer} values too, which must not exceed \code{length( diplos)}.}
\item{  diplos, value}{encoding for the genocalls, usually one of those in \code{define_genotypes} (qv). It's a character vector defining the genocalls to which the raw elements in \code{unclass(x)} will correspond.}
\item{  n_samples, n_loci}{for default constructor. If not specified, these will probably be deduced either from the dimensions of \code{x}, or from \code{info} and \code{locinfo}. Specifying them explicitly can be useful if you merely want an "empty shell" \code{snpgeno}.}
\item{  info}{dataframe with subject-specific usually-non-genetic data (name, date, size, ...). Must contain a field \code{Our_sample} (unique identifier for that sample/"library"/replicate/...); some downstream functions may require other fields too. For the default constructor, a placeholder will be constructed if \code{info} is not supplied.}
\item{  locinfo}{dataframe with locus-specific data. Must include "Locus" (name/definition); some downstream functions may require other fields too. For the default constructor, a placeholder will be constructed if \code{locinfo} is not supplied.}
\item{  allow_nonchar}{whether to allow raw or integer genocalls as input to the default constructor, for possible "efficiency"- but then it's the user's responsibility to ensure they do correspond appropriately to \code{diplos}).}
\item{  ...}{Extra user-defined attributes for the \code{snpgeno} object (default \code{snpgeno} constructor only). For \code{str.snpgeno}, the dots have no effect but must exist Becos R{}.}
\subsection{Note}{For constructing from an existing \code{\link{loc.ar}} object, \code{x$locinfo} must already include \code{x$locinfo$pambig} (a 4-column matrix of ABCO allele frequency estimates) and \code{x$geno_amb} (provisional genocalls). Those are generated as part of CSIRO's CKMR-genocalling pipeline (not yet public), and hopefully are documented in there...
}
}
\value{A \code{snpgeno} object.
}
\seealso{\code{\link{NGS_count_ar}}, \code{\link{loc.ar}}, \code{\link{define_genotypes}}
}
\examples{
# Create a snpgeno from scratch, filled with garbage.
set.seed(1111)
library( mvbutils) # becoz I say so
## Genotyping encoding? Show (current & legacy) possibilities
get_genotype_encoding()#
# ... let's use 4-way genotype encoding (probably commonest)
genotypes4_ambig <- get_genotype_encoding()$genotypes4_ambig
## Random genotypes...
n_loci <- 5
n_samples <- 3
genomat <- matrix(
    rsample( n_samples * n_loci, genotypes4_ambig, replace=TRUE),
    n_samples, n_loci)
## Locus information:
locodat <- data.frame(
      # "Locus" must be present: unique ID strings
    Locus= sprintf( 'L\%i', 1:n_loci),
    lenseq= rsample( n_loci, 100:200),
    chromo= rsample( n_loci, 1:24, replace=TRUE),
    poschro= round( 1e7 * runif( n_loci))
  )
## Sample information (ie covariates AKA "metadata")
sampodat <- data.frame(
      # "Our_sample" must be present: unique ID string
    Our_sample= sprintf( "S\%i", 1:n_samples),
    Year= rsample( n_samples, 2001:2004, replace=TRUE),
    Weight= runif( n_samples, 1, 5)
  )
## Put theem together
snpgarbage <- snpgeno(
    x = genomat,
    diplos = genotypes4_ambig,
    info = sampodat,
    locinfo = locodat
  )
snpgarbage
diplos( snpgarbage) # what encoding?
snpgarbage$info     # sample info
snpgarbage$locinfo  # locus info
str( snpgarbage) # ... confirming @reliability is there
# Subsetting:
mini <- snpgarbage[ 1:2, 1:2] # subset; also deals with $locinfo & $info
mini
mini$info
mini$locinfo
# Fix a "mistake" (manual editing):
snpgarbage[ 1, 1] <- 'OO'
# Some protection against users:
snpgarbage[ 1, 1] <- 'womble'
snpgarbage
table( as.character( snpgarbage))
# Other aspects of each genocall (a "user-defined extra"):
library( atease) # for x@y instead of attr( x, "y") etc
# which is very convenient
snpgarbage@manual <-
    matrix( FALSE, n_samples, n_loci)
snpgarbage@manual[1,1] <- TRUE # we fixed that one...
snpgarbage@subset_like_both <- 'manual'
snpgarbage    # @manul is not printed, but...
snpgarbage@manual # ... it is there
snpgarbage[ 1:2, 1:2]@manual     # ... and it subsets nicely
str( snpgarbage) # ... confirming @manual is there
}
\keyword{misc}
