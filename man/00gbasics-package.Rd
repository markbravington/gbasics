\name{gbasics-package}
\alias{gbasics}
\alias{gbasics-package}
\alias{gbasics}
\docType{package}
\title{Genotype storage ultimately for Close-Kin Mark-Recapture}
\description{
The ultimate purpose is to support Close-Kin Mark-Recapture, and specifically for the finding of close-kin pairs, which is handled by the \pkg{kinference} package. The only thing most (non-CSIRO) users will need is the \code{\link{snpgeno}} class, which stores called genotypes, sample-level covariates, and locus-level metadata. \code{\link{snpgeno}} is the only class accepted by functions in the \pkg{kinference} package.

The \code{\link{snpgeno}} object requires all genotypes to have been called already (that's not my problem!). It is normally created by reading from a SNPGDS or VCF file, or by converting from an existing R{} object of class \code{snpgds} (see package \pkg{SNPRelate}). But you can also construct a \code{\link{snpgeno}} object manually; see \code{\link{snpgeno}} examples.
\subsection{Legacy}{
\code{gbasics} also includes "legacy" support for other types of data and intermediate stages of genotyping, used in several CSIRO pipelines. In fact, there are 4 (!) S3 classes for handling multilocus genotype data (alongside sample-specific and marker-specific data), primarily NGS but also usats at a pinch. There are also a few utilities to help with SaddlePoint Approximation.

If you are starting from "raw genetic data" and using your own genotype-calling (genocalling) process, you may need the other classes as intermediate steps. CSIRO usually does genocalling for CKMR via package \pkg{genocalldart}, which is not publicly available. The other three \code{gbasics} classes are \code{\link{diploido}}, \code{\link{loc.ar}}, and \code{\link{NGS_count_ar}}. They are the results of evolution, specifically of our data-sources and genotyping processes; so they do have a few quirks (i.e. starting from what we know now, we would designed those classes differently), and there might be changes to their APIs in future. \code{\link{NGS_count_ar}} is perhaps the most stable.
}
\subsection{Null alleles and missings and dropout and encodings}{
\code{gbasics} doesn't do anything with your data, it just stores it in a useful format. Therefore, this isn't really the place to go into details about what your data "should" be, because \code{gbasics} itself doesn't care. However, the \pkg{kinference} package is much stricter, and it might help to be aware of some of this at the reading-it-in stage.

Basically, \code{kinference} expects a \emph{definite} genotype call at every locus for every sample. Of course, the call might be wrong, but that's survivable; the main thing is that it should not be missing or "dunno". \code{\link{snpgeno}} object \emph{can} store missing genotype calls, but you will have to impute them or something before you can use the object in \code{kinference}.

\code{kinference} is designed to cope with genuine null alleles, which are common in some loci in some CKMR applications with some genotyping methods (eg based on ddRAD). A null allele is a repeatable (ie across multiple samples from the same animal), heritable mutation at a locus, that is \emph{different} to the major & minor (reference & alternate) alleles and that doesn't show up in when genotyping in any directly-visible way. Nulls are \emph{completely different} from dropout where an allele fails to show up in one sample due to bad DNA or processing glitches. Occasionally, both copies of a locus will be null in one sample; this is a "double null" which should be clear from the genotyping details. More commonly, just one copy is null, so that a naive user might interpret the genotype as a homozygote for whichever (visible) allele is in the other copy. That would be a huge error for kin-finding purposes, even if it doesn't much matter in other genetic applications. Fortunately, the frequency of null alleles can be estimated basically from the excess (if any) of apparent homozygotes (plus any double nulls that are present), and nulls can then be allowed for in kin-finding. The \pkg{kinference} package has details, but

One common phenomenon is for genotyping software to put "empty" or "dunno" calls in when it's not sure. "Dunno" calls are \emph{not} the same as double-null calls; for the latter, the genotyping software should be confident that neither the Major nor the Minor allele is present. If your genotyping software is making many "dunno" calls, then you basically need to tell it to pull itself together and not be so wussy--- just make the damned call! Within limits, the \pkg{kinference} package is robust to some amount of genotyping error, but it won't accept missing data (a very deliberate decision, and a good one!).

To continue the rant: sometimes people use "null" loosely to mean "double-null" call--- so that "there are very few nulls in my dataset" might really mean "there are very few \emph{double} nulls--- or missings". Double-nulls aren't the main problem; \emph{single} nulls are much commoner, not directly visible, and can have a major effect on kin-finding.

Different genotyping methods might handle nulls differently; eg some methods \emph{can} infer the presence of a single null, albeit with uncertainty. Some methods lead to many more nulls than other methods. Also, it's not just about nulls: some methods and datasets might conceivably have a lot of loci with more than two non-null alleles (on the way to "microhaplotypes"), which are potentially powerful for kin-finding. To allow flexibility, the \code{\link{snpgeno}} class can therefore accept different sets of possible genotype-calls, as defined by the "genotype encoding". You will have to decide which encoding to use for your data; see \code{\link{snpgeno}} and \code{\link{get_genotype_encoding}} for more information. Only a few encodings are accepted by \code{kinference}, but it might be possible to convert a non-standard encoding manually before passing it to \code{kinference}; see \code{\link{get_genotype_encoding}}.
}
\subsection{Not currently true}{Package \pkg{gbasics} also includes some utility routines which duplicate stuff in other MVB-written packages, and which are included here to avoid excessive dependencies. These are being tidied up...
}
}
\keyword{misc}
