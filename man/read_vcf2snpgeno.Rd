\name{read_vcf2snpgeno}
\alias{read_vcf2snpgeno}
\alias{read_vcf2NGS_count_ar}
\title{Read VCF genotype data}
\description{
\code{read_vcf2snpgeno} and \code{read_vcf2NGS_count_ar} try to read the genotype part of a VCF file, and load it into a \code{\link{gbasics}} object: either a \code{\link{snpgeno}} (genotype calls) or a \code{\link{NGS_count_ar}} (sequence counts). They discard the "metadata" at the top of the file, and the "data lines" describing parts of the genome.

The \code{\link{snpgeno}} version expects a "GT" subfield containing genotypes. The \code{\link{NGS_count_ar}} expects a subfield "AD" (Allele Depth).

VCF format is complicated, and this isn't fully tested; the goal is to produce something that will get through the next few steps of the kinference process. If you want something more sophisticated, please feel free to hack this code (and if you do a good job of it, please let us know!). There are other R{} packages out there which- to \emph{some} extent- handle VCF, and you may be able to use one of those to create a \code{snpgds}, which can then be converted by \code{sngeno} (to the latter class).
}
\usage{
# BLOCK has same default in both functions
read_vcf2snpgeno(
  vfilename, BLOCK=formals( read_vcf2NGS_count_ar)$BLOCK,
  allow_disordered_unphased= TRUE)
read_vcf2NGS_count_ar(
  vfilename, biforce = TRUE, BLOCK = 1000,
    allow_disordered_unphased= TRUE, majjik_bleeble = "")
}
\arguments{
\item{  vfilename}{Filename, or an existing connection that is already open for reading. If \code{vfilename} is a connection, it is \emph{not} closed on exit.}
\item{  biforce}{if TRUE, condense all multi-allelic loci into just one Ref (as defined in the VCF- it's the first count) and one "Alt", the latter consisting of \emph{all} other sequence counts for that locus.}
\item{  BLOCK}{how many lines to read in at once. Default should be fine.}
\item{  allow_disordered_unphased}{By default, unphased genotypes are allowed in any order, so "1/0" and "0/1" are both tolerated and they mean the same thing. I'm not sure that's actually intended according to VCF specifications, though; my initial reading was that only non-decreasing (e.g. "0/1") was allowed, but I've now relaxed it to save having to write too many explanations like this one. Set the parameter to FALSE to enforce a \code{stop()} if any decreasing pairs are found.}
\item{  majjik_bleeble}{Do not mess with this.}
}
\value{
\code{\link{NGS_count_ar}} (qv) or \code{\link{snpgeno}} object. If the latter, then the genotype encoding (the \code{\link{diplos}} attribute) is currently \code{genotypes4_ambig}, ie AB/AAO/BBO/OO. The \code{info} attribute contains the following fields:
\item{  Our_sample}{column name from VCF}
\item{  Fishtot}{total counts for that sample, or NA if \code{snpgeno}}
\item{  File}{filename, or \code{description} of the connection}
\item{  MD5}{\code{md5sum()} of the file, or NA if a connection NB:"File" and "MD5" are character vectors, but also with class \code{dull} from package \pkg{mvbutils}, so they don't clutter the printout. It's only printing that's affected; and you can see their contents via \code{unclass(x$info$File)}. The \code{locinfo} attribute always contains the following fields:}
\item{  Locus}{a name for the locus, concocted from CHROM and POS in the VCF}
\item{  CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO}{as per VCF (these are mandatory fields) and for \code{NGS_count_ar} objects, \code{locinfo} will also contain these, which are needed by some functions that operate on \code{NGS_count_ar} objects:}
\item{  consensus}{same as ID}
\item{  var_pos}{info on variant position}
\item{  n_alleles}{what it says There are a couple of other housekeeping fields required by \code{NGS_count_ar} (qv), which need not concern us here. For \code{NGS_count_ar} objects, there is also a \code{seqinfo} attribute which, aside from housekeeping fields, contains the field "FullAlleleSeq", which is concocted from the ID and REF and ALT fields. I haven't tried to disentangle ALT (it doesn't look easy!), so if there's more than two allele (variant) at a locus, then the multiple ALTs are distinguished by ".1", ".2", etc after the \emph{entire} ALT field in the VCF. If the locus is biallelic, there's no ".1".}
}
\seealso{\code{\link{snpgeno}}, \code{\link{NGS_count_ar}}
}
\examples{
# package BinaryDosage has some nice small VCFs
# The one below just has genotypes
# But I'm not putting in Suggests; hassle of versions etc
# Workaround for CRANkiness:
re_bloody_quireNamespace <- get(
    sprintf( '\%s\%s', 're', 'quireNamespace'), baseenv()) # anti CRANky
if( re_bloody_quireNamespace( 'BinaryDosage'))\{
  thrubb <- read_vcf2snpgeno( system.file(
      'extdata/set1b.vcf', package='BinaryDosage'))
  print( thrubb)
\}
# Should also have a small public example file for 'read_vcf2NGS_count_ar' but
# must have AD subfield
# I don't have a good one
}
\keyword{misc}
