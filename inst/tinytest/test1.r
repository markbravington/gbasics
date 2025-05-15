## tinytest

library( gbasics)
library( mvbutils) # rsample etc, unless they move...
library( tinytest) # so I can check interactively using 'mdrun()'

set.seed( 999)
gencoding <- c( 'AB', 'AAO', 'BBO', 'OO')
genochars <- matrix( rsample( 12, gencoding, replace=TRUE), 3, 4)

sinfo <- read.table( header=TRUE, textConnection( r"--{
   Our_sample Our_plate
1  Billy       Wedgewood
2  Silly          Silver
3  Milly       Wedgewood
}--"))

linfo <- read.table( header=TRUE, row.names=NULL, textConnection( r"--{
Locus FullSeq Varpos Ref Alt
Sonic  GA.AC       3  T   C
Manic  TTT.        4  G   A
Panic  A.TTG       2  T   A
Tonic  GGGGG.      6  G   C
}--"))

# Manual construction
sg <- snpgeno(
    NULL,
    n_samples= nrow( sinfo),
    n_loci= nrow( linfo),
    diplos= c( 'AB', 'AAO', 'BBO', 'OO')
  )

sg
sg$info <- sinfo
sg$locinfo <- linfo
sg[] <- genochars
sg

# Direct
sg2 <- snpgeno(
    x= genochars,
    diplos= c( 'AB', 'AAO', 'BBO', 'OO'),
    info= sinfo,
    locinfo= linfo
)
expect_equal( sg, sg2)

# Cheaty non-numeric--- not recommended
rawg <- as.raw( match( genochars, diplos( sg2), 0L))
dim( rawg) <- dim( genochars)
expect_error(
  sg3 <- snpgeno(
      x= rawg,
      diplos= c( 'AB', 'AAO', 'BBO', 'OO'),
      info= sinfo,
      locinfo= linfo
    ) # splat!
  )
sg3 <- snpgeno(
    x= rawg,
    diplos= c( 'AB', 'AAO', 'BBO', 'OO'),
    info= sinfo,
    locinfo= linfo,
    allow_nonchar= TRUE
  ) # allowed
expect_equal( sg2, sg3)


# Basic subsetting
# I am not sure how best to formally test these... of course, I can test
# individual elements with 'expect_equal( as.character(...))' but conceivably
# the test should check the whole returned structure, which is a (small) 'snpgeno'
sg[ 1:2, 2:3]
sg[ 1,]
sg[ 0,0]
sg[1,1] <- 'nonsense!' # warning
sg # prints a "?"
sg[1,1] <- 'OO'

# Next gives an "illegal" (or pointless) result, but this is S3:
# ... it's up to you to not do stupid things!
sg[ c(1,1),] # duplicated

# rbind and cbind
rbind( sg, sg) # warning: samples
cbind( sg, sg) # warning: loci

sg2 <- sg
sg2$info$Our_sample <- 'B' %&% sg$info$Our_sample

rbind( sg, sg2)

sg2$info$Our_sample <- sg$info$Our_sample # else can't cbind; different samples!
sg2$locinfo$Locus <- 'B' %&% sg$locinfo$Locus
cbind( sg, sg2)

# Other subsettable attributes
sg$beeble <- matrix( 1:12, 3, 4)
sg$subset_like_both <- 'beeble'

# ... with optional dimnames
dimnames( sg$beeble) <- list( sg$info$Our_sample, sg$locinfo$Locus)
sg$beeble
sg$beeble[1,2] # drop
sg[1,2]$beeble # no drop; IMO that's fine!

# Let's try an array too. 3rd (and more) dimensions should be preserved
sg$fake_counts <- array( 1:24, c( 3, 4, 2)) # ie 2 values per genocall
sg$subset_like_both <- c( sg$subset_like_both, 'fake_counts')

sg[1:2,]$fake_counts
sg$fake_counts[1:2,,]
sg[0,0]$fake_counts

dimnames( sg$fake_counts) <- list( sg$info$Our_sample, sg$locinfo$Locus, NULL)

sg2 <- sg
sg2$locinfo$Locus <- 'B' %&% sg$locinfo$Locus
sg2$fake_counts <- sg2$fake_counts + 100
dimnames( sg2$fake_counts) <- list( sg2$info$Our_sample, sg2$locinfo$Locus, NULL)
testy <- cbind( sg, sg2)
testy
testy$fake_counts

sg2 <- sg
sg2$fake_counts <- sg2$fake_counts + 100
sg2$info$Our_sample <- 'B' %&% sg$info$Our_sample
dimnames( sg2$fake_counts) <- list( sg2$info$Our_sample, sg2$locinfo$Locus, NULL)
testy <- rbind( sg, sg2)
testy
testy$fake_counts

#... luvvly jubbly

