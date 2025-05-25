# This is package gbasics 

"!=.snpgeno" <-
function(e1, e2) {
  if( e1 %is.a% 'snpgeno'){
    if( is.character( e2)){
      e2 <- as.raw( match( e2, e1@diplos, 0))
    }
    return(unclass(e1) != e2)
  } else {
    if( is.character( e1)) {
      e1 <- as.raw(match( e1, e1@diplos, 0))
    }
  }
  return( unclass(e2) != e1)
}


"$.loc.ar" <-
function( x, name) return( attr( x, 'info')[[ name]])


"$.NGS_count_ar" <-
function( x, name) attr( x, name)


"$.snpgeno" <-
function( x, name){
  y <- attr( x, name)
  if( name %in% x@subset_like_both){
    if( !is.null( ridf <- rowid_field( x))){
      rownames( y) <- x@info[[ rowid_field]]
    } else {
      row.names( y) <- NULL
    }
    colnames( y) <- x$locinfo$Locus
  } else {
    # Strip any row names...
    if( (y %is.a% 'data.frame') &&
        (.row_names_info( y, type=1) >= 0)){
      row.names( y) <- NULL
    }
  }
return( y)
}


"$<-.loc.ar" <-
function( x, name, value) {
  attr( x, name) <- value
return( x)
}


"$<-.NGS_count_ar" <-
function( x, name, value) {
  attr( x, name) <- value
return( x)
}


"$<-.snpgeno" <-
function( x, name, value) {
  ## A bit of sanity-checking...
  if( name=='info'){
stopifnot( value %is.a% 'data.frame',
      nrow( value) == nrow( x)
    )
    ridf <- rowid_field( x)
    if( !is.null( ridf)){
stopifnot( ridf %in% names( value))
      value <- with_rowid_field( value, ridf)
    }
    if( .row_names_info( value, type=1) >= 0){
      row.names( value) <- NULL
    }
  } else if( name=='locinfo'){
stopifnot( value %is.a% 'data.frame',
      nrow( value) == ncol( x),
      'Locus' %in% names( value)
    )
    if( .row_names_info( value, type=1) >= 0){
      row.names( value) <- NULL
    }
  } else if( name %in% x@subset_like_both){
stopifnot( length( dim( value)) >= 2, # arrays OK
      nrow( value) == nrow( x),
      ncol( value) == ncol( x)
    )
    rownames( value) <- NULL
    colnames( value) <- NULL
  } # else whatever!

  attr( x, name) <- value
return( x)
}


".onLoad" <-
function( libname, pkgname){
  dedoc_namespace( pkgname) # drop plain-text docu (it's just for maintenance)
}


"[.diploido" <-
function( x, i, j, drop=FALSE) {
#scatn( 'NARGS=%i\n', nargs())
#scatn( 'NARGS-!missing(drop)=%i\n', nargs()-!missing(drop))
#scatn( 'missing(drop,i,j)=%i,%i,%i\n', missing( drop), missing( i), missing( j))
  if( (nargs()-!missing( drop))==2) {
    if( missing( j))
return( structure( unclass( x)[ i,,drop=FALSE], class='diploido'))
    else # missing i
      unclass( x)[,j,drop=TRUE]
  } else
    unclass( x)[ i, j, drop=drop]
}


"[.loc.ar" <-
function( x, i, j, k=1:nals, drop=FALSE){
  # x[i] means just some columns of info-- so i is cols of info
  # x[i,j] means i is rows and j is loci; k=1:nals in that case
  # x[i,j,k] means the works
  # Normally nals==2 but 'loc.ar' can also handle count-type data with
  # > 2 alleles (but fixed number per locus; otherwise, see 'NGS_count_ar')
  # which is a perversion really

  # Subscript checks taken from '[.data.frame'
  mdrop <- missing(drop)
  Narg <- nargs() - (!mdrop)

  info <- x@info
  loci <- named( x@locinfo$Locus) # named( dimnames( x)[[2]])
  nals <- dim(x)[3]

  # Other stuff, eg: print for specialprint; genotypes if x is really read-counts
  other_atts <- attributes( x) %without.name% cq( dim, dimnames)

  # ... requiring similar subsetting, which You can control by setting these attributes.
  subset_like_both <- names(other_atts) %that.are.in% other_atts$subset_like_both

  if( Narg >= 3) {
    oclass <- class( x)
    x <- unclass( x)
    if( !missing( j) && (length( loci[j])==1) && (length( (1:nals)[k])==nals) &&
        (missing( drop) || drop)) {
      # Extract single locus
      x <- x[ i, j, k, drop=FALSE]
      dim( x) <- dim( x)[-2]
      dimnames( x) <- list( NULL, loci[ j] %&% '.' %&% 1:2)

    } else {
      x <- x[ i, j, k, drop=FALSE]
      if( dim( x)[3]==2) { # still looks like a loc.ar
        dimnames( x) <- list( NULL, loci[j], NULL)
        oldClass( x) <- oclass
      } else if( drop && dim(x)[3]==1) {
        dim( x) <- dim( x)[1:2]
        dimnames( x) <- list( NULL, loci[j])
      }
    }
  }

  # Anything else to be subsetted the same way

  # This strange little function ensures that a 3D "like_loci" array will be
  # subscripted thus:
  # thing[ j, , , drop=drop] and a 2D "like_both" array as
  #  thing[ i, j, drop=drop] etc
  # Conceivably more efficient just to construct the string and then
  #  eval(parse(...)) but that feels like surrender...
  make_subset_call <- function( ld, ...) {
    dots <- match.call( expand.dots=FALSE)$...
    callo <- rep( list( formals( sys.function())$ld), ld+3)
    callo[[1]] <- quote( `[`)
    callo[[2]] <- quote( thing)
    callo[[ length( callo)]] <- quote( drop)
    names( callo) <- c( rep( '', length( callo)-1), 'drop')
    callo[ 2 + seq_along( dots)] <- dots
    as.call( callo)
  }

  # except it's not as "simple" as that. Classes with badly-written subset
  # methods, such as 'noquote' don't like recursive missing args (whereas eg
  # 'data.frame' is fine) Hence 'if( missing..)'  below, to insert "physical
  # missings". Sigh ...

  # ... except it's not as "simple" as that. Classes with badly-written subset methods, such as 'noquote'...
  # ... don't like recursive missing args (whereas eg 'data.frame' is fine)
  # ... Hence 'if( missing..)'  below, to insert "physical missings"
  # ... Sigh ...

  subset_like_loci <- 'locinfo'
  if( !missing(j)) {
    for( ioth in subset_like_loci) {
      # "j" will be first subscript, so eg thing[ j,,drop=FALSE] for locinfo
      thing <- other_atts[[ ioth]]
      other_atts[[ ioth]] <- eval( make_subset_call( length( dim( thing)), j))
    }
  }

  subset_like_samples <- 'info'
  if( !missing( i)) {
    for( ioth in subset_like_samples) {
      # "i" will be first subscript
      thing <- other_atts[[ ioth]]
      other_atts[[ ioth]] <- eval( make_subset_call( length( dim( thing)), i))
    }
  }

  if( !missing( i) || !missing( j)) {
    for( ioth in subset_like_both) {
      # "i,j" will be first two subscripts
      thing <- other_atts[[ ioth]]
      callo <- make_subset_call( length( dim( thing)), i, j)
      if( missing( i)) {
        callo[[3]] <- formals( sys.function())$i # a missing
      }
      if( missing( j)) {
        callo[[4]] <- formals( sys.function())$j # a missing
      }
      other_atts[[ ioth]] <- eval( callo)
    }
  }

  attributes( x) <- c( attributes( x), other_atts)
return( x)
}


"[.NGS_count_ar" <-
function( x, i, j, k, ...){
###########################################
# drop not allowed
  # Allow missing k (ie only 1 comma), but nothing else
  stopifnot( nargs() >= 3)
  no_k <- nargs() == 3

 # Dec 2020: *removed* character lookup for 'i' (and ability to lookup against arbitrary 'info' field)
  if( !missing( i)){
stopifnot( is.numeric( i) || is.logical( i))
  }

  if( !missing( j) && is.character( j)) {
    j <- match( j, x@locinfo$Locus, NA)
  }

  info <- x@info
  locinfo <- x@locinfo
  max_n_alleles <- max( locinfo$n_alleles)
  seqinfo <- x@seqinfo


  alli <- seq_len( nrow( x)) # dimnames( x)[[1]]
  loci <- named( locinfo$Locus) # named( dimnames( x)[[2]])

  # Other stuff, eg: print for specialprint; genotypes if x is really read-counts
  other_atts <- attributes( x) %without.name% cq( dim, dimnames, info, seqinfo, locinfo)

  x <- unclass( x)
  if( missing( i)) {
    i <- alli # 1 %upto% nrow( x)
  } else {
    i <- alli[ i] # ensure integer even if i was logical
  }

  if( !missing( j)) {
    if( !is.logical( j) && any( duplicated( j))) {
stop( 'Duplicated loci not allowed')
    }

    j <- structure( 1 %upto% nrow( locinfo), names=locinfo$Locus)[ j] # to integer

    # Need to allow for changes in locus ordering (sigh...)
    msa <- match( seqinfo$Locus, locinfo$Locus[ j], 0) # 222200111
    o <- order( msa)
    seq_allele <- o[ msa[ o] > 0]

    if( !my.all.equal( seq_allele, 1 %upto% nrow( seqinfo))) {
      x <- x[ , seq_allele, drop=FALSE]
      seqinfo <- seqinfo[ seq_allele,,drop=FALSE]
      locinfo <- locinfo[j,,drop=FALSE]
      locinfo$end_col <- cumsum( locinfo$n_alleles)
      # locinfo$start_col <- c( 1L, head( locinfo$end_col, -1)+1L)
    }
  }

  # do we REALLY need to do an i-subset?
  if( !my.all.equal( i, alli)){
    x <- x[ i,,drop=FALSE]

    # 'dull' attr on sub-data-frames will get mucked up
    # unfixable base-R thing
    # not sure I should encourage it here--- sub-data-frames are iffy
    # but...

    if( length( which_is_dull_df <- which( sapply( info, is.data.frame)))) {
      which_is_dull_df <- which_is_dull_df %SUCH.THAT% (info[[.]] %is.a% 'dull')
      if( length( which_is_dull_df)) { # make them temporarily interesting
        info[ which_is_dull_df] <- lapply( info[ which_is_dull_df], as.data.frame)
      }
    }
    info <- info[ i,,drop=FALSE]
    if( length( which_is_dull_df)) {
      info <- make_dull( info, names( info)[ which_is_dull_df])
    }

    # samps <- samps[ i]
  }

  # Under certain circumstances, we *might* subscript by k too
  # Returns a pure 3D array
  # k >= n_alleles for that locus => NA
  if( (!no_k || !missing( k)) && (length( k)>0)) {
    k <- (1 %upto% max_n_alleles)[ k] # logical to integer
# stopifnot( all( k) <= min( locinfo$n_alleles))
    start_col <- c( 0L, head( locinfo$end_col, -1))+1L
    xx <- do.call( 'cbind', FOR( 1 %upto% nrow( locinfo),
        x[ , (start_col[.] : locinfo$end_col[.])[k],drop=FALSE]))
    if( length( k) > 1) {
      dim( xx) <- c( nrow( x), length( k), nrow( locinfo))
      xx <- aperm( xx, c( 1, 3, 2))
      dimnames( xx) <- list( NULL, locinfo$Locus, NULL) # cannot assign names to the alleles, since they'll differ by locus
    } else { # save some work in the common case that length( k)==1
      dim( xx) <- c( nrow( x), nrow( locinfo))
      dimnames( xx) <- NULL # list( NULL, locinfo$Locus)
    }
return( xx)
  }

  x@locinfo <- locinfo
  x@seqinfo <- seqinfo
  x@info <- info
  attributes( x) <- c( attributes( x), other_atts)
return( x)
}


"[.snpgeno" <-
function( x, i, j, ..., drop=FALSE){
  # x[i] means just some columns of info-- so i is cols of info
  # x[i,j] means i is rows and j is loci; k=1:nals in that case
  # x[i,j,k] means the works
  # i defaults to using "Our_sample" column of x@info
  # but you can change that by naming a column, eg
  # x[ TargetID=c( 109679, 103258), 'L1049']
  # Normally nals==2 but 'loc.ar' can also handle count-type data with
  # > 2 alleles (but fixed number per locus; otherwise, see 'NGS_count_ar')
  # which is a perversion really

  # Subscript checks taken from '[.data.frame'
  # Other stuff, eg: print for specialprint; genotypes if x is really read-counts
  other_atts <- attributes( x) %without.name% cq( dim, dimnames)

  # ... requiring similar subsetting, which You can control by setting these attributes.
  # 2023: can only set 'subset_like_both'
  subset_like_loci <- 'locinfo' # names(other_atts) %that.are.in% c( other_atts$subset_like_loci, 'locinfo')
  subset_like_samples <- 'info' # <- names(other_atts) %that.are.in% c( other_atts$subset_like_samples, 'info')
  subset_like_both <- names(other_atts) %that.are.in% other_atts$subset_like_both

  oclass <- class( x)
  x <- unclass( x)

  rowid_field <- x@info@rowid_field # for character lookups
  if( !missing( i)) {
    if( is.character( i) && !is.null( rowid_field)){
      i <- match( i, x@info[[ rowid_field]], NA)
    }
stopifnot( (is.numeric( i) || is.logical( i)))
  }

  if( !missing( j) && is.character( j)) {
    j <- match( j, x@locinfo$Locus, NA)
  }

  x <- x[ i, j, drop=FALSE] # Amazingly, this works even if i / j are missing...
  if( drop && min( dim( x))==1) {
    dimx <- dim( x)
    x <- other_atts$diplos[ as.integer( x)]
    if( dimx[1]>1) {
      names( x) <- rownames( other_atts@info)[i]
    } else if( dimx[2]>1) {
      names( x) <- other_atts$locinfo$Locus[ j]
    } # if 1X1, don't know which you want...
return( x)
  }

  # Anything else to be subsetted the same way

  # This strange little function ensures that a 3D "like_loci" array will be subscripted thus:
  # thing[ j, , , drop=drop] and a 2D "like_both" array as thing[ i, j, drop=drop] etc
  # Conceivably more efficient just to construct the string and then
  # eval(parse(...)) but that feels like surrender...
  make_subset_call <- function( ld, ...) {
    dots <- match.call( expand.dots=FALSE)$...
    callo <- rep( list( formals( sys.function())$ld), ld+3)
    callo[[1]] <- quote( `[`)
    callo[[2]] <- quote( thing)
    callo[[ length( callo)]] <- quote( drop)
    names( callo) <- c( rep( '', length( callo)-1), 'drop')
    callo[ 2 + seq_along( dots)] <- dots
  as.call( callo)
  }

  # ... except it's not as "simple" as that. Classes with badly-written subset methods, such as 'noquote'...
  # ... don't like recursive missing args (whereas eg 'data.frame' is fine)
  # ... Hence 'if( missing..)'  below, to insert "physical missings"
  # ... Sigh ...

  if( !missing(j)) {
    for( ioth in subset_like_loci) {
      # "j" will be first subscript, so eg thing[ j,,drop=FALSE] for locinfo
      thing <- other_atts[[ ioth]]
      other_atts[[ ioth]] <- eval( make_subset_call( length( dim( thing)), j))
    }
  }

  if( !missing( i)) {
    for( ioth in subset_like_samples) {
      # "i" will be first subscript
      thing <- other_atts[[ ioth]]
      other_atts[[ ioth]] <- eval( make_subset_call( length( dim( thing)), i))
    }
  }

  if( !missing( i) || !missing( j)){
    for( ioth in subset_like_both) {
      # "i,j" will be first two subscripts
      thing <- other_atts[[ ioth]]
      callo <- make_subset_call( length( dim( thing)), i, j)
      if( missing( i)) {
        callo[[3]] <- formals( sys.function())$i # a missing
      }
      if( missing( j)) {
        callo[[4]] <- formals( sys.function())$j # a missing
      }
      other_atts[[ ioth]] <- eval( callo)
    }
  }

  attributes( x) <- c( attributes( x), other_atts)
return( x)
}


"[.with_rowid_field" <-
function( x, i, j, ... ){
  rowid_field <- attr( x, 'rowid_field')
  if( !missing( i) && is.character( i)){
    # Avoid recursion cos it's [[ ]] not single
    i <- match( i, x[[ rowid_field]], NA)
  }
  x <- NextMethod( '[')
  attr( x, 'rowid_field') <- rowid_field
x
}


"[<-.diploido" <-
function( x, i, j, value) {
# scatn( 'nargs,miss(i),miss(j): %i,%i,%i', nargs(), missing(i), missing( j))
  x <- unclass( x)
  if( nargs()==4)
    x[ i,] <- value
  else {
    mxi <- max( i)
    if( mxi > nrow( x))
      x <- rbind( x, matrix( NA, mxi-nrow( x), 2))
    x[ i, j] <- value
  }
  oldClass( x) <- 'diploido'
return( x)
}


"[<-.snpgeno" <-
function( x, i, j, value) {
# scatn( 'nargs,miss(i),miss(j): %i,%i,%i', nargs(), missing(i), missing( j))
  x <- unclass( x)
  if( is.character( value)) {
    value <- as.raw( match( value, x@diplos))
  }

  if( nargs()==3) { # just i; logical/matrix subscript, can't be missing
    x[ i] <- value
  } else if( nargs()==2 || (missing( i) && missing( j))) {
    x[] <- value
  } else if( missing( i)) { # i & j used but could be
    x[ , j] <- value
  } else { # i is present; use rowid_field lookup if required
    if( is.character( i)){
      i <- match( i, x@info[[ x@info@rowid_field]], NA)
    }
    if( missing( j)) {
      x[ i, ] <- value
    } else {
      if( is.character( j)){
        j <- match( j, x@locinfo$Locus, NA)
      }
      x[ i, j] <- value
    }
  }
  oldClass( x) <- 'snpgeno'
return( x)
}


"[<-.with_rowid_field" <-
function( x, i, j, value){
  if( !missing( i) && is.character(i)){
    i <- match( i, x[[ x@rowid_field]], NA)
  }
  NextMethod( '[<-')
}


"==.snpgeno" <-
function(e1, e2) {
  if(e1 %is.a% 'snpgeno'){
    if(is.character(e2)){
      e2 <- as.raw(match( e2, e1@diplos, 0))
    }
    return(unclass(e1)==e2)
  }else{
    if(is.character(e1)) {
      e1 <- as.raw(match(e1, e1@diplos, 0))
    }
  }
  return( unclass(e2)==e1)
}


"as.character.snpgeno" <-
function( x, ...) {
  y <- matrix( '', nrow( x), ncol( x))
  # Get character versions
  # pmax() cos 00 behaves differently and causes havoc;
  # whereas too big => NA, which is OK
  y[] <- x@diplos[ pmax( as.integer( x), 1L)]
  dimnames( y) <- list( rownames( x@info), x@locinfo$Locus)
return( y)
}


"as.data.frame.loc.ar" <-
function(x, row.names, optional, ...){
  la <- x
  df <- attr( la, 'info')
  rownames( df) <- NULL
  loci <- dimnames( la)[[2]]
  la <- unclass( la)
  apla <- aperm( la, c( 1,3,2))
  dim( apla) <- c( dim( la)[1], prod( dim(la)[2:3]))
  dimnames( apla) <- list( NULL, c( matrix( c( loci %&% '.1', loci %&% '.2'), nrow=2, byrow=TRUE)))
  df <- cbind( df, apla)
  df
}


"as.diploido" <-
function(x, ...) UseMethod( 'as.diploido')


"as.diploido.character" <-
function(x, ...){
  x[ is.na( x)] <- 'NA/NA'
  x[ !grepl( '/', x)] <- 'NA/NA'
  x <- suppressWarnings( matrix( as.integer( unlist( strsplit( x, '/'))), ncol=2, byrow=TRUE))
  class( x) <- 'diploido'
return( x)
}


"cbind.loc.ar" <-
function( ...){
  l <- list( ...)
  if( length( l) < 2){
return( l) # !
  }

  l1 <- l[[1]]
  info1 <- l1$info
  locinfo1 <- l1$locinfo

  att1 <- attributes( l1) %without.name% cq(
      info, locinfo,
      class, dim
    ) # these are the right elements; those in subset_like_both will get changed

stopifnot(
    all( do.on( l, . %is.a% 'loc.ar')),
    all( do.on( l, listmeq( .$info, info1))),
    # Next: shooould check at least classes of locinfo cols for consistency...
    # ... but this _is_ S3, after all!
    all( do.on( l, my.all.equal( sort( names( .$locinfo)), sort( names( locinfo1))))),
    # Next could certainly be more thorough, eg about object shape... just check existence for now
    all( do.on( l, my.all.equal( sort( .@subset_like_both), sort( l1@subset_like_both) ) ) )
  )

  locinfo1 <- do.call( 'rbind', FOR( l, .@locinfo))

  out <- guts_cbind( l, NULL)

  for( ioth in l1@subset_like_both) {
    att1[[ ioth]] <- guts_cbind( l, ioth)
  }


  if( any( duplicated( locinfo1$Locus))){
warning( "Duplicated 'Locus' field in cbind--- baaaad karma")
  }

  mostattributes( out) <- c(
      att1, list(
      dim= dim( out),
      info= info1, # destroyed by the raw cbind()
      locinfo= locinfo1, # changed
      subset_like_both= l1@subset_like_both # just the names
    ))

  oldClass( out) <- 'loc.ar'
return( out)
}


"cbind.NGS_count_ar" <-
function( ...){
  l <- list( ...)
  if( length( l) < 2){
return( l) # !
  }

  # This is a bit simpler than for loc.ar and snpgeno
  # cos no 'subset_like_both' att for NGS_count_ar (yet)
  l1 <- l[[1]]
  info1 <- l1$info
  locinfo1 <- l1$locinfo
  seqinfo1 <- l1$seqinfo

  # Get unchanging attributes into a list
  att1 <- attributes( l1) %without.name% cq(
      info, locinfo, seqinfo,
      class, dim
    )

stopifnot(
    all( do.on( l, . %is.a% 'NGS_count_ar')),
    all( do.on( l, listmeq( .$info, info1))),
    # Next: shooould check at least classes of locinfo cols for consistency...
    # ... but this _is_ S3, after all!
    all( do.on( l, my.all.equal( sort( names( .$locinfo)), sort( names( locinfo1))))),
    all( do.on( l, my.all.equal( sort( names( .$seqinfo)), sort( names( seqinfo1)))))
  )

  # Appended loci will correspond to different seqinfo cols in the result
  cols_so_far <- 0L
  for( i in seq_along( l)){
    l[[i]]$locinfo$start_col <- NULL
    l[[i]]$locinfo$end_col <- l[[i]]$locinfo$end_col + cols_so_far
    cols_so_far <- cols_so_far + nrow( l[[i]]$seqinfo)
  }

  locinfo1 <- do.call( 'rbind', FOR( l, .$locinfo))
  seqinfo1 <- do.call( 'rbind', FOR( l, .$seqinfo))

  out <- do.call( 'cbind', FOR( l, unclass( .)))
  mostattributes( out) <- c( list(
      dim= dim( out),
      info= info1,
      locinfo= locinfo1,
      seqinfo= seqinfo1),
      att1)
  oldClass( out) <- 'NGS_count_ar'
return( out)
}


"cbind.snpgeno" <-
function( ...){
  l <- list( ...)
  if( length( l) < 2){
return( l) # !
  }

  l1 <- l[[1]]
  ridf1 <- rowid_field( l1)
  diplos1 <- l1@diplos
  info1 <- l1$info
  locinfo1 <- l1$locinfo

  # subset_like_both is (if it exists) a charvec naming attributes
  # Get those actual attributes into a list
  att1 <- attributes( l1) %without.name% cq(
      info, locinfo,
      class, dim
    ) # these are the right elements; those in subset_like_both will get changed

stopifnot(
    all( do.on( l, . %is.a% 'snpgeno')),
    all( do.on( l, my.all.equal( .@diplos, diplos1))),
    all( do.on( l, listmeq( .$info, info1))),
    # Next: shooould check at least classes of locinfo cols for consistency...
    # ... but this _is_ S3, after all!
    all( do.on( l, my.all.equal( sort( names( .$locinfo)), sort( names( locinfo1))))),
    # Next could certainly be more thorough, eg about object shape... just check existence for now
    all( do.on( l, my.all.equal( sort( .@subset_like_both), sort( l1@subset_like_both) ) ) )
  )

  locinfo1 <- do.call( 'rbind', FOR( l, .@locinfo))

  out <- guts_cbind( l, NULL)
  for( ioth in l1@subset_like_both) {
    att1[[ ioth]] <- guts_cbind( l, ioth)
  }

  if( any( duplicated( locinfo1$Locus))){
warning( "Duplicated 'Locus' field in cbind--- baaaad karma")
  }

  mostattributes( out) <- c(
      att1,
      list(
      dim= dim( out),
      info= info1, # destroyed by the raw cbind()
      locinfo= locinfo1, # changed
      subset_like_both= l1@subset_like_both # just the names
    ))

  oldClass( out) <- 'snpgeno'

  if( !is.null( ridf1)){
    out <- with_rowid_field( out, ridf1)
  }
return( out)
}


"cbind.with_rowid_field" <-
function( ...){
  r1 <- rowid_field( list( ...)[[1]])
  y <- NextMethod( 'cbind')
  if( !is.null( r1)){
    y <- with_rowid_field( y, r1)
  }
return( y)
}


"define_genotypes" <-
structure( function( nlocal=sys.parent()) mlocal({
  ABCO <- named( cq( A, B, C, O))
  extract.named( ABCO) # A, B, C, and O are now defined
  NA_geno <- as.raw( 255)

  genotypes <- cq( OO, AO, BO, AB, AA, BB, AAO, BBO, AC, BC, CO, CC, CCO, BBOO) # maybe automate: sorted union of others?
  genotypes_ambig <- cq( OO, AB, AC, BC, AAO, BBO, CCO)
  genotypes4_ambig <- cq( OO, AB, AAO, BBO)
  genotypes6 <- cq( AA, AB, AO, BB, BO, OO)
  genotypes_C <- cq( AA, AB, AC, AO, BB, BC, BO, CC, CO, OO)
  genotypes3 <- cq( AB, AAO, BBOO)

  # Create objects called 'AAO' etc
  for( ig in genotypes) {
    # assign( ig, structure( as.raw( match( ig, genotypes)), class='ABOSNP'))
    assign( ig, structure( ig, class='noquote')) # for nicer printing
  }
})
, doc =  mvbutils::docattr( r"{
get_genotype_encoding       package:gbasics
define_genotypes


Diploid genotype encodings


DESCRIPTION

Many of the functions in 'kinference' etc can accept several different encodings for diploid genotypes, depending on how the genotyeps are distinguished. For example, one encoding specifies that single-nulls are called separately from homozygotes (perhaps with error), whereas another says that single-nulls and homozygotes are not distinguished; some encodings allow a 3rd or 4th etc option for the allele; etc.

Encoding is specified by the 'diplos' (qv) attribute in the 'snpgeno' object. The hope is that it should normally be handled for you automatically, but you can set the encoding manually to one of a few pre-specified options that 'kinference' will understand; the choices are shown by 'get_genotype_encoding()', though only a few will actually work with the 'kinference' package. An encoding is stored as a character vector of recognized genotype-categories; some are described next.


.SOME.POPULAR.ENCODINGS

This is probably the _wrong_ place to describe the various encodings currently used by 'kinference' (a few others are or were used by CSIRO's own genotyping pipeline before the results are 'kinference'-ready). Anyway, the commonest ones are:
 
 - 4-way ('genotypes4_ambig') allows the values "OO", "AB", "AAO", "BBO". This means nulls are allowed, double-nulls ("OO") are taken notice of, but single-nulls and homozygotes are not differentiated, so that "AAO" means "A" was present but "B" was not, and conversely for "BBO". You can use this encoding even if your data is guaranteed null-free; there just won't be any "OO", and the null-allele frequency can be set to 0, which means 'kinference' will always interpret "AAO" as a homozygote.
 
 - 3-way ('genotypes3') comprises "AB", "AAO", and "BBOO" where the latter covers  "BB" homozygotes and "BO" single-nulls _and_ "OO" double-nulls. The motivation is for loci where nulls are rare but not non-existent, and where double-nulls (which _might_ in practice be artefacts, although we assume they aren't) get wildly over-interpreted in a kinship setting. Merging genotype-categories never causes bias, but it does sacrifice a small amount of statistical information to gain (much more important) robustness. Nulls can actually be informative, so "blanket 3-way" is mostly not the default for CSIRO's datasets; problems only come when nulls are rare. So CSIRO has sometimes used 'genotypes4_ambig' for all loci, but setting the 'locinfo$useN' field to 3 for loci where double-nulls are too risky. Various functions in 'kinference' look at 'locinfo$useN' to decide whether an individual locus is treated _during analysis_ as 3-way, 4-way, or perhaps 6-way (next) rather than following the general encoding of the whole dataset.

 - 6-way ('genotypes6') with possibilities "AA", "AB", "AO", "BB", "BO", and "OO". Here, the genotyping process has tried to distinguish between single-nulls and true homozygotes, based on read-depth. The discrimination will never be 100% accurate, even for "good" loci, so explicit allowance must be made for getting it wrong sometimes (if you see references to 'snerr', that's why).  For some loci it's not even worth trying, so for them 'locinfo$useN' is set to 6 if it's worth trying to differentiate single-nulls from homozygotes for that locus, or 4 if nulls are fairly common but hard to differentiate, or 3 if nulls are rare and hard to differentiate. The pipeline for 6-way genotyping uses in-house CSIRO code and is quite fiddly--- so it's likely to _stay_ in-house! CSIRO has used 6-way genotyping for _large_ ongoing studies of SBTuna and school sharks, so this encoding will continue to be supported, but for new projects we would recommend 4-way instead, using more loci to compensate.
 
You will see that 4-way is a coarse-grained version of 6-way, and 3-way is a coarsened version of 4-way. The 'locinfo$useN' field (a column of numbers) can be used to selectively coarsen certain loci, but of course it can't "uncoarsen" any part of a 3-way-encoded dataset if the overall encoding was too coarse. In most cases, 4-way will be the starting point.


.INTERNAL.REPRESENTATION

Internally (though you _shouldn't_ need to know this), the genotypes in a 'snpgeno' are stored in a 'raw'-mode matrix, where each element takes 1 byte. The elements of the encoding correspond in order to the numbers 1,2,3,.... With 'genotypes4_ambig', for example, an element with raw value 1 means "OO", with value 2 means "AB", etc.  Value 0 should never be seen in real datasets, but is sometimes useful as a temporary placeholder to show that "something needs fixing"! It will display as "?".

Generally speaking, 'sngpeno' objects use raw value 'ff' (255) to indicate "missing" (as distinct from double-null). However, the 'kinference' package deliberately does not handle missing data; _all_ genotypes have to get called, even though some of those calls might be in error, and genotyping errors are allowed for (or ignored) in the subsequent statistical analyses.

Special S3 methods exist for printing, comparing, and assigning genotype encodings, so it should at least _look_ clear to you, even if you have no idea what all this is about. By design, you can't directly assign integer values into genotypes; see EXAMPLES.


.CHANGING.ENCODINGS

Mostly, encoding should be handled automatically during the creation of your 'snpgeno'. However, you do sometimes need do access to it, and _occasionally_ you might want to manually recode your data. That takes some care because the internal representation of the encodings won't generally match up; for example, "AB" is encoded as 1 for 'genotypes3' but as 2 in 'genotypes4_ambig'. See EXAMPLES for how to re-encode all genotypes from 4-way to 3-way (which you also could do less painfully via 'useN', without changing the encoding--- I think 'useN' is described more in the 'kinference-vignette'.


.EXPANDABILITY

You can in fact use any character vector you like; there is no need to stick to those provided by 'define_genotypes', and the interpretation of the character-strings is entirely at the discretion of downstream software. For example, you could specify your own encoding scheme for microhaplotypes; from memory, I think up to 7 variants could be accommodated within the 8-bit 'raw' storage. However, most of the 'kinference' functions can only handle one or two of the pre-specified encodings, so you'd need new software.


.INSIDE.CODE

From the command-line, and in many functions, 'get_genotype_encoding' is the most useful function. However, inside some code 'define_genotypes' might be more useful. 
It is meant only for use inside other functions, typically at the start of the function; it creates objects in whatever environment it was called from. Once you've called it, you can refer directly to a genotype '"AB"' as just 'AB' etc (ie no quotes), and you can also refer to genotype-encodings that it knows about by name, eg 'genotypes4_ambig'.

If _you_ run 'define_genotypes()' from the R prompt, those things will be created in '.GlobalEnv'. So, most people shouldn't. There is usually no specific need to run 'define_genotypes' in your own code either, unless you're aiming to extend 'kinference'; the point of documenting it here, is to explain and illustrate how the genotype encodings work. For All Practical Purposes, 'get_genotype_encoding' should be enough.



USAGE

# Really meant to be used early inside another function
define_genotypes(nlocal = sys.parent())
get_genotype_encoding()


ARGUMENTS

  nlocal: the frame number, or environment, to create things in- see 'mlocal'. Leave this alone unless you _really_ know what you're doing...


VALUE

'get_genotype_encoding' returns a list of character vectors. Their names aren't "sacred"; the only thing that matters to 'kinference' is the character-vector of strings in the encoding.

'define_genotypes' returns nothing, but various objects are created; see the code.


SEE.ALSO

'snpgeno'; 'cq' and 'mlocal' in package 'mvbutils'


EXAMPLES

# Use of define_genotypes() inside a function:
library( mvbutils) # just for '%is.a%'
count_hetz <- function( x){
    stopifnot( x %is.a% 'snpgeno')
    define_genotypes() 
    # Now eg genotypes4_ambig is directly available...    
    # and so is AB instead of "AB", etc. Cor blimey.
  return( sum( x==AB)) # thus saving a pair of quotes...
  }

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

}")

)

"dim.diploido" <-
function( x) NULL


"dim.NGS_count_ar" <-
function( x) c( nrow( unclass( x)), nrow( x@locinfo))


"dimnames.diploido" <-
function( x) NULL


"diploido" <-
structure( function( cop1, cop2){
  mm <- matrix( c( cop1, cop2), ncol=2)
  check <- mm==0
  mm[ check] <- NA
  check <- is.na( mm[,1]) & !is.na( mm[,2])
  mm[ check, 1] <- mm[ check, 2]
  check <- is.na( mm[,2]) & !is.na( mm[,1])
  mm[ check, 2] <- mm[ check, 1]
  check <- !is.na( mm[,1]) & (mm[,1] > mm[,2])
  mm[ check,] <- mm[ check,2:1]
  storage.mode( mm) <- 'integer'

return( structure( mm, class='diploido'))
}
, doc =  mvbutils::docattr( r"{
diploido      package:gbasics

Class for one diploid locus


DESCRIPTION

The ancient and probably obsolete class 'diploido' lets a diploid locus be treated as a 1D object (a vector), even though it's really a 2D matrix; no, I can't remember exactly why. The usual generics are available. Contrast with 'loc.ar', which works with entire data.frames most of which consists of a multilocus genotype.


USAGE

diploido(cop1, cop2)


ARGUMENTS

  cop1: integer vector

  cop2: integer vector


SEE.ALSO

loc.ar
}")

)

"diplos" <-
function( x) x@diplos


"diplos<-" <-
function( x, value){
  x@diplos <- value
return( x)
}


"find.root" <-
function(f, start, step, fdirection = ifelse(f(start + step, ...) > f0, "increasing", "decreasing"), target = 0,
    min.x =  -Inf, max.x = Inf, args.to.uniroot = list(), ...)
{
  f0 <- f(start, ...)
  if( f0==target)
return( start)

  step <- abs(step) * ifelse(xor(f0 < target, fdirection == "decreasing"), 1, -1)
  bound <- ifelse(xor.thing <- (step < 0), min.x, max.x)
  repeat {
    new <- start + step
    if(xor(new > bound, xor.thing))
      new <- (start + bound)/2
    f1 <- f(new, ...)
    if(xor(f0 < target, f1 < target))
      break
    start <- new
    f0 <- f1
    step <- step * 2
  }
  o <- order(x <- c(start, new))
  fvals <- c(f0, f1) - target
  dotnames <- paste(c(names(f)[1], names(list(...))), collapse = ",")
  ff <- function(...) f(...) - target

#  ff <- c(f[ - length(f)], list(FFFF = 1, target = 1), parse(text = "FFFF(" %&% dotnames %&% ") - target"))
#  mode(ff) <- "function"

  do.call("uniroot", c(list(ff, x[o]), # f.lower = fvals[o[1]], f.upper = fvals[o[2]]),
      args.to.uniroot, list(...)))$root
}


"format.diploido" <-
function( x, diploido.NA='NA/NA', ...){
  namx <- rownames( x)
  x <- unclass( x)
  nax <- is.na( x[,1])
  x <- x[,1] %&% '/' %&% x[,2]
  x[ nax] <- diploido.NA
  if( !is.null( namx))
    names( x) <- namx
return( x)
}


"get_genotype_encoding" <-
function(){
  define_genotypes()
return( mget( ls( environment(), pattern='^genotypes')))  
}


"guts_cbind" <-
function( l, namatt=NULL){
## cbind (ie along 2nd dimension) a bunch of things in a loc.ar or snpgeno
## the thing itself, or a specific attribute

  thing <- function( li){
      if( is.null( namatt)){
    return( unclass( li))
      } else {
    return( attr( li, namatt))
      }
    }

  a <- thing( l[[ 1]])
  if( length( dim( a))==2){
    newa <- do.call( 'cbind', lapply( l, thing))
    # ... OK for data.frames and matrices
    row.names( newa) <- NULL
  } else { # >= 3Darray. 'rbind' the product of first two dims, then reset dim
    # NB this cannot happen when namatt==NULL cos the main object is 2D
    dima <- dim( a)
    diman <- dimnames( a)
    newa <- do.call( 'rbind', FOR( l, {
        next_thing <- thing( .)
        pd12 <- prod( dim( next_thing)[1:2]) # length of first two dims
        dim( next_thing) <- c( pd12, length( next_thing) %/% pd12)  # 2D result
        next_thing
      }))
    dim( newa) <- c( dima[1], length( newa) / prod( dima[-2]), dima[-(1:2)])
    dimnames( newa) <- c( # all must be lists
        diman[1],
        if( is.null( namatt))
          NULL
        else
          list( unlist( FOR( l, dimnames( attr( ., namatt))[2]))),
        diman[-(1:2)])
  } # if a
return( newa)
}


"guts_rbind" <-
function( l, namatt=NULL){
## rbind (ie along 1st dimension) a bunch of things in a loc.ar or snpgeno
## the thing itself, or a specific attribute

  thing <- function( li){
      if( is.null( namatt)){
    return( unclass( li))
      } else {
    return( attr( li, namatt))
      }
    }

  a <- thing( l[[ 1]])

  if( a %is.a% 'data.frame'){
    # Enforce the more-sensible mvbutils version of rbind.data.frame
    newa <- do.call( 'rbdf', lapply( l, thing))
    row.names( newa) <- NULL
  } else if( length( dim( a))==2){
    newa <- do.call( 'rbind', lapply( l, thing))
  } else { # array. Coerce to matrix, then
    dima <- dim( a)
    diman <- dimnames( a)
    newa <- do.call( 'rbind', FOR( l, {
        subatt <- thing( .)
        pd1 <- dim( subatt)[1]
        dim( subatt) <- c( pd1, length( subatt) %/% pd1)  # 2D result
      subatt
      }))

    dim( newa) <- c( dim( newa)[1], dima[-1])
    dimnames( newa) <- c( # all must be lists
        if( is.null( namatt))
          NULL
        else
          list( unlist( FOR( l, dimnames( attr( ., namatt))[1]))),
        diman[-1])
  } # if data.frame or whatever
return( newa)
}


"head.loc.ar" <-
function (x, n = 6, nl=5, ...) {
stopifnot(
    length(n) == 1L,
    length( nl)==1
  )

  nrx <- nrow( x)
  n <- if (n < 0)
      max(nrx + n, 0)
    else
      min(n, nrx)

  ncx <- ncol( x)
  nl <- if( nl < 0)
      max( ncx+nl, 0)
    else
      min( nl, ncx)
  x[seq_len(n),seq_len(nl),, drop = FALSE]
}


"head.NGS_count_ar" <-
function (x, n = 6, nl=5, ...) {
stopifnot(
    length(n) == 1L,
    length( nl)==1
  )

  nrx <- nrow( x)
  n <- if (n < 0)
      max(nrx + n, 0)
    else
      min(n, nrx)

  ncx <- ncol( x)
  nl <- if( nl < 0)
      max( ncx+nl, 0)
    else
      min( nl, ncx)
  x[seq_len(n),seq_len(nl)]
}


"head.snpgeno" <-
function (x, n = 6, nl=5, ...) {
    stopifnot(length(n) == 1L)
    nrx <- nrow( x)

    n <- if (n < 0)
        max(nrx + n, 0)
      else
        min(n, nrx)

    ncx <- ncol( x)
    nl <- if( nl < 0)
        max( ncx+nl, 0)
      else
        min( nl, ncx)
    x[seq_len(n),seq_len(nl)]
}


"inv_CDF_SPA2" <-
function( p, K, dK, ddK, tol=formals( ridder)$tol, already_vectorized=TRUE) {

######## Invert L-R SPA approx to CDF on "s-scale"
######## Avoids "double iteration" of nonlinearity
######## Should work for vector p (the target) but only if your K etc do

  isqrt_2pi <- 1/sqrt(2*pi)
  x <- sqrt_ddK_s <- Leg_trans <- u <- w <- 0*p # vectorized

  if( !already_vectorized) {
    K <- Vectorize( K) # names( formals( K)) %except% names( list( ...)) kinda thing
    dK <- Vectorize( dK)
    ddK <- Vectorize( ddK)
  }

  CSPA <- function( s) {
      x <<- dK(s)
      Leg_trans <<- s*x - K(s)
      w <<- sign(s) * sqrt( 2*Leg_trans)
      sqrt_ddK_s <<- sqrt( ddK( s))
      u <<- s * sqrt_ddK_s
    return( pnorm( w) + dnorm( w) * (1/w - 1/u) - p_target)
    }

  p_target <- 0
  isK2 <- 1/sqrt(ddK(0))
  seps <- 0.001 * isK2 # fraction of 1 SD; (tol/2) tends to give numeric errors...
  p0 <- (CSPA( seps) + CSPA( -seps)) / 2 # CSPA(0)==NA--- avoid!
  is_lower <- p0 > p
  # K3 <- 3 * (sqr(u)-sqr(w)) / seps^3 # yeh not bad FWIW

  hi <- lo <- 0*p

  # Shifted start...
  q0 <- qnorm( inv.logit( logit( p) - logit( p0)))
  bingo <- q0==0.5 # this case won't converge
  mean_x <- dK( 0)

  # L-R does not work at x==E[X]..!
  if( any( bingo)) { # bingo!
    # This was in the original code, and needed vectorizing...
    # I think p[bingo] corresponds exactly (??) to mean(X)
    if( all( bingo)) {
return( 0*x + mean_x)
    }

    # change p for those cases to something that will converge and carry on
    # sub back correct x (ie mean) on exit
    p[ bingo] <- p[ !bingo][1]
    q0[ bingo] <- q0[ !bingo][1]
  }

  p_target <- p # so CSPA gives 0 at solution

  # Only risk I can see, is that if true s ~= 0, start may be on wrong side... and CSPA calcs go wrong.
  s0 <- 0.9 * q0 *isK2 # "zeroth-order" approx, times 0.9 for guess at lower bound
  while( any(
      bad <- !is.finite(
        C0 <- CSPA( s0)))) {
    s0[ bad] <- s0[ bad]/2
  }

  toobig <- C0 > 0
  multor <- 2 * xor( toobig, s0>0) - 1 # -1 or +1
  bracketed <- s0 != s0
  step <- 1.2 ^ multor
  snext <- s0

  repeat{
    snext[ !bracketed] <- s0[ !bracketed] * step[ !bracketed]
    bracketed[ !bracketed] <- xor( toobig, CSPA( snext) > 0)[ !bracketed]
    if( all( bracketed))
  break
    s0[ !bracketed] <- snext[ !bracketed]
  }

  hi <- pmax( s0, snext)
  lo <- pmin( s0, snext)

  s <- ridder( CSPA, lo, hi, tol=tol, skip_bounds=TRUE) # root finder
  CSPA( s)

  x[ bingo] <- mean_x # any that hit first time
return( x)
}


"is.array.diploido" <-
function( x) FALSE


"is.matrix.diploido" <-
function( x) FALSE


"is.na.diploido" <-
function( x) is.na( unclass( x)[,1])


"is.na<-.diploido" <-
function( x, value) {x[ value,1] <- x[ value,2] <- NA; x}


"length.diploido" <-
function( x) nrow( unclass( x))


"listmeq" <-
function( a, b){
  # Used in rbinds and cbinds to check that data.frames (or lists...) have the same contents
  # even tho order of cols is permutable for rbind.data.frame()
  my.all.equal( a[ order( names(a))], b[ order( names( b))])
}


"loc.ar" <-
structure( function( x, ...)
  UseMethod( 'loc.ar')
, doc =  mvbutils::docattr( r"{
loc.ar      package:gbasics
unloc.ar

Create & manipulate locus-arrays


DESCRIPTION

A 'loc.ar' object is meant for "numerical diploid genotype data". The class is _not_ used by the 'kinference' package, is not particularly well thought out, and may not be properly documented here. We no longer use 'loc.ar's much at CSIRO except as intermediate steps in constructing a 'snpgeno', but in the past we have used them for: (i) microsat genotypes, with each sample-locus being a pair of numbers for allele-length; and (ii) biallelic SNP read-depth pairs, with the genotypes (if decided) stored as a hidden extra. The latter purpose is generalized to "microhaplotypes" (multiple alleles per locus) by 'NGS_count_ar' (qv), which is probably more "future-proofed". 'loc.ar' probably won't go away from 'gbasics', but it is not guaranteed to be actively maintained (whereas 'snpgeno' and 'NGS_count_ar' probably will be); then again, perhaps it will be improved in future. Who knows?

There are methods for the usual things like subsetting (see *Details*), 'rbind', 'cbind', 'print', 'head', 'tail', 'str'. There is also a constructor 'loc.ar'- an S3 generic- but currently there are _no_ methods for it (at least not in package 'gbasics'; a dodgy one has been moved to 'genocalldart' for now). Any 'loc.ar' object has to be built manually, as done by a few functions in package 'genocalldart', such as 'read_count_dart'.

At heart, a 'loc.ar' is a 3D numeric (or integer) array with dimension '(n_samples,n_loci,X)' where X is a fixed maximum number of non-null alleles for each locus (always either 2 or 3 in our applications to date), but with attributes that contains auxiliary data: sample-specific, locus-specific, optionally sample-by-locus, plus whatever other fixed attributes you want to give it. It is broadly similar to 'snpgeno', which is better documented, except for the details of its contents (3D numeric, vs 2D 'raw'-mode interpreted as actual called genotypes). Currently, there are methods for printing, subsetting, rbind/cbind, and for transforming to/from dataframes with two cols per locus. Subsetting can handle the auxiliary data automatically, as described below.

Auxiliary data is stored in attributes which can be accessed via the '$' operator. Sample-specific auxiliary info (which must exist) can be accessed via 'x$info'; it should be a dataframe containing a column "Our_sample". Locus-specific auxiliary info is optional (but you'd be brave not to have it...); it should live in a dataframe containing one column "Locus" (the name), and is stored in 'x$locinfo'. Loci are _rows_ here, whereas loci are _columns_ in 'x' itself. Other than subsetting, there are currently no special operations on 'locinfo'; for example, 'print( x)' does not show it (whereas, for 'snpgeno' objects, the field "Locus" is shown vertically as the column name).

Further auxiliaries, such as SNP genotypes (sample-by-locus), can also be tacked onto the attributes as you please. All attributes are copied by subsetting. By default, the whole attribute is copied, except for (i) 'info' and 'locinfo', which subset automatically in the appropriate way, and (ii) any auxiliaries named in the attribute 'x$subset_like_both' (a character vector that you can set manually), which should be arrays/matrices whose first two dimensions pertain to sample and to locus respectively (just like the main data in a 'loc.ar'); they can have more dimensions, too.

I sometimes add class 'specialprint' from package 'mvbutils' to 'loc.ar' objects, to _suppress_ printing of some auxiliaries. See EXAMPLES.


.INTERNAL.NOTE

I used to use class 'diploido' to hold microsat genotypes inside dataframes, but that was less flexible. The 'genocalldart::SBTlike_loc.ar' creator function can take input from "diploido"-class objects in dataframes, rather than just dataframes with e.g. "D225_A.1" and "D225_A.2" columns. If any column is of class "diploido", all are assumed to be, and "normal" loci are ignored.


USAGE

loc.ar( x, ...) # generic
unloc.ar( la)


ARGUMENTS

  x: in theory, something to be converted into a 'loc.ar' (but there are no methods yet!)

  la: a 'loc.ar' with exactly 2 alleles per locus (thus an 'S*L*2' array) to be turned into a dataframe with the 'info' columns, plus two columns for each locus. Why, you may ask?
  
  ...: just in case


DETAILS

Sample info and locus info are stored as dataframes in 'attr(<x>,"info")' and 'attr(<x>,"locinfo")' respectively. Internally, the code uses the 'atease' package so it can just write e.g. 'x@info' instead, but you should not need to load 'atease' because you can always use the '$' and '$<-' accessors instead.

For subsetting, note that 'k', or 'j' _and_ 'k', can be omitted, in which cases 'i' and 'j', or 'i' alone, will be used as the _first_ subscript; this is _unlike_ dataframes, where a single subscript gets used as the _second_ subscript. If the resulting array contains only one locus, it will by default be collapsed to an n-by-2 matrix, unless 'drop=FALSE' explicitly. Otherwise, the result will stay as a 'loc.ar' provided the 3rd dimension stays at length 2. The first dimension (samples) is never dropped.

'unloc.ar' should probably acquire a "diploido.use" argument, to optionally return loci as 'diploido' objects rather than pairs of columns.


VALUE

  loc.ar, [, rbind: a 'loc.ar'

  unloc.ar: a dataframe with two cols per locus


SEE.ALSO

NGS_count_ar, snpgeno
}")

)

"make_genopairer" <-
structure( function( genotypes){
  n_geno <- length( genotypes)
  res <- matrix( 0L, n_geno, n_geno)
  res[] <- seq_along( res)
  swappo <- row( res) < col( res)
  res[ swappo] <- t( res)[ swappo]
  res[] <- match( res, sort( unique( c( res))))
  dimnames( res) <- rep( list( genotypes), 2)

  # What are the genopairs? These are NOT duplicated, eg...
  # ... AB/AA doesn't imply samp 1 is AB
  what <- outer( genotypes, genotypes, paste, sep='/')
  res@what <- c( what)[ match( 1 %upto% max( res), res)]
return( res)
}
, doc =  mvbutils::docattr( r"{
make_genopairer      package:gbasics

Mapping between pairs of possible genotypes and compressed form


DESCRIPTION

You almost certainly do not need, nor want, to call this. But it needs to be exported so that 'kinference' and other packages can find it.

With diploid genotypes such as "AA" and "AB" at a locus, the pairwise probs under given co-inheritance (0, 1, or 2 copies) are symmetric. It is daft to store the full symm matrix, and causes some boring technical problems too. So, 'make_genopairer' creates a mapping between a vector of pair-names such as "AA/AB" and the underlying pairs. It allows eg a 2000x6x6 array to be stored as 2000x15, and makes passing it to Rcpp faaaaar easier.


USAGE

make_genopairer(genotypes)


ARGUMENTS

  genotypes: character vector, produced by e.g. 'define_genotypes', which you probably shouldn't be calling either.


VALUE

Matrix with rows & columns being the individual genotypes, and entries a lookup from that pair into a character vector of possible _ordered_ pairs, which is stored as the "what" attribute. See *Examples*.


EXAMPLES

make_genopairer( c( 'AA', 'AB', 'BB'))
#   AA AB BB
#AA  1  2  3
#AB  2  4  5
#BB  3  5  6
#attr(,"what")
#[1] "AA/AA" "AB/AA" "BB/AA" "AB/AB" "BB/AB" "BB/BB"
}")

)

"matches.diploido" <-
function( x, y, sloppy=FALSE, ...){
  strict <- is.na( x[,1]) | is.na( y[,1]) | ((x[,1]==y[,1]) & (x[,2]==y[,2]))
  if( !sloppy)
return( strict)

  # Else use loose matching (AA ~ AB), but record these as imperfect
  loose <- is.na( x[,1]) | is.na( y[,1]) | ((x[,1]==y[,1]) & (x[,2]==y[,2])) |
      ((x[,1]==x[,2]) & ((y[,1]==x[,1]) | (y[,2]==x[,1]))) |
      ((y[,1]==y[,2]) & ((x[,1]==y[,1]) | (x[,2]==y[,1])))

  if( any( imperfect <- loose & !strict)) {
    synth <- x[ imperfect,,drop=FALSE]

    xy <- cbind( unclass( x), unclass( y))
    el1 <- apply( xy, 1, min)
    el2 <- apply( xy, 1, max)
    synth[ ,1] <- el1[ imperfect]
    synth[ ,2] <- el2[ imperfect]

    attr( loose, 'partial') <- index( imperfect)
    attr( loose, 'synthesis') <- synth
  }

return( loose)
}


"NGS_count_ar" <-
structure( function( x, ...) UseMethod( 'NGS_count_ar')
, doc =  mvbutils::docattr( r"{
NGS_count_ar      package:gbasics
NGS_count_ar.matrix
print.NGS_count_ar
print

Next-Gen sequencer count data


DESCRIPTION

'NGS_count_ar' is a class for raw counts by sequence-variant within locus, per sample. The number of alleles (variants) can differ between loci. It's not used by any functions in the 'kinference' package, which expects 'snpgeno' objects throughout. However, it's used extensively in CSIRO's private 'genocalldart' package. It's an early stage in processing DartSeq/DartCap/DartTag data--- long before genotypes are called--- and can also be used to hold similar data from other formats such as VCF, via 'read_vcf2NGS_count_ar()'. It's also a generic S3 method for creating such an object, with one low-level method for matrices (which you probably shouldn't call yourself) and another method for converting from 'loc.ar' objects, which will presumably have come from 'genocalldart::read_count_dart' or such. The default method (ie if 'NGS_count_ar' is called on any old rubbish) will 'stop()', by design. See *Details* for internal structure, and what you can add/expect.

There are methods for basic stuff: try 'methods(class="NGS_count_ar")'. 'dim(x)' returns '#samples*#loci'; for the number of alleles, use 'nrow( x$seqinfo)'. The 'print' method tries to make things clearly readable, and has a few extra arguments to help control that- see *Arguments*. Subsetting can have up to 3 indices 'i', 'j', and 'k'; 'j' can be numeric, logical, or character (i.e. locus names, which are matched against 'x$locinfo$Locus'). If 'k' is omitted, another 'NGS_count_ar' object is returned, as controlled by 'i' and/or 'j'. If 'length(k)==1', then the result is a 2D array of sample-by-locus. If 'length(k)>1', the result is a 3D array of sample-by-locus-by-allele, with NA counts added when 'k' refers to an allele "number" that doesn't exist for that locus; obviously, all loci then have the same number of "alleles" in the new structure. 'dimnames' are set for the 2nd (locus) dimension only. These "pure array" results don't preserve the detailed sample or locus information (except via the 'dimnames') but can be useful for subsequent manipulation.


.OBSCURE.NOTE

A similar class is 'loc.ar', which preserves counts but allows max 3 non-null alleles per locus. It's not clear which of 'loc.ar' or 'NGS_count_ar' is "more advanced"- depends on the context. However, most later stages in 'genocalldart' pipeline (eg for bait selection) require 'loc.ar'. To go back the other way, you probably need to call 'pick_ref_alt' (qv); see the pipeline examples.


USAGE

NGS_count_ar( x, ...) # generic
NGS_count_ar( x,  # S3 method for "matrix'
    sampinfo, locinfo, seqinfo,
    strip_numerics = TRUE,  rename = TRUE, ...) # S3 method for matrix
print( x,  # S3 method for "NGS_count_ar'
    trailing_dot = getOption("trailing_dot_NGS_count_ar", FALSE),
    dot_for_0 = getOption("dot_for_zero_NGS_count_ar", FALSE),
    center_dot=centre_dot,
    centre_dot= getOption( 'centre_dot', getOption( 'center_dot', FALSE)),
    ...) # S3 method for NGS_count_ar


ARGUMENTS

  x: for the constructor, an integer matrix with dimensions (allele,sample), i.e. the _transpose_ of final storage. Or a 'loc.ar'. For 'print', an 'NGS_count_ar'

  sampinfo: dataframe of sample info

  locinfo: dataframe of locus info

  seqinfo: dataframe of sequence info (including which locus the sequence belongs to)

  strip_numerics: if TRUE, remove any non-integer columns from 'locinfo'- typically, these are summary data in a CSV that you can live without

  rename: certain names are expected by fun/xctions that handle 'NGS_count_ar' objects; this patches up "alternative" names that are generated by 'read_cluster_dart3'

  trailing_dot: set TRUE if you want to show that count data is non-integer (eg after norming by sample-total-reads) so that all counts end with a period. The post-decimal-place digits are not normally important, but if you really need to see them, you can do so by setting a 'k' subscript. Or set the global option, as per *Usage*.

  dot_for_0: set TRUE if you want zero-counts replaced by a dot. Using the global option is often helpful. Or set the global option, as per *Usage*.

  center_dot, centre_dot: iff 'dot_for_0' is TRUE, whether to use a central dot (Latin-1 and Unicode 0xb7) or a period (if FALSE) to replace leading zeros.

  ...: passed to other methods


DETAILS

The allele counts are stored as a (sample*sequence) matrix. There are also three key attributes, all accessible using the '$' operator, and each a dataframe:

 info: Metadata on each sample, including "Our_sample" and, usually for Dart data at least, "TargetID"

 locinfo: Metadata on each locus; see below

 seqinfo: Metadata on each allele, including "Locus", "FullAlleleSeq", and probably others

The key fields in 'x$locinfo' (there might be others) are

 Locus: Unique strings

 n_alleles: Number of sequence variants (presumably 2 or more)

 consensus: String of the overall sequence, with dots at SNP sites

 var_pos: String showing the positions of SNP sites, eg "17,41"

 end_col: of the sequences found at that locus (needed for lookup into the main matrix)

You can also add other attributes that apply to the whole object. One that's used in 'genocalldart' is 'mean_fish_tot' (used in norming count data from new samples), but it's rather dicey for general use because it does depend on the set of loci; so if you subset by loci, things (should) change...


SEE.ALSO

'loc.ar', 'snpgeno'
}")

)

"NGS_count_ar.default" <-
function( x, ...) stop( "Needs 'matrix', 'loc_ar', or something like that...")


"NGS_count_ar.loc.ar" <-
function( x, ...){
## Code for converting from loc.ar to NGS_count_ar--- even though the former is ultimately
## more useful, and a back-conversion will probably be required (see 'pick_ref_alt')!

  guts <- unclass( x)
  loci <- x@locinfo
  sampi <- x@info
  attributes( guts) <- attributes( guts) %without.name% cq( locinfo, sampinfo)

  ## *Probably* need to reorder the contents with something like this...
  thrubbo <- aperm( guts, c( 1, 3, 2)) # a guess...
  dim( thrubbo) <- c( dim( thrubbo)[1], prod( dim(thrubbo)[2:3]) )

  ## Then disentangle original @locinfo into new @locinfo and @seqinfo and
  ## change the 'class' attribute (see also 'oldClass()').

  loci$end_col <- cumsum(loci$ClusterSize)
  loci$n_alleles <- c(rep(2, nrow(loci)))
  loci$var_pos <- as.character(loci$SnpPosition)
  loci$consensus <- loci$ClusterConsensusSeq

  ## replace the random characters with "." in the consensus seq
  for(i in 1:length(loci$consensus)) {
      seq <- loci$consensus[i]
      point <- as.integer(loci$var_pos[i])
      substring(seq, point + 1, point + 1) <- "X" ## weird OBOE-looking thing here...
      loci$consensus[i] <- seq
  }

  # Keep 'locinfo' lean & mean
  # Shane's code was:   loci <- data.frame(Locus = loci$Locus, n_alleles = loci$n_alleles, ...)
  loci <- loci[ cq( Locus, n_alleles, end_col, consensus, var_pos)]

  # Shane's code duplicated stuff between loci and seqi
  seqi <- data.frame(Locus = rep(loci$Locus, times = ClusterSize),
#                     ClusterSize = rep(ClusterSize, times = ClusterSize),
#                     Allele = rep(AlleleID, times = ClusterSize),
#                     consensus = rep(ClusterConsensusSeq, times = ClusterSize),
                     FullAlleleSeq = rep(NA, sum(loci$ClusterSize)))
  seqi$FullAlleleSeq[seq(1, nrow(seqi)-1, 2)] <- loci$RefAlleleSeq
  seqi$FullAlleleSeq[seq(2, nrow(seqi), 2)] <- loci$AltAlleleSeq
  seqi <- make_dull( seqi, 'FullAlleleSeq')

  seqi$count_sum <- colSums(thrubbo) # for pick_ref_alt

  thrubbo@info <- sampi
  thrubbo@seqinfo <- seqi
  thrubbo@locinfo <- locinfo

  # Used to say this:
  ## need dimnames for the subsetter to work properly...
  # ... but (Dec 2020) I don't think it's true, and...
  # ... Our_sample isn't guaranteed unique anyway (only TargetID is)
  # dimnames(thrubbo)[[1]] <- thrubbo@info$Our_sample
  # dimnames(thrubbo)[[2]] <- make.unique(thrubbo@seqinfo$Allele)

  class(thrubbo) <- "NGS_count_ar"

return(thrubbo)
}


"NGS_count_ar.matrix" <-
function( x, sampinfo, locinfo, seqinfo, strip_numerics=TRUE, rename=TRUE, ...) {
  # Construct 'NGS_count_ar' object
  rownames( x) <- NULL # simplify
  x <- t( x)
  # Used to have this:
  # dimnames( x) <- list( Well=wells[ is.fish.col], NULL)
  # but I don't think we need a dimname, and if we do it shouldn't be "Well" anyway...
  # ... seriously non-unique !

  class( x) <- 'NGS_count_ar'

  if( strip_numerics) {
    # Remove summary properties, eg call-rate
    locinfo <- locinfo %SUCH.THAT% (. %is.not.a% 'numeric')
  }

  if( rename) {
    # Geared to read-in expectations from 'read_cluster_dart3'

    rownames( locinfo) <- NULL
    rownames( seqinfo) <- NULL

    # For compatibility with other funcs that expect ClusterTempIndex column
    names( seqinfo) <- sub( 'ClusterIdx', 'ClusterTempIndex', names( seqinfo), fixed=TRUE)
    names( locinfo) <- sub( 'ClusterIdx', 'ClusterTempIndex', names( locinfo), fixed=TRUE)

    locinfo$Locus <- locinfo$CloneID
    locinfo <- locinfo[ c( 'Locus', names( locinfo) %except%
        cq( Locus, CloneID, ClusterTempIndex))]
    seqinfo$Locus <- seqinfo$CloneID
    seqinfo <- seqinfo[ c( 'Locus', names( seqinfo) %except%
        cq( Locus, CloneID, ClusterTempIndex))]
  }

  sampinfo <- make_dull( sampinfo, cq( File, MD5))
  locinfo <- make_dull( locinfo, c( 'consensus', names( locinfo) %that.match% 'Full.*Seq$'))

  x@info <- sampinfo
  x@locinfo <- locinfo
  x@seqinfo <- seqinfo

return( x)
}


"print.diploido" <-
function( x, ...) {
  x <- format( x, ...)
  NextMethod()
}


"print.loc.ar" <-
function( x, width=0, digits=getOption( 'digits.loc.ar', 0),
    use_cdot=getOption( 'cdot.loc.ar', FALSE), dot_for_0=use_cdot, ...){
  df <- attr( x, 'info')
  rownames( df) <- NULL

  # Generalized to allow >2 alleles (but fixed for all loci)
  # useful for Ref/Alt/Other count format
  # but possibly a mistake to use a class designed for biallelic genotypes...
  loci_names <- x@locinfo$Locus
  x <- unclass( x)
  nals <- dim( x)[3]
  if( width <= 0) suppressWarnings(
    width <- max( ceiling( log( x, 10)), na.rm=TRUE)
  )
  width <- pmax( 2, width) # allow NAs

  fmtstr <- paste( rep( sprintf( '%%%i.%if', width+digits, digits), nals), collapse='/')
  printio <- parse( text=sprintf( 'sprintf( fmtstr, %s)',
      paste( sprintf( 'x[,,%i]', 1:nals), collapse=', ')
      ))[[1]]

  # NAs really slow the character steps. Elim them here. Avoid apply() for speed
  sum_for_na_check <- parse( text=paste( sprintf( 'x[,,%i]', 1:nals), collapse='+'))[[1]]
  tot_is_na <- is.na( eval( sum_for_na_check)  )
  if( any( tot_is_na)) {
    for( ial in 1:nals) {
      x[ tot_is_na, ial] <- 0
    }
  }

  m <- eval( printio)
  if( any( tot_is_na)) {
    NA_str <- paste( rep( formatC( 'NA', width=width), nals), collapse='/')
    m[ tot_is_na] <- NA_str
  }

  if( use_cdot) {
    dot <- rawToChar( as.raw( 183))
    m <- gsub( ' ', dot, m, fixed=TRUE)
    if( dot_for_0) {
      m <- gsub( dot %&% '0', dot %&% dot, m, fixed=TRUE) # width >= 2 so can't get "/0"
    }
  }

  # m <- gsub( ' +', '',  sprintf( fmtstr, x[,,1], x[,,2]))
  # m <- paste( unclass( x)[,,1], unclass( x)[,,2], sep='/')
  dim( m) <- dim( x)[1:2]
  dimnames( m)[[2]] <- loci_names
  x <- cbind( df, m)
  print( x, quote=FALSE, right=TRUE, ...)
invisible( x)
}


"print.NGS_count_ar" <-
function( x,
    trailing_dot= getOption( 'trailing_dot_NGS_count_ar', FALSE),
    dot_for_0= getOption( 'dot_for_zero_NGS_count_ar', FALSE),
    center_dot=centre_dot,
    centre_dot= getOption( 'centre_dot', getOption( 'center_dot', FALSE)),
    ...
){

  attx <- attributes( x)
  df <- x@info
  rownames( df) <- NULL

  x <- unclass( x)
  n_fish <- nrow( x)
  extract.named( x@locinfo)
  start_col <- c( 0L, head( end_col, -1))+1L
  n_loci <- length( start_col)

  # If non-integer, can either use trailing dot, or merely round

  if( is.integer( x)) {
    mi <- as.character( x)
  } else {
    mi <- as.character( round( x))
    if( trailing_dot) {
      mi <- paste0( mi, '.')
    }
  }

  # Centre-dot used to look nice, but thx2 Encodings it is borked
  # NFI how to fix
  dot <- if( centre_dot) rawToChar( as.raw( 183)) else '.'
  if( dot_for_0) {
    zero <- if( trailing_dot) '0.' else '0'
    mi[ mi==zero] <- dot
  }
  dim( mi) <- dim( x)

  m <- do.on( 1:n_loci, {
      z <- mi[, start_col[ .]:end_col[ .], drop=FALSE]
      fw <- max( nchar( z))
      z[] <- format( z, width=fw, justify='right')
      if( fw>1) {
        z[] <- gsub( ' ', 'x', z, fixed=TRUE)
        z[] <- gsub( 'x0', 'xx', z, fixed=TRUE) # a zero
        z[] <- gsub( 'x', dot, z, fixed=TRUE)
      } else {
        z[ z=='0'] <- dot
      }
      sprintf( '%s  ', apply( z, 1, paste, collapse='/'))
    })
  dim( m) <- c( n_fish, n_loci) # in case only 1 fish or 1 locus, whereupon do.on will strip dim
  # dimnames( m)[[2]] <- Locus # FFS R rejects this if m is N*1 ... so
  dimnames( m) <- list( NULL, Locus) # a compulsory column, albeit "made-up by Mark" sometimes

  m[ grepl( 'NA', m, fixed=TRUE)] <- 'NA' # presumably only if binding several datasets?

  # x <- cbind( df, m) # fucken cbind turns strings into fucken factors...
  x <- cbind( df, data.frame( m, stringsAsFactors=FALSE)) # a data.frame FFS
  print( x, quote=FALSE, right=TRUE, ...)

  # Extra atts will not be printed by default (and print.data.frame ignores them)
  # Set 'x@print_atts' to character to
  dont_print_atts <- names( attx) %except% c( attx$print_atts,
      cq( class, dim, dimnames, print_atts, info, locinfo, seqinfo))
  if( length( dont_print_atts)) {
    scatn( '\nOmitting these attributes: %s', paste( dont_print_atts, collapse=', '))
  }
  for( att_to_print in attx$print_atts) {
    scatn( '\nattr( ., %s):', att_to_print)
    print( attx[[ att_to_print]], ...)
  }

  if( 'print_atts' %in% names( attx)) {


  }

  invisible( x)
}


"print.snpgeno" <-
function(x, locus_width= max( c( 4, nchar( x@diplos))),
                            chunk_size= 10000, ...){
  n_fish <- nrow( x)
  attx <- attributes( x)
  df <- attx$info
  if( is.null( df)) {
    df <- data.frame( Sample=sprintf( 'Dummy%i', seq.int( n_fish)))
  }
  rownames( df) <- NULL

  x <- unclass( x)
  n_loci <- ncol( x)
  if( !n_loci) {
    print( df)
return( invisible( x))
  }

  diplos <- c( x@diplos, '?')
  if( is.null( x@locinfo)) {
    x@locinfo <- data.frame( Locus=sprintf( 'FakeL_%i', seq.int( n_loci)) )
  }

  locus_names <- x@locinfo$Locus

  # Diplo genotypes are stored as 'raw'-mode (1 byte)
  tx <- as.vector( x)
  dim( tx) <- dim(x)
  tx <- t( tx)

  # NAs --- which shouldn't really be there, they're unusable en masse with SNPs--- go to "?"
  # OOR values too (Jan 2023)
  # tx[ tx==as.raw( 255)] <- as.raw( length( diplos))
  # diplos already extended with '?' for missign
  tx[ (tx < 1L) | (tx > length( diplos))] <- as.raw( length( diplos))

  # Trying to be efficient...
  raw_lines <- array( as.raw( 0), c( locus_width, n_loci, n_fish))
  diplos <- sprintf( sprintf( '%%%is', locus_width), diplos) # same length
  for( iw in 1 %upto% locus_width) {
    diplo_char <- charToRaw( paste( substring( diplos, iw, iw), collapse=''))

    # Can't subscript directly with raw :(
    # Converting *entire* 'tx' to 'integer' could blow up memory...
    # ... so do it in chunks
    fish_start <- 1
    while( fish_start <= n_fish) {
      fishes <- seq( from=fish_start, to=min( fish_start+chunk_size-1, n_fish))
      raw_lines[ iw,, fishes] <- diplo_char[ as.integer( tx[ ,fishes])]
      fish_start <- fish_start + chunk_size
    }
  }

  per_fish <- apply( raw_lines, 3, rawToChar, multiple=FALSE)

  # Now put locus-names *vertically*
  locus_name_length <- max( nchar( locus_names))
  vlocus_names <- sprintf( sprintf( '%%%is', locus_name_length), locus_names) # same length
  raw_lines <- array( charToRaw( ' '), c( locus_width, locus_name_length, n_loci))
  raw_lines[ locus_width,,] <- charToRaw( paste( vlocus_names, collapse=''))
  vlocus_names <- apply( raw_lines, 2, rawToChar, multiple=FALSE)

  # m <- gsub( ' +', '',  sprintf( fmtstr, x[,,1], x[,,2]))
  # m <- paste( unclass( x)[,,1], unclass( x)[,,2], sep='/')

  # Force single-line printing. No, you do not get a choice
  owidth <- options( width=9999)
  df$.... <- per_fish
  pdf <- try( capture.output( print( df, quote=FALSE, right=TRUE, ...)))
  options( owidth)
  if( pdf %is.a% 'try-error') {
stop( "Couldn't print 'info' columns")
  }

  pdf_width <- max( nchar( pdf)) # was: nchar( pdf[ 1]))
  wdf <- pdf_width - nchar( vlocus_names[1])
  substring( pdf[ 1], wdf+1, nchar( pdf[1])) <- vlocus_names[ locus_name_length]
  # strrep() is neater but only appears in R >= 3.3.1
  # vert_loc_names <- strrep( ' ', wdf) %&% vlocus_names[-locus_name_length]
  vert_loc_names <- paste( rep( ' ', wdf), collapse='') %&% vlocus_names[-locus_name_length]
  pdf <- c( vert_loc_names, pdf)

  # Lots of single lines. Might need to split. Row, or row-name, is always left-justified at the start
  longest_rowname <- min( grep( ' ', pdf, fixed=TRUE) %such.that% (.>1))-1
  longest_line <- max( nchar( pdf))
  field_name_row <- length( vert_loc_names)+1
  row_heads <- substring( pdf, 1, longest_rowname)
  pdf <- substring( pdf, longest_rowname+1)

  # ... Where to split columns, if lines are too long
  split_points <- c( 1,
      gregexpr( '(?<! ) ', pdf[ field_name_row], perl=TRUE)[[1]][-1],
      longest_line)
  breako <- longest_rowname+1
  so_far <- 0
  this_split_point <- 1
  repeat{
    next_split_point <- findInterval( so_far + owidth$width, split_points)
    if( next_split_point==this_split_point) {
      next_split_point <- this_split_point+1
    }

    so_far <- split_points[ next_split_point]
    if( is.na( so_far)) {
      so_far <- longest_line
    }

    these_lines <- row_heads %&%
        substring( pdf, split_points[ this_split_point], so_far)
    print( as.cat( these_lines))

    if( next_split_point >= length( split_points))
  break

    # ... if here, will need another group of lines
    cat( '\n')
    this_split_point <- next_split_point
  }

invisible( x)
}


"rbind.loc.ar" <-
function( ...){
## This 2023 version is completely different from an earlier, weirder
# version that used fancy merging and more elaborate checks, and which has
# been moved to 'genocalldart'
# This one is fairly vanilla...

  l <- list( ...)
  l1 <- l[[1]]

  info1 <- l1$info
  nal1 <- dim( unclass( l1))[3]
  locinfo1 <- l1$locinfo
  att1 <- attributes( l1) %without.name% cq( dim, class)

  # ... requiring similar subsetting, which You can control by setting these attributes.

stopifnot(
    all( do.on( l, . %is.a% 'loc.ar')),
    all( do.on( l, my.all.equal( .@dim[3], nal1))),
    all( do.on( l, listmeq( .$locinfo, locinfo1))),
    # Next: shooould check at least classes of locinfo cols for consistency...
    # ... but this _is_ S3, after all!
    all( do.on( l, my.all.equal( sort( names( .$info)), sort( names( info1))))),
    # Next could certainly be more thorough, eg about object shape... just check existence for now
    all( do.on( l, my.all.equal( sort( .@subset_like_both), sort( l1@subset_like_both) ) ) )
  )

  out <- guts_rbind( l, NULL) # should nicely rearrange
  for( ioth in c( 'info', l1@subset_like_both)) {
    att1[[ ioth]] <- guts_rbind( l, ioth)
  }

  if( any( duplicated( att1$info$Our_sample))){
warning( "Duplicated 'Our_sample'--- baaaaad karma")
  }

  mostattributes( out) <- c( list(
      dim= dim( out)),
      att1 # includes locinfo (unchanged) and info (changed)
    )

  oldClass( out) <- 'loc.ar'

return( out)
}


"rbind.NGS_count_ar" <-
function( ..., deparse.level=NA){
  # Hasty fudge...
  nextra <- ...length()
  x <- ..1
  xinfo <- x@info
  xli <- x@locinfo
  xsi <- x@seqinfo
  xmft <- x@mean_fish_tot
  xbod <- unclass( x)
  attributes( xbod) <- list( dim=xbod@dim)

  for( i in 2 %upto% nextra){
    y <- unclass( ...elt( i))
stopifnot( my.all.equal( dim( y@info), dim( xinfo)),
        my.all.equal( y@mean_fish_tot, xmft)
      )
    if( !my.all.equal( y@info, xinfo)){
      orig_names <- unique( c( names( xinfo), names( y@info)))
      # then drop columns that don't match. This is not the best way
      botho <- names( y@info) %that.are.in% names( xinfo)
      xinfo <- xinfo[ botho]
      yinfo <- y@info[ botho]
      # all.equal() does NOT just return TRUE or FALSE... sigh
      goodo <- mapply( function( a, b) isTRUE( all.equal( a, b)),
          xinfo, yinfo)
      xinfo <- xinfo[ goodo]
stopifnot( length( xinfo) > 0) # at least 1 field needed!
      warning( sprintf( "Dropping info column(s): %s",
          paste( orig_names %except% names( xinfo), collapse=',')))
    }

    xsi <- rbind( xsi, y@seqinfo)
    xli <- rbind( xli, y@locinfo)
    xbod <- cbind( x, y)
  }
  xbod@locinfo <- xli
  xbod@seqinfo <- xsi
  xbod@info <- xinfo
  xbod@mean_fish_tot <- xmft
  class( xbod) <- 'NGS_count_ar'
return( xbod)
}


"rbind.snpgeno" <-
function( ..., deparse.level=NA) {
  l <- list( ...)
  l1 <- l[[1]]
  ridf1 <- rowid_field( l1)

  diplos1 <- l1@diplos
  info1 <- l1$info
  locinfo1 <- l1$locinfo
  att1 <- attributes( l1) %without.name% cq( dim, class)

  # ... requiring similar subsetting, which You can control by setting these attributes.

      ### SMB to MVB -- is this what you want listmeq to do in the stopifnot?
    listshaneq <- function(a,b) { my.all.equal(names(a[order(names(a))]), names(b[order(names(b))]) ) }
stopifnot(
    all( do.on( l, . %is.a% 'snpgeno')),
    all( do.on( l, my.all.equal( .@diplos, diplos1))),
    ## all( do.on( l, listmeq( .$locinfo, locinfo1))),
    all( do.on( l, listshaneq( .$locinfo, locinfo1))),
    # Next: shooould check at least classes of locinfo cols for consistency...
    # ... but this _is_ S3, after all!
    all( do.on( l, my.all.equal( sort( names( .$info)), sort( names( info1))))),
    # Next could certainly be more thorough, eg about object shape... just check existence for now
    all( do.on( l, my.all.equal( sort( .@subset_like_both), sort( l1@subset_like_both) ) ) )
  )

  out <- guts_rbind( l, NULL)
  for( ioth in c( 'info', l1@subset_like_both)) {
    att1[[ ioth]] <- guts_rbind( l, ioth)
  }

  raws <- do.call( 'rbind', lapply( l, unclass)) # matrices
  mostattributes( raws) <- c(
      att1, # includes locinfo (unchanged) and info (changed)
      list(
      dim= dim( raws)
    ))

  oldClass( raws) <- 'snpgeno'
  if( !is.null( ridf1)){
    raws <- with_rowid_field( raws, ridf1)
  }

return( raws)
}


"rbind.with_rowid_field" <-
function( ...){
  rowids <- lapply( list( ...), rowid_field)
  rowids <- unique( unlist( rowids[[ lengths( rowids)>0]]))
  if( length( rowids) > 1){
stop( "Inconsistent rowid_fields")
  }
  x <- NextMethod( 'rbind')
  x <- with_rowid_field( x, rowids)
  if( (length( rowids)==1) &&
      (any( is.na( x[[ rowids]])) || any( duplicated( x[[ rowids]])))
    ){
warning( sprintf( 'rowid_field (%s) messed up...', rowids))
  }
return( x)
}


"read_snpgds2snpgeno" <-
structure( function(
    filename,
    locusID,
    sampleID,
    infoFrame = NULL,
    infoFields = NULL,
    locinfoFrame = NULL,
    locinfoFields = NULL,
    plateField = NULL
){
#############
stopifnot(
    requireNamespace( 'SNPRelate'),
    !missing( locusID),
    !missing( sampleID)
  )

  # Avoid having to specify every function via SNPRelate::blah
  e <- new.env( parent=asNamespace( 'SNPRelate'))
  # Can't do "multiple inheritance" from several packages, but this is sole mvbutils usage:
  e$FOR <- mvbutils::FOR 
  
  e$snpgeno <- snpgeno
  # Must also copy in parameters...
  for( x in names( formals( sys.function()))){
    # In theory, there might be a delayedAssign() version that avoids forcing
    # ... but life is too short
    e[[ x]] <- get( x)
  }

  snpg <- evalq( envir=e, {
    genofile <- snpgdsOpen(filename, allow.duplicate = TRUE)
    genos <- snpgdsGetGeno(genofile, sample.id = NULL, snp.id = NULL,
        snpfirstdim = FALSE, .snpread = NA, with.id = FALSE, verbose = TRUE)

    info <- if(!is.null(infoFields)) {
      data.frame( FOR( infoFields, read.gdsn( index.gdsn( genofile, .))),
          row.names=NULL, check.names=FALSE)
    } else if(!is.null(infoFrame)) {
      read.gdsn(index.gdsn(genofile, infoFrame))
    } else {
      data.frame( Our_sample= read.gdsn(index.gdsn(genofile, sampleID)))
    }
    names( info)[ match( sampleID, names( info), 0)] <- 'Our_sample'

    if(!is.null(plateField)) {
      names( info)[ match( plateField, names( info), 0)] <- 'Our_plate'
    }

    locinfo <- if(!is.null(locinfoFields)) {
      data.frame( FOR( locinfoFields, read.gdsn( index.gdsn( genofile, .))),
          row.names=NULL, check.names=FALSE)
    } else if(!is.null(locinfoFrame)) {
      read.gdsn(index.gdsn(genofile, locinfoFrame))
    } else {
      data.frame( Locus= read.gdsn( index.gdsn( genofile, locusID)))
    }
    names( locinfo)[ match( locusID, names( locinfo), 0)] <- 'Our_sample'

    if( ! "Locus" %in% names(locinfo)) {
        Locus <- paste( "L", 1:nrow( locinfo))
        locinfo <- cbind( Locus, locinfo)
    }

    snpgdsClose(genofile)

    cgenos <- matrix( '', nrow( genos), ncol( genos))
    cgenos[ is.na( genos)] <- 'OO'
    cgenos[ genos==0] <- 'BBO'
    cgenos[ genos==1] <- 'AB'
    cgenos[ genos==2] <- 'AAO'

    # Return just from evalq()
  return( snpgeno( cgenos, 
      diplos= get_genotype_encoding()$genotypes4_ambig, 
      info= info, locinfo= locinfo))
  })

return(snpg)
}
, doc =  mvbutils::docattr( r"{
read_snpgds2snpgeno      package:gbasics

Create a 'snpgeno' from 'snpgds' file or object


DESCRIPTION

'read_snpgds2snpgeno' is a wrapper allowing one-line conversion from a 'snpgds'-format _file_ to a 'snpgeno'. Needs the 'SNPRelate' package.


USAGE

read_snpgds2snpgeno(
    filename,
    locusID,
    sampleID,
    infoFrame = NULL,
    infoFields = NULL,
    locinfoFrame = NULL,
    locinfoFields = NULL,
    plateField = NULL)


ARGUMENTS

  filename: string giving the snpgds file name, contents formatted as per 'snpgdsExampleFileName' in package 'SNPRelate'.

  locusID: string namgiving the name for locus IDs. Cannot be blank

  sampleID: string naming the column for sample IDs. Cannot be blank

  infoFrame: optional string naming the sample metadata dataframe

  infoFields: optional character vector naming the sample metadata variables

  locinfoFrame: optional string naming the locus metadata dataframe

  locinfoFields: optional character vector naming the locus metadata variables

  plateField: optional string naming the column for sample-specific plate ID


DETAILS

Locus-specific metadata and sample-specific metadata may each be supplied as either a single dataframe (using a non-NULL 'infoFrame' and/or 'locinfoFrame'), or as a vector of field names (using 'infoFields' and/or 'locinfoFields'). If all are NULL, no metadata other than the locus ID and sample ID are read in.


EXAMPLES

# From SMB:
#if (!requireNamespace("BiocManager", quietly=TRUE))
#    install.packages("BiocManager")
# BiocManager::install("SNPRelate")

if( requireNamespace( 'SNPRelate')){
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
}
}")

)

"read_vcf2NGS_count_ar" <-
function( 
  vfilename, 
  biforce=TRUE, 
  BLOCK=1000, 
  allow_disordered_unphased= TRUE, 
  majjik_bleeble=''
){
## Reads counts from VCF file (in AD format) into 'NGS_count_ar'
## ... or reads (biallelic) genotypes 
## directly into 'snpgeno'; ref=A, any-alt=B, nulls acknowledged
## NB I did try package 'vcfR', but (1) it has huge dependencies
## and (2) it doesn't actually extract what we want anyway
## So read the file ourselves...

  # Only if called from vcf2snpgeno...
  just_genotypes <- majjik_bleeble == 'jUst_geNOtypes'
  if( just_genotypes){
    define_genotypes()
    # Create OO etc
    for( i in seq_along( genotypes4_ambig)){
      assign( genotypes4_ambig[i], as.raw( i))
    }
  }

  if( vfilename %is.not.a% 'connection'){ # then open one
    f <- file( vfilename, open='r')
    on.exit( close( f))
  } else {
    f <- vfilename
  }

  lines <- character()
  repeat{
    new_lines <- readLines( f, n=BLOCK)
    if( !length( new_lines)){
stop( "File ended before good stuff found :(")
    }

    lines <- c( lines, new_lines)
    if( any( grepl( '^\\s*#CHROM\\>', new_lines))){
  break
    }
  }

  gthead <- grep( '^\\s*#CHROM\\>', lines)
  gtlines <- lines[ -(1 %upto% gthead)]

  # There is nothing much of direct interest to *me* in the prelim stuff
  # *You* might feel differently :) in which case feel free to add stuff!

  # Mandatory fields are apparently these:
  # #CHROM  POS ID  REF ALT QUAL  FILTER  INFO
  # and then FORMAT followed by sample names
  # No options! Phew :)
  fields <- strsplit( sub( '^\\s+', '', lines[ gthead]), '\\s+')[[1]]
  mandatory_fields <- c( "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
stopifnot( my.all.equal( fields[1:9], mandatory_fields))

  n_samps <- length( fields) - 9
  samps <- fields[ -(1:9)]

  # We have No Idea how many loci there are...
  # biforce means condensing each locus into its Ref and *one* ALT, made up of all non-Refs
  # NB All VCF loci are supposed to have >1 variant
  loci <- matrix( '', 0, 8) # don't need FORMAT

  if( biforce || just_genotypes){
    # See if we can pre-allocate space
    peek_count_lines <- try( eval( parse( text='fpeek:' %&% ':peek_count_lines')))
    n_loci <- if( peek_count_lines %is.not.a% 'try-error'){
      peek_count_lines( vfilename) - gthead
      } else {
        0
      }
    n_sofar <- 0
    if( just_genotypes){
      geno <- matrix( as.raw( 0), n_loci, n_samps)
    } else { # biallelic counts
      reffo <- alto <- matrix( 0L, n_loci, n_samps)
    }
  } else {
    # NFI how many variants; must create count-holder on-the-fly
    varcounts <- matrix( 0L, 0, n_samps)
    nvar <- numeric() # per locus
  }

  AD_GT <- if( just_genotypes) 'GT' else 'AD'

  repeat{
    # Might start with 0 lines if we read just up to the FORMAT line!
    if( n_new_loci <- length( gtlines)){
      # "Lex" it
      fields <- matrix( unlist(
          strsplit( sub( '\\s+$', '', gtlines), '\\s+'), use.names=FALSE),
          nrow=n_new_loci, ncol=n_samps+9, byrow=TRUE)

      ## Locus info... for later
      # Prolly orta not-grow this if nlines known
      # but main memory use is in the counts etc, so don't bother
      loci <- rbind( loci, fields[,1:8])

      ## Counts
      FORMAT <- fields[,9] # eg GT:AD:DP:GQ:PL for one locus

      # Locate the AD (or GT) entry (per locus)
      colons_before <- nchar( gsub( '[^:]', '', sub( AD_GT %&% '.*', '', FORMAT))) # 0 means at the start

      # For efficiency, do all loci with AD in the same place at once--- since sub() requires ONE pattern
      AD <- matrix( '', nrow( fields), n_samps)
      for( iii in unique( colons_before)){
        AD[ colons_before==iii,] <- sub( sprintf( '([^:]+:){%i}([^:]+).*', iii), '\\2',
            fields[ colons_before==iii,-(1:9)])
      }
      rm( fields)

      if( just_genotypes){
        GT <- AD # clear? good!

        # Completely different bloody meaning whether PHASED vs UNPHASED
        phased <- any( grepl( '|', GT, fixed=TRUE))
        unphased <- any( grepl( '/', GT, fixed=TRUE))
        if( !xor(phased, unphased)){
stop( "Are the genotypes phased or not? No idea...")
        }

        new_geno <- matrix( BBO, n_new_loci, n_samps) # default
        if( unphased){
          if( !allow_disordered_unphased){
r"--{
 Check order. When alleles differ, lower number should always be first AFAICS from manual. However, others might disagree. if such a heresy does occur, allow user to override via arg }--"

            suppressWarnings({
              inord <- do.on( 
                  strsplit(unique( c( GT)), '/', fixed = TRUE), 
                  all( diff( as.numeric( .)) >= 0))
              inord[is.na(inord)] <- TRUE
            })

            if( !all( inord)){
  stop( sprintf( 
              'Dicey unphased "genotypes" (order of variants): %s',
                paste( unique( c(GT) )[ !inord], collapse=', ')))
            }
          } else { # 1/0 etc allowed
            # Switch any (IMO) bad orderings. 
            # Only 2 variants BTW: 0 and non-0
            GT[ grepl( '[0-9]+/0', GT)] <- '0/1'
          }

          new_geno[ GT=='./.'] <- OO
          new_geno[ startsWith( GT, '0')] <- AB
          new_geno[ GT=='0/0'] <- AAO
        } else { # phased...
          new_geno[ GT=='0|0'] <- AAO
          GT <- gsub( '[0-9][0-9]+', '9', GT) # in case 10+ alleles, and so next lines work OK
          new_geno[ grepl( '[^0][|]0', GT)] <- AB   # 1|0, 2|0, etc
          new_geno[ grepl( '^0[|][^0]', GT)] <- AB
          new_geno[ grepl( '.', GT, fixed=TRUE)] <- OO
        }

        if( n_loci==0){
          geno <- rbind( geno, new_geno)
        } else {
          geno[ n_sofar + (1 %upto% n_new_loci),] <- new_geno
        }
        n_sofar <- n_sofar + n_new_loci
      } else if( biforce){
        if( n_loci==0){ # grow reffo
          reffo <- rbind( reffo, matrix( 0L, n_new_loci, n_samps))
          alto <- rbind( alto, matrix( 0L, n_new_loci, n_samps))
        } # else it was already created at full size

        # VCF's "double null" has dots everywhere, but fewer commas than it logically should
        AD[ AD=='.'] <- '0' # this works...

        reffo[ n_sofar + (1 %upto% n_new_loci),] <-
            matrix( as.integer( sub( ',.*', '', AD)), # first number, before comma
            ncol=n_samps)

        # Add up counts for all remaining numbers in DP--- ie Alts. Number-of-alts might vary per locus
        # so... be sneaky :) and make them into (sub)expressions which R can add
        # ie "7,4,6" loses the "7," (that's the ref) and then becomes "4+6" and then...
        alto[ n_sofar + (1 %upto% n_new_loci),] <-
            matrix( as.integer( sapply(
            parse( text=gsub( ',', '+', sub( '^[^,]+,', '', AD))), # expression() with lots of elements
            eval) ), ncol=n_samps)

        n_sofar <- n_sofar + n_new_loci
      } else { # multiple variants allowed; can differ by locus
        # For speed, try to avoid looping over loci (or samples)...
        # ... via some pretty nasty vectorization
        vprev <- length( nvar)
        ncommas <- nchar( gsub( '[^,]+', '', AD))

        # All samps per locus will have the same number of commas...
        new_nvar <- rowSums( ncommas) / rowSums( ncommas > 0) + 1
        new_nvar[ is.na( new_nvar)] <- 1 # all double-nulls for *these* samples!
        nvar <- c( nvar, new_nvar)
        whichOO <- which( ncommas==0, arr.ind=TRUE)
        repOO <- do.on( new_nvar, paste( rep( '0', .), collapse=','))
        AD[ whichOO] <- repOO[ whichOO[,1]] # matrix repl; whichOO[,1] is row

        # Now every floco has the right number of comma-separated
        varcounts <- rbind( varcounts,
            matrix( as.integer( unlist( strsplit( AD, ',', fixed=TRUE), use.names=FALSE)),
            ncol=n_samps))
      } # if genotypes/biforce/mulltivar
    } # if we have lines to process

    # Get more lines...
    gtlines <- readLines( f, n=BLOCK)
    if( !length( gtlines)){
  break # we are done
    }
  } # repeat until read all

  rm( AD, gtlines)

  if( (just_genotypes || biforce) && (n_sofar < n_loci)){
    # then blank lines at the end
    n_loci <- n_sofar
    if( just_genotypes){
      geno <- geno[ 1:n_loci,]
    } else {
      alto <- alto[ 1:n_loci,]
      reffo <- reffo[ 1:n_loci,]
    }
  }

  n_loci <- nrow( loci) # just to make sure!

  # Everything will be multiallelic from now on
  if( !just_genotypes && biforce){
    varcounts <- array( 0L, c( n_loci, 2, n_samps))
    varcounts[,1,] <- reffo
    varcounts[,2,] <- alto
    varcounts <- aperm( varcounts, c( 2, 1, 3))
    dim( varcounts) <- c( n_loci * 2, n_samps)
    nvar <- rep( 2L, n_loci)
    rm( reffo, alto)
  }

  ## Locus names: boil down CHROM & POS
  # Take out all identical prefixes
  n_loci <- nrow( loci) # we might drop some later
  CHROM <- loci[,1]
  chmax <- min( nchar( CHROM))
  for( ich in 1 %upto% chmax){
    thrubb <- substring( CHROM, ich, ich)
    if( length( unique( thrubb)) > 1){
      chmax <- ich
  break
    }
  }

  rm( CHROM) # stop myself accidentally using non-loci$version
  colnames( loci) <- gsub( '\\W+', '', mandatory_fields %except% 'FORMAT')
  loci <- data.frame( loci)
  loci <- within( loci, {
    CHROM <- substring( CHROM, chmax, nchar( CHROM))
    POS <- as.integer( loci[ ,'POS']) # integer *just about* OK in theory, but...
    Locus <- sprintf( 'C%s:P%i', CHROM, POS)
  })
  loci <- loci[ c( 'Locus', names( loci) %except% 'Locus')]

  sinfo <- data.frame(
      Our_sample= samps,
      Fishtot= if( just_genotypes) rep( NA, n_samps) else colSums( varcounts),
      File= summary( f)$description[1],
      MD5= if( vfilename %is.not.a% 'connection') tools::md5sum( vfilename) else NA_integer_,
      check.rows= FALSE, row.names= NULL)
  sinfo <- make_dull( sinfo, cq( File, MD5))

  if( just_genotypes){
    geno <- t( geno) # transpose!!!
    geno@diplos <- genotypes4_ambig
    geno@locinfo <- loci
    geno@info <- sinfo
    oldClass( geno) <- 'snpgeno'
    geno@call <- deparse1( sys.call())
return( geno)
  }

  ## Scrunge it all into 'NGS_count_ar': each col is a variant within a locus...
  # ... and all seqs in a locus are grouped together
  nvar <- as.integer( nvar)
  loci$n_alleles <- nvar
  loci <- within( loci, {
    end_col <- cumsum( nvar)
    start_col <- 1L + c( 0L, head( end_col, -1))

    # Some more fields *may* be required...
    consensus <- loci$ID
    var_pos <- as.integer( sub( ':.*', '', sub( '^[^:]+:', '', loci$ID)))
  })

  ## Not genotypes; deffo want NGS_count_ar

  # Names for "full allele seqs":
  # multiallelic is ID:REF, ID:ALT.1, ID:ALT.2, etc
  # biallelic is ID:REF, ID:ALT
  # Code could be much simpler for biforce==TRUE
  # but I've gotta write the multi version anyhow :/
  loci2 <- loci[ rep( 1:nrow( loci), loci$n_alleles), cq( ID, REF, ALT, n_alleles)]
  loci2$suffix <- '.' %&% unlist( FOR( loci$n_alleles, (1:.)-1))
  loci2 <- within( loci2, {
    REFALT <- ifelse( suffix=='.0', REF, ALT)
    suffix[ n_alleles==2] <- ''
    suffix[ suffix=='.0'] <- ''
  })
  fullseq <- with( loci2, ID %&% ':' %&% REFALT %&% suffix)

  seqi <- data.frame(
      Locus=rep( loci$Locus, times=nvar),
      FullAlleleSeq=fullseq,
      count_sum= rowSums( varcounts)
    )

  NGS <- NGS_count_ar( varcounts, # must be transposed! and is!
      sampinfo= sinfo,
      locinfo= loci,
      seqinfo= seqi,
      strip_numerics= FALSE,
      rename= FALSE
    )
  NGS@call <- deparse1( sys.call())

return( NGS)
}


"read_vcf2snpgeno" <-
structure( function( 
  vfilename, 
  BLOCK=formals( read_vcf2NGS_count_ar)$BLOCK,
  allow_disordered_unphased= TRUE
){
## Just read genotypes from VCF file, not the counts
## Uses machinery of function below
## Default for BLOCK is same as function below
  res <- read_vcf2NGS_count_ar( 
      vfilename, 
      BLOCK=BLOCK,   
      allow_disordered_unphased= allow_disordered_unphased,
      majjik_bleeble='jUst_geNOtypes')
  res@call <- deparse1( sys.call())
return( res)
}
, doc =  mvbutils::docattr( r"{
read_vcf2snpgeno    package:gbasics
read_vcf2NGS_count_ar      


Read VCF genotype data


DESCRIPTION

'read_vcf2snpgeno' and 'read_vcf2NGS_count_ar' try to read the genotype part of a VCF file, and load it into a 'gbasics' object: either a 'snpgeno' (genotype calls) or a 'NGS_count_ar' (sequence counts). They discard the "metadata" at the top of the file, and the "data lines" describing parts of the genome.

The 'snpgeno' version expects a "GT" subfield containing genotypes. The 'NGS_count_ar' expects a subfield "AD" (Allele Depth).

VCF format is complicated, and this isn't fully tested; the goal is to produce something that will get through the next few steps of the kinference process. If you want something more sophisticated, please feel free to hack this code (and if you do a good job of it, please let us know!). There are other R packages out there which- to _some_ extent- handle VCF, and you may be able to use one of those to create a 'snpgds', which can then be converted by 'sngeno' (to the latter class).


USAGE

# BLOCK has same default in both functions
read_vcf2snpgeno( 
  vfilename, BLOCK=formals( read_vcf2NGS_count_ar)$BLOCK,
  allow_disordered_unphased= TRUE)
read_vcf2NGS_count_ar(
  vfilename, biforce = TRUE, BLOCK = 1000, 
    allow_disordered_unphased= TRUE, majjik_bleeble = "")


ARGUMENTS

  vfilename: Filename, or an existing connection that is already open for reading. If 'vfilename' is a connection, it is _not_ closed on exit.

  biforce: if TRUE, condense all multi-allelic loci into just one Ref (as defined in the VCF- it's the first count) and one "Alt", the latter consisting of _all_ other sequence counts for that locus.

  BLOCK: how many lines to read in at once. Default should be fine.

  allow_disordered_unphased: By default, unphased genotypes are allowed in any order, so "1/0" and "0/1" are both tolerated and they mean the same thing. I'm not sure that's actually intended according to VCF specifications, though; my initial reading was that only non-decreasing (e.g. "0/1") was allowed, but I've now relaxed it to save having to write too many explanations like this one. Set the parameter to FALSE to enforce a 'stop()' if any decreasing pairs are found.
  
  majjik_bleeble: Do not mess with this.


VALUE

'NGS_count_ar' (qv) or 'snpgeno' object. If the latter, then the genotype encoding (the 'diplos' attribute) is currently 'genotypes4_ambig', ie AB/AAO/BBO/OO. The 'info' attribute contains the following fields:

  Our_sample: column name from VCF

  Fishtot: total counts for that sample, or NA if 'snpgeno'

  File: filename, or 'description' of the connection

  MD5: 'md5sum()' of the file, or NA if a connection NB: "File" and "MD5" are character vectors, but also with class 'dull' from package 'mvbutils', so they don't clutter the printout. It's only printing that's affected; and you can see their contents via 'unclass(x$info$File)'. The 'locinfo' attribute always contains the following fields:

  Locus: a name for the locus, concocted from CHROM and POS in the VCF

  CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO: as per VCF (these are mandatory fields) and for 'NGS_count_ar' objects, 'locinfo' will also contain these, which are needed by some functions that operate on 'NGS_count_ar' objects:

  consensus: same as ID

  var_pos: info on variant position

  n_alleles: what it says There are a couple of other housekeeping fields required by 'NGS_count_ar' (qv), which need not concern us here. For 'NGS_count_ar' objects, there is also a 'seqinfo' attribute which, aside from housekeeping fields, contains the field "FullAlleleSeq", which is concocted from the ID and REF and ALT fields. I haven't tried to disentangle ALT (it doesn't look easy!), so if there's more than two allele (variant) at a locus, then the multiple ALTs are distinguished by ".1", ".2", etc after the _entire_ ALT field in the VCF. If the locus is biallelic, there's no ".1".


SEE.ALSO

'snpgeno', 'NGS_count_ar'


EXAMPLES

# package BinaryDosage has some nice small VCFs
# The one below just has genotypes
# But I'm not putting in Suggests; hassle of versions etc
# Workaround for CRANkiness:
re_bloody_quireNamespace <- get( 
    sprintf( '%s%s', 're', 'quireNamespace'), baseenv()) # anti CRANky
if( re_bloody_quireNamespace( 'BinaryDosage')){
  thrubb <- read_vcf2snpgeno( system.file( 
      'extdata/set1b.vcf', package='BinaryDosage'))
  print( thrubb)
}
# Should also have a small public example file for 'read_vcf2NGS_count_ar' but
# must have AD subfield
# I don't have a good one
}")

)

"rel.delta" <-
function(old, new) {
  infs <- !is.finite( abs(old)+abs(new))
  qq <- pmax(abs(old), abs(new))
  qq <- ifelse(qq > 0, qq, 1)
  qq <- abs(old - new)/qq
  qq[ infs] <- ifelse( old[ infs]==new[ infs], 0, 1)
  max( qq)
}


"renorm_SPA" <-
structure( function(
    K, dK, ddK,
    return_what=c( 'func', 'mulfuncby'),
    tol=formals( ridder)$tol,
    sd_half_range = 10,
    already_vectorized= TRUE,
    limits = c( -Inf, Inf),
    try_reducing_range_if_NA = TRUE
  # , ... ; should really allow extra args to K & co, and build them in...
){
  isqrt_2pi <- 1/sqrt( 2*pi)
  isqrt_ddK_0 <- 1/sqrt( ddK( 0))
  absmax <- sd_half_range * isqrt_ddK_0 # x between +/- 10SD of mean

  # mc <- as.list( match.call( expand.dots=FALSE)$...)
  if( !already_vectorized) {
    K <- Vectorize( K) # names( formals( K)) %except% names( list( ...)) kinda thing
    dK <- Vectorize( dK)
    ddK <- Vectorize( ddK)
  }

  itotto <- 1
  sfunc <- function( s, ddK_s=ddK(s)) {
      x <- dK( s)
      itotto * exp( K( s) - s * x) * sqrt( ddK_s) * isqrt_2pi
    }

  repeat{
    itotto <- 1
    itotto <- suppressWarnings( try( 1 / integ( sfunc( x),
        max( min( limits), -sd_half_range * isqrt_ddK_0),
        min( max( limits), +sd_half_range * isqrt_ddK_0)),
        silent=TRUE))
    if( itotto %is.a% 'try-error') {
      if( try_reducing_range_if_NA) {
        sd_half_range <- 0.9 * sd_half_range
        if( sd_half_range > 4) {
          warning( sprintf( "renorm_SPA: NA met, reducing integration range to %5.2f SD", sd_half_range))
        } else {
stop( "renorm_SPA: still getting NAs at +/-4SD, not good")
        }
      } else {
stop( "renorm_SPA: NA met, integration range reduction not allowed")
      }
    } else {
  break
    }
  }

  return_what <- match.arg( return_what)
  if( return_what=='mulfuncby')
return( itotto)

  # Otherwise we need the x-ready version, with (vectorized) root-finding
  xfunc <- function( x) {
    # Lower & upper bounds for s

    # Aim for exact, miss, then try to get equal dist the other side
    iddK0 <- 1/ddK( 0)
    s1 <- (x-dK(0)) * iddK0

    dK1 <- dK( s1)
    # Paranoia: might be perfect!
    bingo <- dK1==x
    if( any( bingo)) { # otherwise, if K()==Vectorize(...), it fucks up and turns it all into a list FFS
      s1[ bingo] <- 1.99 * s1[ bingo]
      dK1[ bingo] <- dK( s1[ bingo])
    }

    ddK1 <- ddK( s1)
    s2 <- s1 - 2*(dK1-x) / ddK1
    dK2 <- dK( s2)
    while( any( same_sign <- (dK1-x)*(dK2-x)>1) ) { # unlikely; try a bit further
      s2 <- s2 - same_sign * (dK1-x)/ddK1
      dK2[ same_sign] <- dK( s2[ same_sign])
    }

    # 'ridder' wants a vector func with no args
    dK_min_x <- function( s) dK(s)-x
    s <- ridder( dK_min_x, pmin( s1, s2), pmax( s1, s2), tol=tol) # root finder
    ddK_s <- ddK( s)

    # Undo the reordering
  return( sfunc( s, ddK_s) / ddK_s)
  }

return( xfunc)
}
, doc =  mvbutils::docattr( r"{
renorm_SPA      package:gbasics
renorm_SPA_cumul
inv_CDF_SPA2

Saddlepoint approximation support


DESCRIPTION

These return functions with which you can then evaluate various renormalized univariate SPAs of PDF, CDF and inverse CDF. "All" you have to do is provide the KGF and derivatives, then run one of these routines once; then you will have a single function which you can repeatedly call easily and fairly cheaply to get your SPA values, without knowing what in hell you are doing or what on earth a saddlepoint approximation is. I like to preserve a little mystique, after all.


USAGE

# PDF:
renorm_SPA(K, dK, ddK, return_what = c("func", "mulfuncby"),
    tol = formals(ridder)$tol,
    sd_half_range = 10, already_vectorized=TRUE,
    limits = c( -Inf, Inf),
    try_reducing_range_if_NA = TRUE
    )
# CDF and inverse CDF, via integration of the SPA PDF and...
# ... the monotone interpolating spline in stats::splinefun( method="hyman")
renorm_SPA_cumul(K, dK, ddK,
    sd_half_range = 10, n_pts = 2001, already_vectorized=TRUE)
# Inverse CDF directly via Lugannini-Rice-type formula
inv_CDF_SPA2( p, K, dK, ddK,
    tol = formals(ridder)$tol, already_vectorized=TRUE)


ARGUMENTS

  K, dK, ddK: KGF single-argument functions for the 1D KGF and its 1D derivatives.

  already_vectorized: Use TRUE (the default) if you have prepared 'K' etc to take vector arguments (i.e. multiple values of the intrinsically-scalar KGF parameter, to be computed "in parallel"). FALSE means that 'Vectorize' will be called to do it for you. If you can reasonably code 'K' etc yourself in a vectorized way, then that will run faster than 'already_vectorized=FALSE'. But you don't have to.

  return_what: ('renorm_SPA' only) "func" gets you a renormalized PDF SPA that you can just, like, call. "mulfuncby" gives you the renormalization constant, i.e. the scalar you should multiply the unnormalized PDF SPA by to get the renormalized version.

  sd_half_range: (not 'inv_CDF_SPA2') how many Standard Deviations each side of the mean to span in the numerical integration required for renormalization. The default of 10 is pretty massive and should be OK for "reasonable" distros, but may be too big (leading to infinities/NAs) for some distros, or too small (not capturing enough of the probability mass) for other distros. So it is up to you to either check or gamble.

  n_pts: ('renorm_SPA_cumul') number of points to base the renormalization on. The default of 2001 is often excessive; 201 has been fine in my applications (and fewer means faster).

  tol: tolerance for root-finding the SPA tilt for each desired PDF argument. (Not required in 'renorm_SPA_cumul' since no root-finding is used there.)

  limits: (currently only 'renorm_SPA') if absolute limits are available, you can pass these and it won't try integrating outside them. EG if the true distro is discrete with finite range of support, go slightly within that range. You may well get away without setting this, but there is then a risk that the default 'sd_half_range' will try to integrate in a range where the SPA breaks down altogether (e.g. beyond the support).

  try_reducing_range_if_NA: (currently only 'renorm_SPA') If the renormalizing-integration step hits an NA, then 'integrate()' will barf; it's _probably_ just because the integration range was ludicrously wide (the default is 10). So the default TRUE for this parameter tries repeatedly reducing 'sd_half_range' to 90% of previous value until either (i) success or (ii) it hits 4 SD, at which point it will 'stop()'. If you're getting NAs that "close" to the mean, it's a bad sign. FALSE means 'stop()' if the original 'sd_half_range' doesn't work.
  
  p: What quantile to calculate.


VALUE

Function(s) with a special environment within which things like the renormalization constant are embedded. 'renorm_SPA_cumul' returns a list of two functions, 'CDF' and 'inv_CDF'; the others return just one.
}")

)

"renorm_SPA_cumul" <-
function( K, dK, ddK, sd_half_range=10, n_pts=2001, already_vectorized=TRUE) {
  if( !already_vectorized) {
    K <- Vectorize( K) # names( formals( K)) %except% names( list( ...)) kinda thing
    dK <- Vectorize( dK)
    ddK <- Vectorize( ddK)
  }

  x <- 0
  SPA_s_dxds <- function( s) {
    x <<- dK( s)
    # SPA is (K(x) - s*x) / sqrt( 2*pi*ddK( s))
    # but dx/ds = ddK( s)
    exp( K(s) - s*x) * sqrt( ddK( s) / (2*pi))
  }

  sd <- sqrt( ddK( 0))
  # To norm, we'd integ X over say +/- 10 sd
  # but s ~= (x-mu) / (sd)^2 hence s-range below

  spoints <- seq( -sd_half_range, sd_half_range, length=n_pts) / sd
  pdf_s <- SPA_s_dxds( spoints) * diff( spoints[1:2]) # diff() gets the integral about right
  C <- 1 / sum( pdf_s) # should be close to 1--- and NB sum() is very accurate: robust to sorting etc

  # Renormalize. Note we have edge/middle issues--- evaluating *at* the breaks between "steps"
  pdf_s <- pdf_s * C
  cdf1 <- cumsum( pdf_s)
  cdf2 <- rev( 1-cumsum( rev( pdf_s)))  # Should be same, but appreciably different because of edge/middle

  # At each end, give more weight to the smaller-so-far sum, where rounding error has bitten less
  wt1 <- 1-cdf1
  cdf <- wt1*cdf1 + (1-wt1)*cdf2

  # NOT SURE that ties=mean is good here...
  inv_CDF_bod <- splinefun( cdf, x,  method='hyman', ties=mean)
  inv_CDF <- function( p) inv_CDF_bod( p) # sensibler arg name, and no deriv arg
  CDF <- splinefun( x, cdf, method='hyman')
returnList( CDF, inv_CDF)
}


"ridder" <-
structure( function( FUN, lo=NULL, hi=NULL, tol=1e-6, skip_bounds=FALSE) {
  # Ridder's root-finder for all loci at once, given constant error rate e
  SIGN <- function( a, b) abs( a) * sign( b)

  if( !skip_bounds) {
    # Start by bracketing root
    xl <- (lo*2+hi)/3
    xh <- (lo+2*hi)/3
  } else { # why bother to narrow the range, if already OK?
    xl <- pmin( lo, hi) # paranoia!
    xh <- pmax( lo, hi)
  }

  fl <- FUN( xl)
  fh <- FUN( xh)
  fincr <- fl < fh

  repeat{
    toobig <- xor( fincr, fl < 0)
    if( !any( toobig))
  break
    xh[ toobig] <- xl[ toobig]
    fh[ toobig] <- fl[ toobig]
    xl[ toobig] <- (lo+xl)[ toobig] / 2
    fl <- FUN( xl)
  }

  repeat{
    toobig <- xor( fincr, fh > 0)
    if( !any( toobig))
  break
    xh[ toobig] <- (hi+xh)[ toobig] / 2
    fh <- FUN( xh)
  }

  dun <- xh < xl # FALSE
  ans <- xl - (xh-xl)/10
  repeat{
    xm <- 0.5 * (xl+xh)
    fm <- FUN( xm)
    s <- sqrt( fm*fm - fl*fh)
    dun <- dun | (s==0)
    s[ dun] <- 1 # avoid sillies
    fm[ dun] <- 0

    # Avoid pointless Inf when s==0, since xnew will not be used in that case anyway...
    xnew <- xm + (xm-xl)*sign( fl-fh)*(fm/(s+(s==0)))

    # if abs(xnew-ans) <= xtol...
    dun <- dun | (abs( xnew-ans) <= tol)
    ans[ !dun] <- xnew[ !dun]
    fnew <- FUN( ans)

    # if fnew==0 ...
    dun <- dun | (fnew==0)
    case1 <- !dun & ( SIGN(fm,fnew) != fm)
    xl[ case1] <- xm[ case1]
    fl[ case1] <- fm[ case1]
    xh[ case1] <- ans[ case1]
    fh[ case1] <- fnew[ case1]

    case2 <- !dun & !case1 & ( SIGN(fl,fnew) != fl)
    xh[ case2] <- ans[ case2]
    fh[ case2] <- fnew[ case2]

    case3 <- !dun & !(case1 | case2) & ( SIGN( fh, fnew) != fh)
    xl[ case3] <- ans[ case3]
    fl[ case3] <- fnew[ case3]

    dun <- dun | (abs( xh-xl) <= tol)
    if( all( dun))
  break
  } # until Ridder converges

return( ans)
}
, doc =  mvbutils::docattr( r"{
ridder      package:gbasics

Parallel root-finder


DESCRIPTION

Ridder's root-finder for many univariate functions in parallel. Faster than looping over each function in turn, IME. 

Ridder's method is a very good 1D root-finding algorithm; see 'https://doi.org/10.1109%2FTCS.1979.1084580' for original, or the section in "Numerical Recipes" (probably chapter 9, for the 2007 edition) by Press, Teukolsky, Vetterling, Flannery.

One weakness of this implementation is that all components of 'FUN' get evaluated in each iteration, until the very last component has converged. It would be faster if I allowed 'FUN' to have a logical-index argument, saying which components to actually evaluate each time. 


USAGE

ridder(FUN, lo = NULL, hi = NULL, tol = 0.000001, skip_bounds = FALSE)


ARGUMENTS

  FUN: function with one _vector_ argument. If your underlying function has other args, you need to wrap it first so that you can pass a one-argument function to 'ridder'.

  lo, hi: bounds for search that should bracket the roots (can be vectors)

  tol: for roots (absolute tolerance- a scalar, duhhh)

  skip_bounds: if TRUE, start the searches exactly at the bounds; quicker, but a bad idea if infinite values result. If FALSE, the search will start inside the bounds and automatically do some bracketing, being sure to avoid the bounds.


VALUE

The roots (a vector).
}")

)

"rowid_field" <-
function( x)  UseMethod( 'rowid_field')


"rowid_field.default" <-
function( x) x@rowid_field


"rowid_field.snpgeno" <-
function( x) rowid_field( x@info)


"snpgeno" <-
structure( function( x, ...) UseMethod( 'snpgeno')
, doc =  mvbutils::docattr( r"{
snpgeno      package:gbasics
snpgeno.default
diplos
diplos<-
set_rowid_field

Class for SNP genotypes


DESCRIPTION

'snpgeno' is an S3 class for storing _already-called_ genotypes from multiple samples and loci, as well as associated information about those samples and loci. Pretty much all the functions in the 'kinference' package expect a 'snpgeno' input, in many cases requiring extra information about the loci, such as allele frequency estimates, which has usually been added by other functions in 'kinference'. You can extend the per-sample and per-locus information by adding fields in the usual R manner, as well as adding extra information that is per-sample-and-locus (such as number-of-reads); there is generally no need to make an "inherited class" for such ad-lib extensions. Printing looks IMO nice, and is succinct. See how elegantly the locus-names are shown (vertically), to save space!

There are a couple of ways to create a 'snpgeno'. One is to call 'read_vcf2snpgeno' (qv) on a VCF file that already contains called genotypes, or 'read_snpgds2snpgeno' (qv) for a SNPDGS file. Another is via the constructor 'snpgeno', which is an S3 generic with a default method (plus a method for the obscure class 'loc.ar'; see *Note*). The normal way to call the default, would be to give it a character matrix of genotypes, as well as sample-specific and locus-specific data. You can also use the default to give you an "empty" 'snpgeno' of the right size but without genocalls. That might be useful if and only if you are planning to call your own genotypes and want to create the 'snpgeno' manually, in which case the code of 'snpgeno.default' and 'snpgeno.loc.ar' may be informative.

A 'snpgeno' can be subsetted by samples and/or loci, using numeric, logical, or character indices; see *Subsetting*.

'snpgeno' specific methods currently exist for: 'rbind', 'cbind', dollar and dollar-assign (direct access to attributes such as sample-specific covariates and locus-specific information), subset (see *Subsetting*), (in)equality ('==' and '!='), 'as.character', 'print', 'str'.


.GENOTYPE.ENCODING.AND.STORAGE

(See also 'get_genotype_encoding' for the same information, written slightly differently.)

The word "genotype" can be a bit ambiguous, so we here use "genocall" to describe an already-called genotype (which may be ambiguous, null, etc) specifically for _one_ sample and _one_ locus. Genocalls are stored internally in a matrix of mode 'raw' (1 byte per entry), for efficiency; if you do 'unclass(x)[1:5, 1:3]' you can see the guts. As the next paragraph explains, the raw values are automatically converted to characters for printing, assignment, and equality-testing (about the only test you can apply to a genocall). Thus, things like 'x[1,1]=="AA"' or 'x[2,3] <- "BBO" or 'x[4,5]==x[17,5]' should work fine. To do anything more sophisticated, first call 'as.character(x)' to convert the raw values into a character matrix; don't try to handle the raw elements directly.

Normally all the genocall-encoding stuff (the mapping between raws and characters) will be set up for you automatically, through the constructor call or via 'read_vcf2snpgeno' etc. You should almost never need to worry about the raw values themselves- and it's dangerous to mess with them. But for the record: the encoding, which applies to _all_ loci, is via the 'diplos' attribute, a character vector which you can obtain by 'diplos(x)'. It is indexed by the raw values, starting at '01'. So, if your 'diplos' attribute was 'c("AB","AA","BB")', then a raw value of '02' would denote an "AA" homozygote (and would print as such). This allows for flexible encodings according to the nature of the data, including null alleles (via deliberately ambiguous genocalls such as"AAO"- either single-null, or reference homozygote) and multallelic loci with up to 6 alleles (though multiallelics are not used by package 'kinference' version 1.x).

The interpretation of the encoding is determined entirely by whatever code processes it; there's no intrinsic meaning attached to the character representations. Most 'kinference' functions assume some specific encoding, amongst one of those pre-defined by 'define_genotypes' (qv), and they will check. You can convert manually between encodings, if you are feeling brave; see 'get_genotype_encodings' for an example.

Missing genocalls should be encoded as raw value 'ff'. They (along with any raw values outside the range '1:length(diplos(x)))', such as 'as.raw(0)') will print as '"?"' and are converted to NA by 'as.character(x)', but there's otherwise no guarantee they're properly supported. Missing genocalls are specifically _not_ allowed in our 'kinference' applications, so I haven't made much effort.

You _can_ do some dodgy things, such as assigning non-character items directly to elements of 'x', but you almost certainly _shouldn't_. I'm not going to list them. Just don't.


.OTHER.DATA

Data about individual samples ("metadata" for geneticists; "covariates" for statisticians) is stored in the 'info' attribute, a dataframe which you can access via e.g. 'my_sngeno$info' (thanks to S3 magic, because 'snpgeno' is internally _not_ a list; it's just a matrix, with extra attributes). There must be as many rows in 'info' as there are in the main genocall matrix. You can add whatever sample-specific fields you like to 'x$info'. As SUBSETTING explains, you can later nominate one field as "sample ID".

The "locinfo" attribute, which must have as many rows as there are columns in 'x', can be accessed via 'x$locinfo'. One column _must_ be called "Locus". Again, you can add whatever locus-specific fields you like (including matrices) to 'x$locinfo'.

'snpgeno' methods do their damnedest to ensure that 'rownames' are NULL for 'info', 'locinfo', and any other data named in '<snpgeno>$subset_like_both' (see *User defined extras*). This might not always work... nevertheless, 'rownames' on dataframes are pretty disastrous IMO (and I suspect they're an early design decision much regretted by guRus) and if you _try_ to use them with 'snpgeno' objects, you _will_ hit some kind of trouble. But, see next subsection for a solution...


.SUBSETTING

You can subscript via 'x[i,]' or 'x[,j]' or 'x[i,j]'- but not as 'x[single_subscript]', nor via matrix-subscripting. As usual in R, the subscripts can be integers, logicals, missing, or character. If the 'j' subscript (for loci) is character, then it should match elements of 'x$locinfo$Locus'. To use characters as 'i' subscripts (i.e. by sample), see 'with_rowid_field'. However, note that many functions in 'kinference' require an 'Our_sample' field (which should certainly be unique) and will not pay attention to 'rowid_field' even if it nominates a different field. 

%%#
my_snpg <- with_rowid_field( mysnpg, 'UniqueSampleID')
my_snpg[ c( 'Abner12', 'Zoe9'),] # assuming...
# ... that those names appear in 'my_snpg$info$UniqueSampleID'
%%#

Subsetting gets applied automatically to the 'x$info' and 'x$locinfo' attributes, and also to any attributes whose names are given in 'x$subset_like_both'. Subset-replacement also works as usual.


..USER.DEFINED.EXTRAS

Aside from things that belong in 'x$info' or 'x$locinfo', you can add arbitrary extra attributes via the "$<-" operator, e.g. 'x$extra <- stuff'. These will not be printed by 'print(x)', but of course you can extract and manipulate them yourself if you want.

Occasionally, you might want to add something that is per-genocall, so also has dimensions of sample * locus (and possibly extra dimensions afterwards, too): per-allele counts would be an example. In that case, just add its name to the character-vector attribute 'subset_like_both' (which will be non-existent if there are none such). For example:

%%#
x$qualitee <- array( 0, c( nrow( x), ncol( x), 3))
x$subset_like_both <- c( x$subset_like_both, 'qualitee')

Thereafter, '$qualitee' should behave sensibly when 'x' is subsetted (and allegedly also with 'rbind' and 'cbind', though I'm skeptical). Note that it's up to you to make sure the something really has the its first two dimensions correct; this is S3, there is no safety net. There are no 'dimnames' built into these attributes, but you can add them yourself; if so, presumably the first two should always equal 'x$info$Our_sample' and 'x$locinfo$Locus'.

For extra attributes that are not part of 'x$subset_like_both': 'cbind(x,y,...)' and 'rbind(x,y,...)' should preserve their values from 'x' ie the first argument, but will discard all other copies (in 'y' and so on).


USAGE

snpgeno( x, ...) # generic
snpgeno( x, diplos, n_samples, n_loci, info, locinfo,  allow_nonchar, ...) # S3 method for default
diplos( x)
diplos( x) <- value


ARGUMENTS

  x: thing to convert. For the _default_ constructor, this would normally be a character matrix of genocalls, but NULL is also allowed for an empty result that you can fill in manually later; iff you set 'allow_nonchar=TRUE', you can also pass in 'raw' or 'integer' values too, which must not exceed 'length( diplos)'.

  diplos, value: encoding for the genocalls, usually one of those in 'define_genotypes' (qv). It's a character vector defining the genocalls to which the raw elements in 'unclass(x)' will correspond.

  n_samples, n_loci: for default constructor. If not specified, these will probably be deduced either from the dimensions of 'x', or from 'info' and 'locinfo'. Specifying them explicitly can be useful if you merely want an "empty shell" 'snpgeno'.

  info: dataframe with subject-specific usually-non-genetic data (name, date, size, ...). Must contain a field 'Our_sample' (unique identifier for that sample/"library"/replicate/...); some downstream functions may require other fields too. For the default constructor, a placeholder will be constructed if 'info' is not supplied.

  locinfo: dataframe with locus-specific data. Must include "Locus" (name/definition); some downstream functions may require other fields too. For the default constructor, a placeholder will be constructed if 'locinfo' is not supplied.

  allow_nonchar: whether to allow raw or integer genocalls as input to the default constructor, for possible "efficiency"- but then it's the user's responsibility to ensure they do correspond appropriately to 'diplos').

  ...: Extra user-defined attributes for the 'snpgeno' object (default 'snpgeno' constructor only). For 'str.snpgeno', the dots have no effect but must exist Becos R.


.NOTE

For constructing from an existing 'loc.ar' object, 'x$locinfo' must already include 'x$locinfo$pambig' (a 4-column matrix of ABCO allele frequency estimates) and 'x$geno_amb' (provisional genocalls). Those are generated as part of CSIRO's CKMR-genocalling pipeline (not yet public), and hopefully are documented in there...


VALUE

A 'snpgeno' object.


SEE.ALSO

'NGS_count_ar', 'loc.ar', 'define_genotypes'


EXAMPLES

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
    Locus= sprintf( 'L%i', 1:n_loci),
    lenseq= rsample( n_loci, 100:200), 
    chromo= rsample( n_loci, 1:24, replace=TRUE),
    poschro= round( 1e7 * runif( n_loci))
  )

## Sample information (ie covariates AKA "metadata")
sampodat <- data.frame(
      # "Our_sample" must be present: unique ID string
    Our_sample= sprintf( "S%i", 1:n_samples), 
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

  
}")

)

"snpgeno.default" <-
function(
    x=NULL,
    diplos,
    n_samples= if( !is.null( x) & !is.raw( x)) nrow( x) else nrow( info),
    n_loci= if( !is.null( x) & !is.raw( x)) ncol( x) else nrow( locinfo),
    info= data.frame( Our_sample='S' %&% seq_len( n_samples)),
    locinfo= data.frame( Locus='L' %&% seq_len( n_loci)),
    allow_nonchar= FALSE,
    ...){
## Constructor
stopifnot(
    nrow( info) == n_samples,
    nrow( locinfo) == n_loci,
    'Locus' %in% names( locinfo)
  )

  r <- matrix( as.raw( 0L), n_samples, n_loci)
  r@diplos <- diplos
  row.names( info) <- NULL
  r@info <- info
  row.names( locinfo) <- NULL
  r@locinfo <- locinfo
  attributes( r) <- c( attributes( r), list( ...))
  oldClass( r) <- 'snpgeno'

  if( !is.null( x) ){  ## numeric or raw
    # Normally want character input, but if the user insists...
    if( !is.character( x) && allow_nonchar){
        if( is.numeric( x)){
          stopifnot( my.all.equal( dim( x), dim( r)))
          ow <- options( warn=2) # disallow OOR conversions...
          on.exit( options( ow))
          x[ is.na( x)] <- 255L
          x <- as.raw( x)
#          stopifnot( (min( as.integer( x)) >= 1) & (max( as.integer( x)) <= length(diplos)) )
      }
stopifnot( is.raw( x),
     all( (unique( x) %except% as.raw( 0)) <= length( diplos))
     )
    } else { # character
stopifnot( all( (unique( x) %except% NA) %in% diplos))
    }

      r[] <- x # will coerce from character (safest) or directly use raw
      if( ! (all( as.integer( unique(c(r))) >= 1) & all( as.integer( unique(c(r))) <= length(diplos))) ) {
          warning( "some genotypes are not in 'diplos' - have you correctly specified your genotpes?")
      } ## think this is only possible if the user is giving integer genotypes, otherwise
      ## it will hit a stopifnot
  }
return( r)
}


"snpgeno.loc.ar" <-
function( x, ...){
## snpgeno constructor method for a loc.ar
  ga <- unclass( x@geno_amb)
  if( is.null( ga)) {
stop( "Needs 'geno_amb' to exist")
  }

  pambig <- x@locinfo$pambig
  if( is.null( pambig) ||
      (pambig %is.not.a% 'matrix') ||
      (ncol( pambig) != 4) ) {
stop( "Needs valid 'pambig' matrix in 'x$locinfo'")
  }

  ga@info <- x@info
  li <- x@locinfo
  li$pbonzer <- pambig
  li$pambig <- NULL
  li$useN <- 4
  li$snerr <- matrix( 0, nrow=nrow( li), ncol=4, dimnames=list( NULL, cq( AA2AO, AO2AA, BB2BO, BO2BB)))
  # hopefully don't need perr...

  ga@locinfo <- li
  class( ga) <- 'snpgeno'
returnList( ga)
}


"splug_transform" <-
function( delta) {
  # Define a controlled monotone mapping from [0,1] to [0,1]: 0 -> 0, 1 -> 1
  # Mapping will be determined by delta, a vec of vals in [-inf,inf]
  # delta=rep(0,n) should return linear map
  # Thing returned is a map function that can be applied to arby
  #  input vec (all elts in [0,1])

  n <- length( delta)
  stopifnot( n>0)

  p <- rep( 0, n+2)
  for( i in 1:n) {
    p[ i+1] <- p[ i] + (1-p[ i]) * inv.logit( logit( 1 / (n+2-i)) + delta[ i])
  }
  p[ n+2] <- 1

  # Sneakily, return the INVERSE, which allows more drastic tmfns
  #return( approxfun( p, seq( 0, 1, length=n+2)))
  # Smooth version is nice
  return( splinefun( p, seq( 0, 1, length=n+2), method='hyman'))
}


"str.loc.ar" <-
function( object, loci=TRUE, keys=TRUE, ...){
  cat( "'loc.ar' object with", nrow( object), "samples and", ncol( object), "loci and",
      ncol( object@info), "sample fields\n")

  if( keys) {
    scatn( '  Key fields: %s', str( names( object@info)))
  }

  if( loci) {
    scatn( '  Loci: %s', str( colnames( object)))
  }

  invisible()
}


"str.NGS_count_ar" <-
function( object, loci=TRUE, keys=TRUE, ...){
  cat( "'NGS_count_ar' object with", nrow( object), "samples and",
      nrow( object@locinfo), "loci and",
      nrow( object@seqinfo), "alleles and",
      ncol( object@info), "key fields\n")

  if( keys) {
    scatn( '  Key fields: %s', str( names( object@info)))
  }

  if( loci) {
    scatn( '  Loci: %s', str( object@locinfo$Locus))
  }

  invisible()
}


"str.snpgeno" <-
structure( function( object, loci=TRUE, keys=TRUE, ...){
  scatn( "'snpgeno' object with %i samples and %i loci and %i sample fields",
      nrow( object), ncol( object), ncol( object@info))
      
  if( keys) {
    cat( '  Key fields:\n')
    str( names( object@info))
  }

  if( loci) {
    cat( '  Loci:\n')
    str( object$locinfo$Locus)
  }

  slb <- object$subset_like_both
  if( length( slb)){
    cat( '  Extra per-genotype info\n')
    str( attributes( object)[ slb])
  }

  atts <- atts( object) %except% cq( diplos, info, locinfo, subset_like_both)
  atts <- atts %except% slb
  if( length( atts)){
    cat( '  Other extra attributes:\n')
    str( attributes( object)[ atts])
  }

  invisible()
}
, doc =  mvbutils::docattr( r"{
str.snpgeno      package:gbasics
str.loc.ar
str.NGS_count_ar
str

Summaries for various genotype classes


DESCRIPTION

Default 'str' causes horrible crashes on these objects: too long, or something. So, use these instead.


USAGE

str( object, loci = TRUE, keys = TRUE, ...) # S3 method for snpgeno
str( object, loci = TRUE, keys = TRUE, ...) # S3 method for loc.ar
str( object, loci = TRUE, keys = TRUE, ...) # S3 method for NGS_count_ar


ARGUMENTS

  object: of whatever class

  loci: TRUE or FALSE to show (some) locus names

  keys: TRUE or FALSE to show (most or all) sample-specific fields

  ...: ignored AFAIK
}")

)

"tail.loc.ar" <-
function (x, n = 6, nl=5, ...) {
stopifnot(
    length(n) == 1L,
    length( nl)==1
  )
  nrx <- nrow( x)
  n <- if (n < 0L)
      max(nrx + n, 0L)
    else
      min(n, nrx)

  ncx <- ncol( x)
  nl <- if( nl < 0)
      max( ncx+nl, 0)
    else
      min( nl, ncx)

  x[ seq.int(to= nrx, length.out= n),
      seq.int( to= ncx, length.out= nl),, drop = FALSE]
}


"tail.NGS_count_ar" <-
function (x, n = 6, nl=5, ...) {
stopifnot(
    length(n) == 1L,
    length( nl)==1
  )

  nrx <- nrow( x)
  n <- if (n < 0L)
    max(nrx + n, 0L)
   else
    min(n, nrx)

  ncx <- ncol( x)
  nl <- if( nl < 0)
    max( ncx+nl, 0)
   else
    min( nl, ncx)

  x[seq.int(to = nrx, length.out = n), seq.int( to=ncx, length.out=n)]
}


"tail.snpgeno" <-
function (x, n = 6, nl=5, ...) {
stopifnot(
    length(n) == 1L,
    length( nl)==1
  )

  nrx <- nrow( x)
  n <- if (n < 0)
      max(nrx + n, 0)
    else
      min(n, nrx)

  ncx <- ncol( x)
  nl <- if( nl < 0)
      max( ncx+nl, 0)
    else
      min( nl, ncx)

  x[ seq.int(to= nrx, length.out= n),
    seq.int( to= ncx, length.out= nl)]
}


"undiploido" <-
function( df){
#######
# data.frame with columns as diploids gets turned into pairs-of-cols
  diploids <- sapply( df, function( x) any( class(x) == 'diploido'))
  if( any( diploids)) {
    diploids <- names( df)[ diploids]
    for( idip in  diploids) {
      df[[ idip %&% '.1']] <- df[[ idip]][,1]
      df[[ idip %&% '.2']] <- df[[ idip]][,2]
    }
    df <- df %without.name% diploids
  }
return( df)
}


"unloc.ar" <-
function( la){
  df <- attr( la, 'info')
  rownames( df) <- NULL
  loci <- dimnames( la)[[2]]
  la <- unclass( la)
  apla <- aperm( la, c( 1,3,2))
  dim( apla) <- c( dim( la)[1], prod( dim(la)[2:3]))
  dimnames( apla) <- list( NULL, c( matrix( c( loci %&% '.1', loci %&% '.2'), nrow=2, byrow=TRUE)))
  df <- cbind( df, apla)
  df
}


"with.loc.ar" <-
function( data, expr, ...)
  eval( substitute( expr), attr( data, 'info'), enclos=parent.frame())


"with_rowid_field" <-
structure( function( x, rowid_field) UseMethod( 'with_rowid_field')
, doc =  mvbutils::docattr( r"{
with_rowid_field      package:gbasics
rowid_field

Facilitate sample-based subscripting of sngpenos


DESCRIPTION

'with_rowid_field' can be applied to a 'snpgeno' to specify which "sample ID field" to use when subsetting a 'snpgeno' by row (ie sample), using character "sample ID" of your choice instead of numeric or logical index. This is often a good idea. If the rowid field has been set, then because 'find_HSPs' etc will by default label its kin-pairs with that field, rather than with row-numbers; then, if you subset the main dataset, it's easy to just subset the kin-pairs too (you don't have to keep track of changing row numbers).

You can already subscript a 'snpgeno' by locus with a character index, which is looked up in '<snpgeno>$locinfo$Locus'. The choice of "Locus" for field name is pretty obvious and uncontroversial! But for samples, there's no corresponding uniquely obvious field name so you can use 'with_rowid_field' to specify the one you want.

Specifically, that field in '<snpgeno>$info' will be used to look up character-mode 'i' in '<snpgeno>[i,...]' subsetting. Then you can subset to samples using their "identifiers", rather than having to use logical or numeric lookup. Also, if 'with_rowid_field' has been called, then subsequent calls to 'find_HSPs' (qv) and friends will return the rowid_fields in 'i' and 'j' rather than row-numbers, which is much less error-prone for downstream use.

You certainly don't have to call 'with_rowid_field' on your 'snpgeno' (it didn't exist until 2023...), but unless you do, you'll be stuck with logical/character subsetting by sample. Of course, you can still use those afterwards, as well.

'rowid_field' is a convenience lookup function to remind yourself which field you specified...


.DETAILS

'with_rowid_field' is actually a generic which can mark two (currently) classes of object: 'snpgeno' and 'data.frame' (or descendents of the latter). If you call it on a 'snpgeno', it actually gets applied to '<snpgeno>$info' which is a 'data.frame'; this means that you can subset not just '<snpgeno>' itself, but also '<snpgeno>$info', with character sample info. Note that this might theoretically _contradict_ the default behaviour of a 'data.frame', which is to look up character row indices in the highly-unreliable 'rownames(<data.frame>)'; I force rownames _off_ (ie NULL) in a '<snpgeno>$info' so this shouldn't be a problem.

'with_rowid_field' augments the S3 class of a 'data.frame', but not of a 'snpgeno' (the functionality of using rowid is built into the latter, and the rowid label is attached to '<snpgeno>$info' not to '<snpgeno>' itself). I know it does say "DETAILS" here but there are limits, so I'm not going to explain why.

There are methods for '[', '[<-', 'rbind', and 'cbind', which attempt to do the right thing...


USAGE

with_rowid_field( x, rowid_field) # generic
rowid_field( x) # generic


ARGUMENTS

  x: thing that you'll want to do character-based row subsets on

  rowid_field: which column name to look up.


VALUE

'with_rowid_field' returns the original object marked (somewhere in its internals) with a 'rowid_field' attribute. 'rowid_field' returns a string, or NULL.


EXAMPLES

## Not run:

sn <- snpgeno(...) # with "sampID" as a field in 'info'
try( sn[ 'samp_X12',]) # crash
sn <- with_rowid_field( sn, 'sampID')
sn[ 'samp_X12',] # goodo
rowid_field( sn) # [1] "sampID"
rowid_field( sn$info) # ditto
sn$info[ 'samp_X12', ] # whatever
## End(Not run)
}")

)

"with_rowid_field.data.frame" <-
function( x, rowid_field){
## Method for data.frame and descendencts; there's no default method
stopifnot(
    x %is.a% 'data.frame',
    is.character( rowid_field),
    length( rowid_field)==1,
    rowid_field %in% colnames( x)
  )

  # let it _possibly_ work on descendents of data.frame
  attr( x, 'rowid_field') <- rowid_field
  oldClass( x) <- unique( c( 'with_rowid_field', oldClass( x)))
x
}


"with_rowid_field.snpgeno" <-
function( x, rowid_field){
  x@info <- with_rowid_field( x@info, rowid_field)
return( x)
}


"xfactor" <-
function( x, exclude=if( is.factor( x) && any( is.na( levels( x)))) NULL else NA) {
  if( is.factor( x)) {
    levs <- levels( x) %except% exclude
    if( is.null( exclude) && any( is.na( x))) {
      ax <- attributes( x)
      x <- unclass( x)
      x[ is.na( x)] <- length( ax$levels)+1
      ax$levels <- c( ax$levels, '\001')
      # Bloody R nannies have made this hard
      x <- factor( ax$levels[ x], levels=ax$levels)
    }
  } else
    x <- factor( x, exclude=exclude)
  x
}

