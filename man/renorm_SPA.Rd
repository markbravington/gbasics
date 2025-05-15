\name{renorm_SPA}
\alias{renorm_SPA}
\alias{renorm_SPA_cumul}
\alias{inv_CDF_SPA2}
\title{Saddlepoint approximation support}
\description{These return functions with which you can then evaluate various renormalized univariate SPAs of PDF, CDF and inverse CDF. "All" you have to do is provide the KGF and derivatives, then run one of these routines once; then you will have a single function which you can repeatedly call easily and fairly cheaply to get your SPA values, without knowing what in hell you are doing or what on earth a saddlepoint approximation is. I like to preserve a little mystique, after all.
}
\usage{
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
}
\arguments{
\item{  K, dK, ddK}{KGF single-argument functions for the 1D KGF and its 1D derivatives.}
\item{  already_vectorized}{Use TRUE (the default) if you have prepared \code{K} etc to take vector arguments (i.e. multiple values of the intrinsically-scalar KGF parameter, to be computed "in parallel"). FALSE means that \code{Vectorize} will be called to do it for you. If you can reasonably code \code{K} etc yourself in a vectorized way, then that will run faster than \code{already_vectorized=FALSE}. But you don't have to.}
\item{  return_what}{(\code{renorm_SPA} only) "func" gets you a renormalized PDF SPA that you can just, like, call. "mulfuncby" gives you the renormalization constant, i.e. the scalar you should multiply the unnormalized PDF SPA by to get the renormalized version.}
\item{  sd_half_range}{(not \code{inv_CDF_SPA2}) how many Standard Deviations each side of the mean to span in the numerical integration required for renormalization. The default of 10 is pretty massive and should be OK for "reasonable" distros, but may be too big (leading to infinities/NAs) for some distros, or too small (not capturing enough of the probability mass) for other distros. So it is up to you to either check or gamble.}
\item{  n_pts}{(\code{renorm_SPA_cumul}) number of points to base the renormalization on. The default of 2001 is often excessive; 201 has been fine in my applications (and fewer means faster).}
\item{  tol}{tolerance for root-finding the SPA tilt for each desired PDF argument. (Not required in \code{renorm_SPA_cumul} since no root-finding is used there.)}
\item{  limits}{(currently only \code{renorm_SPA}) if absolute limits are available, you can pass these and it won't try integrating outside them. EG if the true distro is discrete with finite range of support, go slightly within that range. You may well get away without setting this, but there is then a risk that the default \code{sd_half_range} will try to integrate in a range where the SPA breaks down altogether (e.g. beyond the support).}
\item{  try_reducing_range_if_NA}{(currently only \code{renorm_SPA}) If the renormalizing-integration step hits an NA, then \code{integrate()} will barf; it's \emph{probably} just because the integration range was ludicrously wide (the default is 10). So the default TRUE for this parameter tries repeatedly reducing \code{sd_half_range} to 90\% of previous value until either (i) success or (ii) it hits 4 SD, at which point it will \code{stop()}. If you're getting NAs that "close" to the mean, it's a bad sign. FALSE means \code{stop()} if the original \code{sd_half_range} doesn't work.}
\item{  p}{What quantile to calculate.}
}
\value{Function(s) with a special environment within which things like the renormalization constant are embedded. \code{renorm_SPA_cumul} returns a list of two functions, \code{CDF} and \code{inv_CDF}; the others return just one.
}
\keyword{misc}
