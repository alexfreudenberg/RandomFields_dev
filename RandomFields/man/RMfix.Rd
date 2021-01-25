\name{RMfixcov}
\alias{RMfixcov}
\title{Fixed Covariance Matrix}
\description{
  \command{\link{RMfixcov}} is a user-defined (spatially constant or
   locally constant) covariance according to
   the given (globally or locally given) covariance matrix.
   
  It extends to the space through a Voronoi tessellation.
}
\usage{
RMfixcov(M, x, y=NULL, z=NULL, T=NULL, grid, var, proj, raw, givenM)%, norm)
}
\arguments{
 \item{M}{a numerical matrix defining the user-defined  covariance for a
   random field; the \eqn{k \times k}
   matrix should be positive definite, symmetric.
   For a \eqn{m-}variate random field
   its dimension \eqn{k} should equal
   (number of observations or simulation locations) times \eqn{m}.}
 \item{x,y,z,T,grid}{optional.
   The usual arguments as in \command{\link{RFsimulate}} to define the
   locations where the covariates are given.
 }
 \item{var,proj}{optional arguments; same meaning for any
   \command{\link{RMmodel}}. If not passed, the above
   covariance function remains unmodified.}
 
 \item{raw}{\argRaw}
 
 \item{givenM}{numerical matrix.
   In case of kriging or conditional simulation, either
   \code{raw = FALSE} or the second matrix for the conditioning
   locations is passed here (and \code{raw=TRUE}).
   Implicitely, it is assumed that the
   corrleations defined by \command{RMfixcov} are zero between the given
   locations and the locations defined by \code{x}. 
 }
% \item{norm}{\argNorm}

 % \item{vdim}{an integer value; defining the response dimension.}
}

\note{
  Starting with version 3.0.64, the former argument \code{element}
  is replaced by the \code{general} option \code{set} in
  \command{\link{RFoptions}}. 
}
\details{
  The covariances passed are implemented for the given locations.
  Within any Voronoi cell (around a given location) the correlation is
  assumed to be one.

 In particular, it is used in \command{\link{RFfit}} to define neighbour
 or network structure in the data.

%  \NormDetails
}
\value{
 \command{\link{RMfixcov}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.
}
\references{
 \itemize{
 \item Ober, U., Ayroles, J.F., Stone, E.A., Richards, S., Zhu, D., Gibbs, R.A., Stricker, C., Gianola, D., Schlather, M., Mackay, T.F.C., Simianer, H. (2012): \emph{Using Whole Genome Sequence Data to Predict Quantitative Trait Phenotypes in Drosophila melanogaster}. PLoS Genet 8(5): e1002685.
 }
}
\me
\seealso{
 \command{\link{RMcovariate}},
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}},
 \command{\link{RMuser}}.
}


\keyword{spatial}
\keyword{models}
\examples{\dontshow{StartExample(reduce=FALSE)}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again


n <- 10
C <- matrix(runif(n^2), nc=n)
(C <- C \%*\% t(C))
RFcovmatrix(RMfixcov(C), 1:n)

## The following shows that the covariance structure is correctly implemented
(M <- RFcovmatrix(RMfixcov(C), 1:n))
stopifnot(all(C == M))

## The following shows that the covariance structure is interpolated
(M <- RFcovmatrix(RMfixcov(C, 1:n), c(2, 2.1, 2.5, 3)))
idx <- c(2, 2, 3, 3)
stopifnot(all(C[idx, idx] == M))


## Example 3 showing the use in a separable space-time model
model <- RMfixcov(C, 1:n, proj="space") * RMexp(s=200, proj="time")
(z <- RFsimulate(model, x = seq(0,12, 0.5), T=1:1000))
plot(z)

\dontshow{FinalizeExample()}}