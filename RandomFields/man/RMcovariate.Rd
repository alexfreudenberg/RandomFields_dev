\name{RMcovariate}
\alias{RMcovariate}
\alias{RM_COVARIATE}
\title{Model for covariates}
\description{
 The model makes covariates available. In most cases, its use is
 implicite by the user and is made explicit by the package.
}
\usage{
RMcovariate(formula=NULL, data, x, y=NULL, z=NULL, T=NULL, grid, raw, %norm,
            addNA, factor, extra_data)
}
\arguments{
  \item{formula, data}{formula and by which the data should be modelled,
    similar to \link[stats]{lm}.

    If \code{formula} is not given, the the linear model is given by the
    data themselves.
  }
  \item{x,y,z,T,grid}{optional.
    The usual arguments as in \command{\link{RFsimulate}} to define the
    locations where the covariates are given.
  }
 % \item{var}{optional arguments; same meaning for any
%    \command{\link{RMmodel}}. If not passed, the above
%    covariance function remains unmodified.}
%  \item{factor}{vector or matrix of numerical values.
%    The length of the vector must match the given number of locations.
%  }
% 
%  \item{var}{variance, i.e. factor multiplied to the data, which can be
%    estimated through ML
%  }
  \item{raw}{\argRaw }
%  \item{norm}{\argNorm}
  \item{addNA}{
    If \code{addNA} is \code{TRUE}, then an additional
    (linear) factor is estimated in an estimation framework.
    This parameter must be set in particular when \command{RMcovariate}
    passes several covariates.
  }
  \item{factor}{real value. From user's point of view very much the same
    as setting the argument \code{var}.}
  \item{extra_data}{logical. If covariates are used and kriging is
    performed a set of covariates is needed for the given data and for
    the locations where prediction is performed. In this case a matrix
    must be passed with two columns where the second column refers to
    the given data and must have \code{".."} as column names. (The
    shorter of the two covariate vectors has to be filled with \code{NA}s.)
    In the multivariate case, columns for the given data and the for the
    other locations alternate.  
  }
}
\note{
  \itemize{
    \item   \code{c}, \code{x} 
    also accept lists of data. However, its use is not in an advanced
    stage yet.
   }
}


\details{
  The function interpolates (nearest neighbour) between the values if
  \code{raw=FALSE}.

  If \code{x} is not given, it is assumed that the globally given
  coordinates

%  \NormDetails
}
\value{
  Formally, \command{\link{RMcovariate}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.
}
\me
\seealso{
 \command{\link{RMfixcov}},
 \command{\link{RMmodel}},
 \command{\link{RMshape}}
}
\examples{\dontshow{StartExample(reduced = FALSE)}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

z <- 0.2 + (1:10)
RFfctn(RMcovariate(z), 1:10)
RFfctn(RMcovariate(data=z, x=1:10), c(2, 2.1, 2.5, 3))

\dontshow{FinalizeExample()}}

\keyword{spatial}
\keyword{models}
