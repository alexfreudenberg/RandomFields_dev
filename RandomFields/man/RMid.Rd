\name{RMidcoord}
\alias{RMidcoord}
\alias{RMid} % obsolete !!
\title{Identical Model}
\description{
  \command{RMidcoord} is the identical function \eqn{f(x) = x}
  where \eqn{x} is a vector of coordinates and \eqn{f(x)} is a model value.
 }
 \usage{
RMidcoord()%var, scale, Aniso, proj)
RMid() ## obsolete
}

%\arguments{
% \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
%   \command{\link{RMmodel}}. If not passed, the above
%   covariance function remains unmodified.}
%}

 
\value{
 \command{\link{RMidcoord}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.

}


\me

\seealso{
  \command{\link{RMmodel}}.
  \command{RMidcoord},
  \command{RMtrafo},
  \command{RMprod}
}

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

## C(x,y) = < x, y >
RFcov(RMprod(RMidcoord()), as.matrix(1:10), as.matrix(1:10), grid=FALSE)
\dontshow{FinalizeExample()}}


\keyword{spatial}
\keyword{models}

