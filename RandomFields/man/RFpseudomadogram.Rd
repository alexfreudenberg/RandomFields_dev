\name{RFpseudomadogram}
\alias{RFpseudomadogram}
\alias{pseudomadogram}
\alias{cross-pseudomadogram}
\title{Pseudomadogram}
\description{
 The pseudomadogram of two random fields \eqn{X}{X} and \eqn{Y}{Y} is given by 
 \deqn{\gamma_\alpha(r):=\frac{1}{2} \EE
   |(X(t_{i})--Y(t_{j}))|^\alpha}

 For \eqn{X=Y} and \eqn{\alpha=1} the pseudo-madogram in the narrow
 sense is obtained.
 }
\usage{
RFpseudomadogram(model, x, y=NULL, z=NULL, T=NULL, grid, params, distances,
           dim, ..., data, bin=NULL, phi=NULL, theta = NULL,
           deltaT = NULL, vdim=NULL, alpha=1)
}
\arguments{
 \item{model,params}{\argModel }
 \item{x}{\argX}
 \item{y,z}{\argYz}
 \item{T}{\argT}
 \item{grid}{\argGrid}
 \item{distances,dim}{\argDistances}
 \item{...}{\argDots}
 \item{data}{\argData}
 \item{bin}{\argBin}
 \item{phi}{\argPhi} 
 \item{theta}{\argTheta} 
 \item{deltaT}{\argDeltaT}
 \item{vdim}{\argVdim}
}
\details{ \command{\link{RFpseudomadogram}} computes the empirical
 pseudomadogram for given (multivariate) spatial data. 

 
 The spatial coordinates \code{x}, \code{y}, \code{z}
   should be vectors. For random fields of
 spatial dimension \eqn{d > 3} write all vectors as columns of matrix x. In
 this case do neither use y, nor z and write the columns in
 \code{gridtriple} notation.

 If the data is spatially located on a grid a fast algorithm based on
 the fast Fourier transformed (fft) will be used.
 As advanced option the calculation method can also be changed for grid
 data (see \command{\link{RFoptions}}.)
 
 It is also possible to use \command{\link{RFpseudomadogram}} to calculate
 the pseudomadogram (see \command{\link{RFoptions}}).

}

\note{
  In the multivariate case, the direction matters.
  If no angles are given to calculate a directional pseudomadgram,
  the points are ordered with respect to the last dimension.
}



\value{
 \command{\link{RFpseudomadogram}} returns objects of class
 \command{\link[=RFempVariog-class]{RFempVariog}}. 
}
\references{
 Gelfand, A. E., Diggle, P., Fuentes, M. and Guttorp,
 P. (eds.) (2010) \emph{Handbook of Spatial Statistics.}
 Boca Raton: Chapman & Hall/CRL.

 Stein, M. L. (1999) \emph{Interpolation of Spatial Data.}
 New York: Springer-Verlag 
 }


\author{Jonas Auel; Sebastian Engelke; Johannes Martini; \martin}
 
\seealso{
 \command{\link{RMstable}},
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}},
 \command{\link{RFcov}},
 \command{\link{RFmadogram}}.
\command{\link{RFvariogram}}.
}

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
model <- RMbiwm(nudiag=c(1, 2), nured=1, rhored=1, cdiag=c(1, 5), 
                s=c(1, 1, 2))
n <- 2
x <- seq(0, 20, 0.1)
z <- RFsimulate(model, x=x, y=x, n=n)
emp.vario <- RFpseudomadogram(data=z)
plot(emp.vario)

\dontshow{FinalizeExample()}}

\keyword{spatial}
\keyword{models}
