\name{Coins}
\alias{RPcoins}
\alias{RPaverage}
\alias{Coins}
\alias{Average}
\title{Random coin method}
\description{
 The random coin method (or dilution method) is a simulation method for
 stationary Gaussian random fields. It is based on the following procedure:
 For a stationary Poisson point process on \eqn{{\bf R}^d}{R^d}
 consider the random field
 \deqn{Y(y) = \sum_{x\in X} f(y-x)}{Y(y) = \sum_{x\in X} f(y-x)}
 for a function \eqn{f}{f}. The covariance of \eqn{Y}{Y} is
 proportional to the convolution
 \deqn{C(h) = \int f(x)f(x+h) dx }{C(h) = \int f(x)f(x+h) dx}
 If the intensity of the Poisson point process increases, the
 random field \eqn{Y}{Y} approaches a Gaussian random field
 with covariance function \eqn{C}{C}.
}

\usage{
RPcoins(phi, shape, boxcox, intensity, method) 

RPaverage(phi, shape, boxcox, intensity, method) 

}

\arguments{
  \item{phi}{object of class \code{\link[=RMmodel-class]{RMmodel}};
    specifies the covariance function of the Poisson process;
    either \code{phi} or \code{shape} must be given.
  }
  \item{shape}{object of class \code{\link[=RMmodel-class]{RMmodel}};
    specifies the function which is attached to the Poisson points;
    note that this is not the covariance function of the simulated
    random field.}
 \item{boxcox}{the one or two parameters of the box cox transformation.
  If not given, the globally defined parameters are used.
  See \command{\link{RFboxcox}} for details.
 }
 \item{intensity}{positive number, intensity of the underlying Poisson
    point process.
  }
  \item{method}{integer.
    \describe{
      \item{0}{dilution of simulation window (smallest rectangle that
	covers the points) by the support of the coin. In case of random
	coins the height of the coin is corrected accordingly.
	This method only works for coins that have compact support for
	each realisation.}
      \item{1}{Idea is to dilute the set of locations
	by the support of the coin. This idea is realized by deviding
	the space into tiles and to delete those tiles that cannot
	contribute. 
  	This method only works for coins that have compact support for
	each realisation.
      }
      \item{2}{
	intensity chosen according the Berry-Esseen criterion of
	convergence, see Ballani & Schlather (2020).	
      }
    }
    The default value may change according to the current knowledge and
    state of programming.

    Default: 0
  }
}
 

 

%\details{}



\value{
 \command{\link{RPcoins}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.

}

\references{
  \itemize{    
   \item Ballani, F, Schlather, M (2020). In preparation,
 \item Lantuejoul, C. (2002)
 \emph{Geostatistical Simulation: Models and Algorithms.}
 Springer.
 }}
\me
\seealso{ \link{Gaussian},
  \link{RP},
 \command{\link{RPhyperplane}},
 \command{\link{RPspectral}},
 \command{\link{RPtbm}}.
}

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

x <- seq(0, 10, 0.1)
z <- RFsimulate(RPcoins(RMcircular(), intensity=0.1), x, x, seed=0)
plot(z)

\dontshow{\dontrun{
x <- seq(0,25, 0.02)
model <- RPcoins(RMspheric())
z <- RFsimulate(model, x, x, spConform=FALSE) # takes 20 sec
Print(sd(as.vector(z)), mean(z))
image(z)


### Gaussian field approximated by Poisson fields
x <- seq(0,10, 0.02)
for (intensity in c(1, 10, 100)) {
 z <- RFsimulate(x=x, model=RPcoins(RMspheric(), intensity = intensity))
 plot(z)
}
}}

\dontshow{FinalizeExample()}}


\keyword{methods}
