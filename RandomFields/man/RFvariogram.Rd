\name{RFvariogram}
\alias{RFvariogram}
\alias{RFmadogram}
\alias{RFcov}
\alias{RFcovariance}
\alias{RFpseudovariogram}
\alias{RFpseudomadogram}
\alias{pseudomadogram}
\alias{cross-pseudomadogram}
\title{Empirical (Cross-)Variogram, (Generalized) Madograms, Covariance functions}
\description{
  This bundle of functions handels theoretical model values as well as
  empirical values.
}
\usage{
RFvariogram(model, x, y=NULL, z = NULL, T=NULL, grid, params,
            distances, dim, data, vdim=NULL, ...)

RFpseudovariogram(model, x, y=NULL, z = NULL, T=NULL, grid, params,
                  distances, dim, data, vdim=NULL, ...)

RFcovariance(model, x, y = NULL, z = NULL, T=NULL, grid, params,
             distances, dim, data, vdim=NULL, ...)
RFcov(model, x, y = NULL, z = NULL, T=NULL, grid, params,
      distances, dim, data, vdim=NULL, ...)


RFmadogram(model, x, y=NULL, z=NULL, T=NULL, grid, params,
           distances, dim, data, vdim = NULL,  ..., alpha=1)

RFpseudomadogram(model, x, y=NULL, z=NULL, T=NULL, grid, params,
                 distances, dim, data, vdim=NULL, ..., alpha=1)
}

\arguments{
\item{model,params}{\argModel }
 \item{x}{\argX}
 \item{y,z}{\argYz}
 \item{T}{\argT}
 \item{grid}{\argGrid}
 \item{distances,dim}{\argDistances}
 \item{...}{See sections Details and Dots}
 \item{data}{\argData}
 \item{vdim}{\argVdim}
 \item{alpha}{positive parameter value of the generalized madogram. See
   Definitions.}
}


\section{Definitions}{
  The functions are mathematically defined as
  \describe{
    \item{\command{RFvariogram}}{
      is the uncentered semi-variogram. In the multivariate case the
      cross-variogram.
    }
    \item{\command{RFcovariance} (synonym \command{RFcov})}{
      are the (cross-) covariance functions
    }
    \item{\command{RFmadogram}}{
      is the (univariate, generalized) madogram
      \deqn{\gamma_\alpha(r):=\frac12 E  |(X(t_{i})-X(t_{j}))|^\alpha}
      so that the standard madogram is given by \eqn{\alpha=1}.
      
      There is no cross-madrogram.
    }
    \item{\command{RFpseudovariogram}}{
      is the pseudo-crossvariogram. In the univariate case it is
      identical to the variogram.
    }
    \item{\command{RFpseudomadogram}}{
      of two multivariate random fields \eqn{X}{X} and \eqn{Y}{Y} is given by 
      \deqn{\gamma_\alpha(r):=\frac{1}{2} E
	|(X(t_{i})--Y(t_{j}))|^\alpha}.
    }
  }
}


\details{
  Any of the above defined functions show the following behaviour (as
  long it is defined):

  \describe{
    \item{\code{data} and \code{model} are given}{
      Then the command is equivalent to calling \command{\link{RFfit}}:
      That is, also some coordinates must be given and
      the given model, see \command{\link{RFformula}} for defining it,
      is fitted to the data. See also
      \code{\link[=RFoptionsFit]{emp_alpha}}
    }
    \item{\code{data} are given, but not a model}{
      Then also coordinates must be given and
      the empirical version is calculated. For instance,
      the empirical
      (cross-)variogram of two random fields \eqn{X}{X}
      and \eqn{Y}{Y} is given by 
      \deqn{\gamma(r):=\frac{1}{2N(r)} \sum_{(t_{i},t_{j})|t_{i,j}=r} (X(t_{i})-X(t_{j}))(Y(t_{i})-Y(t_{j}))}{\gamma(r):=1/2N(r) \sum_{(t_{i},t_{j})|t_{i,j}=r} (X(t_{i})-X(t_{j}))(Y(t_{i})-Y(t_{j}))}
      where \eqn{t_{i,j}:=t_{i}-t_{j}}{t_{i,j}:=t_{i}-t_{j}}, and where \eqn{N(r)}{N(r)} denotes the number of pairs of data points with distancevector 
      \eqn{t_{i,j}=r}{t_{i,j}=r}.


% The empirical
% (cross-)covariance function of two random fields \eqn{X}{X} and \eqn{Y}{Y} is given by 
% \deqn{C(r):=\frac{1}{N(r)} \sum_{(t_{i},t_{j})|t_{i,j}=r}
%   (X(t_{i})Y(t_{j})) - m_{X} m_{Y}}{\gamma(r):=1/N(r)
%   \sum_{(t_{i},t_{j})|t_{i,j}=r} (X(t_{i})Y(t_{j})) - m_{X} m_{Y}}
%where \eqn{t_{i,j}:=t_{i}-t_{j}}{t_{i,j}:=t_{i}-t_{j}}, \eqn{N(r)}{N(r)}
%denotes the number of pairs of data points with distancevector
%\eqn{t_{i,j}=r}{t_{i,j}=r} and where
%\eqn{m_{X} := \frac{1}{N(r)} \sum_{(t_{i},t_{j})|t_{i,j}=r}
%  X_{t_{i}}}{m_{X} := \frac{1}{N(r)} \sum_{(t_{i},t_{j})|t_{i,j}=r}
%  X_{t_{i}}} and \eqn{m_{Y} := \frac{1}{N(r)} \sum_{(t_{i},t_{j})|t_{i,j}=r}
%  Y_{t_{i}}}{m_{Y} := 1/N(r) \sum_{(t_{i},t_{j})|t_{i,j}=r}
%  Y_{t_{i}}} denotes the mean of data points with distancevector
%\eqn{t_{i,j}=r}{t_{i,j}=r}.

       
    }
    \item{no \code{data} are given, but a model and coordinates are
      given}{
      Then the theoretical model is evaluated at the given coordinates.
    }
    \item{only the model, but no \code{data} and no coordinates are
      given}{
      Then an object is returned that allows for plotting the model.
    }
  }
}

\section{Coordinates}{ 
  The spatial coordinates \code{x}, \code{y}, \code{z}
  should be vectors. For random fields of
  spatial dimension \eqn{d > 3} write all vectors as columns of matrix x. In
  this case do neither use y, nor z and write the columns in
 \code{gridtriple} notation.
 
 If the data is spatially located on a grid a fast algorithm based on
 the fast Fourier transformed (fft) will be used (in 2021 temporarily
 disabled).
 As an advanced option, the calculation method can also be changed for grid
 data (see \command{\link{RFoptions}}.)
}

\note{
  In the multivariate case, the direction matters.
  If no angles are given to calculate a directional pseudomadgram,
  for instance,
  the points are ordered with respect to the last dimension.
}


\section{Dots}{
  The \code{\dots} have different meaning in the four situations described in
  the details:
  
  \describe{
    \item{\code{data} and \code{model} are given}{
      \link{RFoptionsFit} gives an overview over the adjusting screws for
      fitting.
      The most important ones in case of anisotropic or space-time data
      are
      \describe{
	\item{bin}{\argBin}
	\item{phi}{\argPhi} 
	\item{theta}{\argTheta} 
	\item{deltaT}{\argDeltaT}
      }
    }
    \item{\code{data} are given, but not a model}{
      Then the argument \code{\link[=RFoptionsFit]{bins}}
      might be important.
    }
    \item{no \code{data} are given, but a model and coordinates are
      given}{
      When the theoretical values are calculated, the \code{\ldots} will
      be rarely used.
    }
    \item{only the model, but no \code{data} and no coordinates are
      given}{
      In case of preparing a plot,
      \code{\ldots} only a allow
      parameters for \command{\link[graphics]{par}}, see
      \command{\link{RFpar}}
    }
  }
}


\value{
  Depending of on the specified arguments,
  the function returns
   \describe{
    \item{\code{data} and \code{model} are given}{
 	an  object of class \code{\link[=RFfit-class]{RFfit}}
        }
      \item{\code{data} are given, but not a model}{
	an  object of class \code{\link[=RFempVariog-class]{RFempVariog}}
      }
    \item{no \code{data} are given, but a model and coordinates are
      given}{
      a vector, a matrix or an array, depending on the model and the
      specified coordinates.
    }
    \item{only the model, but no \code{data} and no coordinates are
      given}{
      an  object of class \code{\link[=RFplot-class]{RFplot}}.
    }
  }
}

\references{
  Gelfand, A. E., Diggle, P., Fuentes, M. and Guttorp,
  P. (eds.) (2010) \emph{Handbook of Spatial Statistics.}
  Boca Raton: Chapman & Hall/CRL.
  
  Stein, M. L. (1999) \emph{Interpolation of Spatial Data.}
  New York: Springer-Verlag 
}


\author{Sebastian Engelke; Jonas Auel; Johannes Martini; \martin}
 
\seealso{
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.
}

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

n <- 1 ## use n <- 2 for better results

## isotropic model
model <- RMexp()
x <- seq(0, 10, 0.02)
z <- RFsimulate(model, x=x, n=n)
emp.vario <- RFvariogram(data=z)
plot(emp.vario, model=model)


## anisotropic model
model <- RMexp(Aniso=cbind(c(2,1), c(1,1)))
x <- seq(0, 10, 0.05)
z <- RFsimulate(model, x=x, y=x, n=n)
emp.vario <- RFvariogram(data=z, phi=4)
plot(emp.vario, model=model)


## space-time model
model <- RMnsst(phi=RMexp(), psi=RMfbm(alpha=1), delta=2)
x <- seq(0, 10, 0.05)
T <- c(0, 0.1, 100)
z <- RFsimulate(x=x, T=T, model=model, n=n)
emp.vario <- RFvariogram(data=z, deltaT=c(10, 1))
plot(emp.vario, model=model, nmax.T=3)


## multivariate model
model <- RMbiwm(nudiag=c(1, 2), nured=1, rhored=1, cdiag=c(1, 5), 
                s=c(1, 1, 2))
x <- seq(0, 20, 0.1)
z <- RFsimulate(model, x=x, y=x, n=n)
emp.vario <- RFvariogram(data=z)
plot(emp.vario, model=model)


## multivariate and anisotropic model
model <- RMbiwm(A=matrix(c(1,1,1,2), nc=2),
                nudiag=c(0.5,2), s=c(3, 1, 2), c=c(1, 0, 1))
x <- seq(0, 20, 0.1)
dta <- RFsimulate(model, x, x, n=n)
ev <- RFvariogram(data=dta, phi=4)
plot(ev, model=model, boundaries=FALSE)

\dontshow{\dontrun{
######################################################
# NOTE: distinguish the different uses of x and y
x <- c(1,2,1)
y <- c(4,5,6)
# coordinate space 1-dimensional, evaluated at 3 points:
RFcov(model=model, x=as.matrix(x), y=as.matrix(y))
# coordinate space is 3-dimensional, evaluated at a pair of points
RFcov(model=model, x=t(x), y=t(y))
}}


\dontshow{FinalizeExample()}}

\keyword{spatial}
\keyword{models}
