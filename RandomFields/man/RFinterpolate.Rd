\name{RFinterpolate}
\alias{RFinterpolate}
\alias{kriging}
\title{Interpolation methods}
\description{
  The function allows for different methods of interpolation.
  Currently, only various kinds of kriging are installed.
}
\usage{
RFinterpolate(model, x, y = NULL, z = NULL, T = NULL, grid=NULL,
              distances, dim, data, given=NULL, params, err.model, err.params,
              ignore.trend = FALSE, ...)
}
\arguments{
 \item{model,params}{\argModel}
 \item{x}{\argX}
 \item{y,z}{\argYz}
 \item{T}{\argT}
 \item{grid}{\argGrid}
 \item{distances,dim}{\argDistances}

 \item{data}{\argData \argDataGiven
   
  If the argument
 \code{x} is missing,
 \code{data} may contain \code{NA}s, which are then replaced through imputing.
 
 }
 \item{given}{\argGiven}

 \item{err.model,err.params}{For conditional simulation and random imputing
 only. \cr\argErrmodel} 

 \item{ignore.trend}{logical. If \code{TRUE} only the
   covariance model of the given model is considered, without the trend
   part.
 }

 \item{...}{\argDots}
 
}

\note{Important options are
  \itemize{
    \item \code{method} (overwriting the automatically detected variant
    of kriging)
    \item \code{return_variance} (returning also the kriging variance)
    \item \code{locmaxm} (maximum number of conditional values before
    neighbourhood kriging is performed)
    \item \code{fillall} imputing estimates location by default
    \item \code{varnames} and \code{coordnames} in case
    \code{data.frame}s are used to tell which column contains the data
    and the coordinates, respectively.
}}

\details{
  If the argument \code{x} is missing,
  \code{data} may contain \code{NA}s, which are then replaced by
  the kriged values (imputing); 

  Different data sets behave differently from repeated data without missings
  if and only if the model is not fully specified,
  i.e., if the model contains NAs in the linear part or an NA
  in the global variance. Then, different linear parts are fitted to the
  different data sets (but the same to repeated data).
%  NOT TRUE ANYMORE:
%  The \command{\link[RFoptions]{option}} \code{apart_linear_fit}
 % forces repeated data to be treated as different data sets.
  Handling of different data sets instead of repeated data costs
  always (considerably) more time.

%  NOT TRUE ANYMORE:
%  In case of missing values and \command{\link[RFoptions]{option}}
 % \code{na_rm_lines = FALSE} (which is the default value) or imputing,
%  repeated measurements are always treated as if they were from
%  different data sets. 

 In case of intrinsic cokriging (intrinsic kriging for multivariate
 random fields) the pseudo-cross-variogram is used (cf. Ver Hoef and
 Cressie, 1991).
}
\value{
 The value depends on the additional argument \code{variance.return},
 see \command{\link{RFoptions}}. 
 
 If \code{variance.return=FALSE} (default), \code{Kriging} returns a
 vector or matrix of kriged values corresponding to the
 specification of \code{x}, \code{y}, \code{z}, and
 \code{grid}, and \code{data}.\cr
 
 \code{data}: a vector or matrix with \emph{one} column\cr
 * \code{grid=FALSE}. A vector of simulated values is
 returned (independent of the dimension of the random field)\cr
 * \code{grid=TRUE}. An array of the dimension of the
 random field is returned (according to the specification
 of \code{x}, \code{y}, and \code{z}).\cr
 
 \code{data}: a matrix with \emph{at least two} columns\cr
 * \code{grid=FALSE}. A matrix with the \code{ncol(data)} columns
 is returned.\cr
 * \code{grid=TRUE}. An array of dimension
 \eqn{d+1}{d+1}, where \eqn{d}{d} is the dimension of
 the random field, is returned (according to the specification
 of \code{x}, \code{y}, and \code{z}). The last
 dimension contains the realisations.

 If \code{variance.return=TRUE}, a list of two elements, \code{estim} and
 \code{var}, i.e. the kriged field and the kriging variances,
 is returned. The format of \code{estim} is the same as described
 above. The format of \code{var} is accordingly.
 }
\references{
\litChiles

\litCressie

\litGoovaerts

\litHC

\litWackernagel
}

\author{
  \martin; \marco
  \subsection{Author(s) of the code:}{ \martin; Alexander Malinowski; \marco}
}

\seealso{
 \command{\link{RMmodel}},
 \command{\link{RFvariogram}},
 \code{\link[=RandomFields-package]{RandomFields}},
}

\examples{\dontshow{StartExample()} %  library(RandomFields)
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

## Preparation of graphics
dev.new(height=7, width=16) 

## creating random variables first
## here, a grid is chosen, but does not matter
p <- 3:8
points <- as.matrix(expand.grid(p,p))
model <- RMexp() + RMshape(mean=1)
dta <- RFsimulate(model, x=points)
plot(dta)
x <- seq(0, 9, 0.25)
\dontshow{if (!interactive()) x <- seq(0, 5, 1.2)}

## Simple kriging with the exponential covariance model
model <- RMexp()
z <- RFinterpolate(model, x=x, y=x, data=dta)
plot(z, dta)

## Simple kriging with mean=4 and scaled covariance
model <- RMexp(scale=2) + RMshape(mean=4)
z <- RFinterpolate(model, x=x, y=x, data=dta)
plot(z, dta)

## Ordinary kriging
model <- RMexp() + RMshape(mean=NA)
z <- RFinterpolate(model, x=x, y=x, data=dta)
plot(z, dta)



## Co-Kriging
n <- 100
x <- runif(n=n, min=1, max=50)
y <- runif(n=n, min=1, max=50)
\dontshow{if (!interactive()) n <- 2}


rho <- matrix(nc=2, c(1, -0.8, -0.8, 1))
model <- RMparswmX(nudiag=c(0.5, 0.5), rho=rho)

## generation of artifical data
data <- RFsimulate(model = model, x=x, y=y, grid=FALSE)
## introducing some NAs ...
print(data)
len <- length(data)
data@data$variable1[1:(len / 10)] <- NA
data@data$variable2[len - (0:len / 100)] <- NA
print(data)
plot(data)

## co-kriging
x <- y <- seq(0, 50, 1)
\dontshow{if (!interactive()) x <- y <- seq(0, 5, 1)}
k <- RFinterpolate(model, x=x, y=y, data= data)
plot(k, data)

## conditional simulation
z <- RFsimulate(model, x=x, y=y, data= data) ## takes some time
plot(z, data)



\dontshow{\dontrun{

## alternatively 

## Intrinsic kriging
model <- RMfbm(a=1)
z <- RFinterpolate(krige.meth="U", model, x, x, data=dta)
screen(scr <- scr+1); plot(z, dta)


## Interpolation nicht korrekt
## Intrinsic kriging with Polynomial Trend
model <- RMfbm(a=1) + RMshape(polydeg=2)
z <- RFinterpolate(model, x, x, data=dta)
screen(scr <- scr+1); plot(z, dta)
}}

\dontshow{\dontrun{
## Universal kriging with trend as formula
model <- RMexp() + RMshape(arbit=function(X1,X2) sin(X1+X2)) +
 RMshape(mean=1)
z <- RFinterpolate(model, x=x, y=x, data=dta)
screen(scr <- scr+1); plot(z, dta)

## Universal kriging with several arbitrary functions
model <- RMexp() + RMshape(arbit=function(x,y) x^2 + y^2) +
 RMshape(arbit=function(x,y) (x^2 + y^2)^0.5) + RMshape(mean=1)
z <- RFinterpolate(model, x=x, y=x, data=dta)
screen(scr <- scr+1); plot(z, dta)
}}

% folgender Befehl muss unbedingt drin bleiben
close.screen(all = TRUE)

\dontshow{while (length(dev.list()) >= 2) dev.off()}
\dontshow{FinalizeExample()}}


\keyword{spatial}%-- one or more ...
