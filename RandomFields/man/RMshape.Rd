\name{RMshape}
\alias{RMshape}
\alias{shape}
\alias{RM_SHAPE}
\title{Shape Model}
\description{
 \command{\link{RMshape}} is a pure shape model with covariance 0.
}
\usage{
RMshape(mean) %, plane, polydeg, polycoeff, arbitraryfct, fctcoeff)
}
\arguments{ % to do orthogonale Polynome
  \item{mean}{\link{RMmodel}, or numeric in case of \command{RMshape}.
    If it is numerical, it should be a vector of length \eqn{p}{p}, where
    \eqn{p}{p} is the number of variables taken into account by the
    corresponding multivariate random field
    \eqn{(Z_1(\cdot),\ldots,Z_p(\cdot))}{(Z_1(.),\ldots,Z_p(.))}; 
    the \eqn{i}{i}-th component of \code{mean} is interpreted as constant
    mean of \eqn{Z_i(\cdot)}{Z_i(.)}.
 }
% \item{plane, polydeg, polycoeff,arbitraryfct,fctcoeff}{
%   argument currently not maintained anymore.} 
}

\details{ 
 This function refers to trend surfaces in the geostatistical
 framework. Fixed effects in the mixed models framework are also
 implemented, see \command{\link{RFformula}}.
}

\note{
 Any function that is defined inside \command{RMshape} is considered as
 a single trend surface (or covariate).
 All NAs defined there are always treated as arbitrary
 parameters.

 Using uncapsulated subtraction to build up a covariance
 function is ambiguous, see the examples below.
 Best to define the trend separately, or to use
 \command{\link{R.minus}}.

}
\value{
 \command{\link{RMshape}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.
}
\references{Chiles, J. P., Delfiner, P. (1999) \emph{Geostatistics:
 Modelling Spatial Uncertainty.} New York: John Wiley & Sons.
}
\author{\marco; \martin}
\seealso{
 \command{\link{RMmodel}},
 \command{\link{RFformula}},
 \command{\link{RFsimulate}},
 \command{\link{RMplus}}
}

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

## first simulate some data with a sine and a mean as trend
repet <- 100
\dontshow{if(RFoptions()$internal$examples_reduced){warning("reduced 'repet'"); repet <- 3}} 
x <- seq(0, pi, len=10)
trend <- 2 * sin(R.p(new="isotropic")) + 3
model1 <- RMexp(var=2, scale=1) + trend
dta <- RFsimulate(model1, x=x, n=repet)


## now, let us estimate variance, scale, and two parameters of the trend
model2 <- RMexp(var=NA, scale=NA) + NA * sin(R.p(new="isotropic")) + NA
\dontshow{if(RFoptions()$internal$examples_reduced){warning("reduced 'repet'"); model2 <- RMexp(var=NA) + NA * sin(R.p(new="isotropic")) + NA}}
print(RFfit(model2, data=dta))

## model2 can be made explicit by enclosing the trend parts by
## 'RMshape'
model3 <- RMexp(var=NA, scale=NA) + NA *
          RMshape(sin(R.p(new="isotropic"))) + RMshape(NA)
print(RFfit(model2, data=dta))


## IMPORTANT:  subtraction is not a way to combine definite models
##             with trends
trend <- -1
model <- RMexp(var=0.4)
(model0 <- model + trend) ## exponential covariance with mean -1
(model1 <- model + -1)    ## same as model0
(model2 <- model + RMshape(-1)) ## same as model0
(model3 <- model - 1) ## !! this is a purely deterministic model
                      ## that decays exponentially towards -1

plot(RFsimulate(model=model0, x=x, y=x)) # expon. covariance and mean -1
plot(RFsimulate(model=model, trend=trend, x=x, y=x)) ## dito
plot(RFsimulate(model=model1, x=x, y=x)) ## dito
plot(RFsimulate(model=model2, x=x, y=x)) ## dito
plot(RFsimulate(model=model3, x=x, y=x)) ## purely deterministic model!



\dontshow{\dontrun{
##################################################
# Example 1: # 
# Simulate from model with a plane trend surface #
##################################################

#trend: 1 + x - y, cov: exponential with variance 0.01
model <- ~ RMshape(mean=1, plane = c(1,-1)) + RMexp(var=0.04)
#equivalent model:
model <- ~ RMshape(polydeg=1,polycoeff=c(1,1,-1)) + RMexp(var=0.4)
#Simulation
x <- 0:10
simulated0 <- RFsimulate(model=model, x=x, y=x)
plot(simulated0)
}}


\dontshow{\dontrun{

## PLOT SIEHT NICHT OK AUS !!

####################################################################
#
# Example 2: Simulate and fit a multivariate geostatistical model
#
####################################################################
 
# Simulate a bivariate Gaussian random field with trend
# m_1((x,y)) = x + 2*y and m_2((x,y)) = 3*x + 4*y
# where m_1 is a hyperplane describing the trend for the first response
# variable and m_2 is the trend for the second one;
# the covariance function is the sum of a bivariate Whittle-Matern model
# and a multivariate nugget effect
x <- y <- 0:10
x <- y <- seq(0, 10, 0.1)
model <- RMshape(plane=matrix(c(1,2,3,4), ncol=2)) + 
         RMparswm(nu=c(1,1)) + RMnugget(var=0.5)
multi.simulated <- RFsimulate(model=model, x=x, y=y)
plot(multi.simulated)

}}

\dontshow{\dontrun{
# Fit the Gaussian random field with unknown trend for the second
# response variable and unknown variances
model.na <- RMshape(plane=matrix(c(1, 2, NA, NA), ncol=2)) + 
            RMparswm(nu=c(1,1), var=NA) + RMnugget(var=NA)
fit <- RFfit(model=model.na, data=multi.simulated)
}}

\dontshow{\dontrun{
##################################################
#
# Example 3:  Simulation and estimation for model with 
#             arbitrary trend functions
#
##################################################

#Simulation
# trend: 2*sin(x) + 0.5*cos(y), cov: spherical with scale 3
model <- ~ RMshape(arbitraryfct=function(x) sin(x),
 fctcoeff=2) +
 RMshape(arbitraryfct=function(y) cos(y),
 fctcoeff=0.5) +
 RMspheric(scale=3)
x <- seq(-4*pi, 4*pi, pi/10)
simulated <- RFsimulate(model=model, x=x, y=x)
plot(simulated)

################# ?? !!
#Estimation, part 1
# estimate coefficients and scale
model.est <- ~ RMshape(arbitraryfct=function(x) sin(x), fctcoeff=1) +
 RMshape(arbitraryfct=function(y) cos(y), fctcoeff=1) +
 RMspheric(scale=NA)
estimated <- RFfit(model=model.est, x=x, y=x,
 data=simulated@data, mle.methods="ml")


#Estimation
# estimate coefficients and scale
model.est <- ~ RMshape(arbitraryfct=function(x) sin(x)) +
 RMshape(arbitraryfct=function(y) cos(y)) +
 RMspheric(scale=NA)
estimated <- RFfit(model=model.est, x=x, y=x,
 data=simulated@data, mle.methods="ml")



##################################################
#
# Example 4: Simulation and estimation for model with 
#            polynomial trend 
#
##################################################

#Simulation
# trend: 2*x^2 - 3 y^2, cov: whittle-matern with nu=1,scale=0.5
model <- ~ RMshape(arbitraryfct=function(x) 2*x^2 - 3*y^2,
 fctcoeff=1) + RMwhittle(nu=1, scale=0.5)
# equivalent model:
model <- ~ RMshape(polydeg=2, polycoeff=c(0,0,2,0,0,-3))
x <- 0:20		 
simulated <- RFsimulate(model=model, x=x, y=x)
plot(simulated)

#Estimation
# estimate nu and the trend term assuming that it is a polynomial
# of degree 2
model.est <- ~ RMshape(polydeg=2) + RMwhittle(nu=NA, scale=0.5)
estimated <- RFfit(model=model.est, x=x, y=x,
 data=simulated@data, mle.methods="ml")
}}
\dontshow{FinalizeExample()}}
		 
\keyword{spatial}
\keyword{models}



