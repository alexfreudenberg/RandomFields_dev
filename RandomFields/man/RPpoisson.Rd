\name{RPpoisson}
\alias{RPpoisson}

\title{Simulation of Poisson Random Fields}
\description{
Shot noise model, which is also called moving
 average model, trigger process, dilution random field, and by several
 other names.
}

\usage{
RPpoisson(shape, phi, intensity, support)
}

\arguments{
 \item{shape}{the model, \command{\link{RMmodel}}, gives the shape
   function to be used}
 \item{phi}{alternative to define a Poisson field through the covariance
 function. Note that the field is not uniquely given by \code{phi}.}
 \item{intensity}{the intensity of the underlying stationary Poisson
   point process}
 \item{support}{positive real value. Particularly in random token models
 (shot noise models), the maximum radius of support must be set, if the
 shape function does not have compact support.
 }
 }
 
\me

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

# example 1: Posson field based on disks with radius 1
x <- seq(0,25, 0.02)
model <- RMball()
z <- RFsimulate(RPpoisson(model, intensity = 2), x)
plot(z)
par(mfcol=c(2,1))
plot(z@data[,1:min(length(z@data), 1000)], type="l")
hist(z@data[,1], breaks=0.5 + (-1 : max(z@data)))


# another way of defining the above model is via the covariance function
model <- RMtent(s=2)
z1 <- RFsimulate(RPpoisson(phi=model, intensity = 2), x)
stopifnot(all.equal(z, z1))


# note that the argument 'support' must be given if
# the scale is random, e.g. in the following example
model <- RMspheric()
z <- RFsimulate(RPpoisson(phi=model, intensity = 2, support=1), x)

# example 2: Poisson field based on the normal density function
# note that
# (i) the normal density as unbounded support that has to be truncated
# (ii) the intensity is high so that the CLT holds
x <- seq(0, 10, 0.01)
model <- RMtruncsupport(radius=5, RMgauss())
z <- RFsimulate(RPpoisson(model), x, intensity = 100)
plot(z)

# example 2: Poisson field based on the double exponential density (RMexp)
# function with random scale of standard exponential distribution (exp)
x <- seq(0, 10, 0.01)
model <- RMtruncsupport(radius=1.5, 0.2 + 2 * RMexp(s=exp()))
z <- RFsimulate(RPpoisson(model), x, intensity = 1)
plot(z)

\dontshow{FinalizeExample()}}

\seealso{
  \command{\link{RMmodel}},
  \link{RP},
 \command{\link{RPcoins}}.
}




\keyword{spatial}


