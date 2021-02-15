\name{Trend Modelling}
\alias{RMmodelsTrend}
\alias{trend modelling}
\title{Trend Modelling}
\description{
  The coding of trends, in particular multivariate trends,
  will be described here.
}

\details{
  See \link{RFcalc}, \link{RMshape} and also the examples below
  for some insight on the possibilities of trend modelling.
}

%\references{Ribeiro}
\seealso{
  \link{RFcalc}, \link{RM}, \link{RMmodels}, \link{RMshape}, \link{RMmodelsMultivariate}.
}

\examples{\dontshow{StartExample();if (!interactive()) RFoptions(modus="careless")} %ok
%RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
%##                   RFoptions(seed=NA) to make them all random again

data(ca20) ## data set originally from geoR
%ca20.df <- as.data.frame(ca20)
head(ca20.df)
RFoptions(coordnames=c("east", "north"), varnames="data") 

## covariance model with variance, scale and nugget to be estimated;
## just to abbreviate later on
M <- RMexp(var=NA, scale=NA) + RMnugget(var=NA)
\dontshow{if (!interactive()) M <- RMexp(var=NA, scale=NA)} % ok

## short definition of a trend using the fact that ca20.df is a
## data.frame
ca20.RFmod02 <- ~ 1 + altitude + M
(ca20.fit02.RF <- RFfit(ca20.RFmod02, data=ca20.df))

## long definition which also allows for more general constructions
ca20.RFmod02 <- NA + NA*RMcovariate(ca20.df$altitude) + M 
(ca20.fit02.RF <- RFfit(ca20.RFmod02, data=ca20.df))


## Note that the following also works.
## Here, the covariance model must be the first summand
ca20.RFmod02 <- M + NA + ca20.df$altitude 
print(ca20.fit02.RF <- RFfit(ca20.RFmod02, data=ca20.df))

### The following does NOT work, as R assumes (NA + ca20.df$altitude) + M
### In particular, the model definition gives a warning, and an error: 
try(ca20.RFmod02 <- NA + ca20.df$altitude + M) ### error ...%## OK
%try(ca20.fit02.RF <- RFfit(ca20.RFmod02, data=ca20.df)) ### error ...%## OK

## factors:
ca20.RFmod03 <- ~ 1 + area + M ### 
(ca20.fit03.RF <- RFfit(ca20.RFmod03, data=ca20.df))

\dontshow{FinalizeExample(); RFoptions(coordnames=NULL, varnames=NULL)}
}

\me
\keyword{spatial}
