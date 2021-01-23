\name{RMdeclare}
\alias{RMdeclare}
\alias{RM_DECLARE}
\title{Declaration of dummy variables for statistical inference}
\description{
  The only purpose of this function is the declaration of dummy
  variables for defining more complex relations between parameters that
  are to be estimated. The later is used in conjunction with the
  argument \code{params} of many \link{RFfunctions}.

  Its formal value as this model is identically zero, independently of
  the variables declared. So, this model can be added to an existing
  model as if it were a part of the trend.
}
\usage{
RMdeclare(...)
}
\arguments{
  \item{...}{
    the name of any
    additional parameter, which can be a scalar, vector, or matrix.
    In the \code{params} list of \link{RFfunction}, this name must appear on
    the left hand side of a definition of an object
    (frequently with NAs inside). Arguments must all be of the form
    \code{arg=arg}, e.g. \code{RMdeclare(u=u)}.
  }
}

\value{
 \command{\link{RMdeclare}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RFformula}},
 \command{\link{RFformulaAdvanced}},
 \command{\link{RMmodel}}.
}


\keyword{spatial}
\keyword{models}
\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

## The following two examples illustrate the use of RMdeclare and the
## argument 'params'. The purpose is not to give nice statistical models 

x <- seq(1, 3, 0.1)
## note that u is of no use in the simulation,
## but it will be used when estimating
model <- ~ RMexp(sc=sc1, var=var1) + RMgauss(var=var2, sc=sc2) + RMdeclare(u=u)
p <- list(sc1=2, var1=3, sc2=4, var2=5, u=0)
z <- RFsimulate(model = model, x=x, y=x, params=p)
plot(z)

## note that the model remains the same, only the values in the
## following list change. Here, sc1, var1, sc2 and u are estimated
## and var2 is given by a forula.
p.fit <- list(sc1 = NA, var1=NA, var2=~2 * u, sc2 = NA, u=NA)
lower <- list(sc1=20, u=5)
upper <- list(sc2=1.5, sc1=100, u=15)
f <- RFfit(model, data=z, params=p.fit, lower = lower, upper = upper)
print(f)


## The second example shows that rather complicated constructions are
## possible, i.e., formulae involving several variables, both known ('abc')
## and unknown onesz ('sc', 'var'). Note that there are two different
## 'var's a unknown variable and an argument for RMwhittle
\dontrun{%hier ein speicherfehler in Randomfieldsutils - falsch
	 %erkanntes lin modell??
model2 <- ~ RMexp(var=sc) + RMwhittle(var = g, nu=Nu) + RMnugget(var=nugg) +
            RMexp(var=var, Aniso=A) + RMdeclare(CCC=CCC, DD=DD)
p.fit <- list(g=~sc^1.5,
              nugg=~sc * var * abc,
              sc=NA,
              var=~DD,
              Nu=NA,
              abc=123,
	      A = ~ matrix(c(1, 2, DD * CCC, CCC), nc=2),
              CCC = NA, DD=NA)
lower <- list(sc=1, CCC=1, DD=1)
upper <- list(sc=100, CCC=100, DD=100)
f2 <- RFfit(model2, data=z, params=p.fit, lower = lower, upper = upper)
print(f2)
}

\dontshow{FinalizeExample()}
}