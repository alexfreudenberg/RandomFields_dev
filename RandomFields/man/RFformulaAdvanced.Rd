\name{RFformulaAdvanced}
\alias{RFformulaAdvanced}

\title{Advanced RFformula}

\description{Here, examples for much more advanced formula are given

}

\details{
  Formulae are used for \code{RFfit}ting. The following notes refer
  particularly to that situation.
  
  The following terms are evaluated \bold{once}, namely \bold{before}
  any parameter is \code{estimated}:
  \itemize{
    \item definitions, which are
    expressions in the list of \code{params}
    that are not formulae.  Usually, definitions
    contain \code{NA}s on the right hand side.
    \item any formulae and expressions within the \code{RMmodel}
  }

  Then, it will be determined which parameter values have to estimated
  (i.e., those who have value \code{NA}). 
  

  \bold{By way of contrast}, the following terms are evaluated
  \bold{during the optimization process}:
  \itemize{
    \item formulae given within the list of \code{params}
    \item arguments in \code{RMmodel} that are identical to
    a variable of a formula in \code{params}
  }

  \bold{Example, broken in parts for comments}
  \preformatted{
library(RandomFields)
## The following example is tailored to show the behaviour of 'params'.
## It is a silly example from a statistical modelling point of view.

## first we define the model:
model <- ~ # '~' is necessary if 'params' is used for model description
  RMexp(scale=4, # '4' will be immediately evaluated
        var=W) + # 'W', defined in RFfit, will be immediately evaluated
  RMgauss(scale=X, # evaluated immediately
          ## var = ~Y # would give an error, '~' is not allowed in 
          ##          conjuction with the variable Y defined in 'params'
	  ##          below
          var=Y) + 
  RMspheric(scale = Z, # see in params
            ## var = T # leads to an error, see definition of T below
            var = A # given by a complicated formula, see params
            ) +
  RMdeclare(# a dummy RMmodel just added to the genuine model used
      B     # B ist used in params and must be estimated. but is not
      ##      anywhere else in the model. So it must be declared here
      ## B = B # is the only allowed, alterative form
      ## B = ~B # not allowed
      ## B = 3+B# not allowed
      ## Bx = B # not allowed
      )

## now, we define the corresponding 'params':
V <- 5          # any value or even NA; is used in this example in 
    ##           'params' and RFfit
params <- list( # params is always a list
    X = NA,     # evaluated immediatedly, i.e. a value to be estimated
    Y = 1 + NA, # evaluated immediately! 1+NA equals NA !
    ## Z = 2 + Y# would give an error, since 'Y' is not known at the 
    ##            time 'params' is defined!
    Z = ~ W + Y,# a formula evaluated during optimisation process
    T = ~ W + NA,#only evaluated during the optimization process,
    ##            i.e. T equals NA. During the optimization process,
    ##            an NA is not considered as indicators for values 
    ##            to be estimated. It is just the value NA.
    ## U = ~ V + Z, # error as V is unknown, at time when the formula 
    ##            is evaluated
    A = ~ (1 + sin(B)) * B, # complicated formula based on a value B
    B = NA      # since this variable must be estimated it must be
    ##            declared here. Sicne it is not used anywhere in the
    ##            model, RMdeclare must be added, see above
    )

## now, lower and upper bounds are defined:
lower <- list( ## here, lower bounds known by the system can be 
    ##         overwritten; unknown lower bounds must be given.
    ## Y = 0 # not necessary, as 'var=Y' implies that X must be
    ##         non-negative; 'X' also need not be given, as 'scale=X'
    ##         implies that X must be positive    
    B = 0.1 # the system cannot have an idea of the range 'B'. So the 
    ##        lower bound must be given
    ##
    ## Note that bounds for formulae can never be given.
)
upper <- list(B=10,
              X=1.234 ## no harm to overwrite the system's upper bound
              )

## parameter estimation:
data <- data.frame(x = 1:100, data=runif(100))
(RFfit(model=model, data=data,
      params=params,
      lower=lower, upper=upper, # 'B' in RMdeclare must be given bound
      W=V, V=V)) # global variable 'V' will be unknown when the 
##                formulae in 'params' will be evaluated. Transform V
##                to a local one (or several ones as in this example)


(RFgetModelInfo(RFfit)) # NA indicates parameter values to be optimized
##                'NaN' indicates parameter values to be calculated
##                from values that are optimized
  }
}


\seealso{\link{RFformula},
  \command{\link{RMdeclare}}}


\examples{\dontshow{StartExample()} % library(RandomFields)
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

#################################################################
### the following definitions are needed in all the examples  ###
#################################################################
V <- 10
S <- 0.3
M <- 3
x <- y <- seq(1, 3, 0.1)


#################################################################
###       Example 1: simple example                          ###
#################################################################

## the following to definitions of a model and call of RFsimulate
## give the same result:
model <- RMexp(var=V, scale=S) + M
RFplot(model)
z1 <- RFsimulate(model = model, x=x, y=y)
RFplot(z1)


model <- ~ M + RMexp(var=var, scale=sc)
p <- list(var=V, sc=S, M=M)
z2 <- RFsimulate(model = model,x=x, y=y, param=p)
RFplot(z2)



#################################################################
###       Example 2: formulae within the parameter list       ###
#################################################################

## free parameters (above 'var' and 'sc') can be used
## even within the definition of the list of 'param'eters
model <- ~ RMexp(var=var, sc=sc) + RMnugget(var=nugg)
p <- list(var=V, nugg= ~ var * abs(cos(sc)), sc=S) ## ordering does not matter!
RFplot(model, params=p)
z1 <- RFsimulate(model, x, y, params=p)
RFplot(z1)
RFgetModel(RFsimulate) ## note that V * abs(cos(S) equals  9.553365

## so the above is equivalent to
model <- ~ RMexp(var=var, sc=sc) + RMnugget(var=var * abs(cos(sc)))
z2 <- RFsimulate(model, x, y, params=list(var=V, sc=S))
RFplot(z2)



#################################################################
###     Example 3: formulae for fitting (i.e. including NAs)  ###
#################################################################
## first generate some data
model <- ~ RMexp(var=var, sc=sc) + RMnugget(var=nugg)
p <- list(var=V, nugg= ~ var * abs(cos(sc)), sc=S) 
RFplot(model, params=p)
z <- RFsimulate(model, x, y, params=p, n=10)

## estimate the parameters
p.fit <- list(sc = NA, var=NA, nugg=NA)
print(f <- RFfit(model, data=z, params=p.fit))

## estimation with a given boundaries for the scale
p.fit <- list(sc = NA, var=NA, nugg=NA)
lower <- list(sc=0.01)
upper <- list(sc=0.02)
print(f <- RFfit(model, data=z, params=p.fit, lower = lower, upper = upper))




#################################################################
###   Example 4 (cont'd Ex 3): formulae with dummy variables  ###
#################################################################
V <- 10
S <- 0.3
M <- 3
x <- y <- seq(1, 3, 0.1)

model <- ~  RMexp(sc=sc1, var=var) + RMgauss(var=var2, sc=sc2) +
            RMdeclare(u=u) ## introduces dummy variable 'u'
p.fit <- list(sc1 = NA, var=NA, var2=~2 * u, sc2 = NA, u=NA)
lower <- list(sc1=20, u=5)
upper <- list(sc2=1.5, sc1=100, u=15)
print(f <- RFfit(model, data=z, params=p.fit, lower = lower, upper = upper))

\dontshow{FinalizeExample()}

}


