\name{RFformula}
\alias{RFformula}
\alias{RMformula}

\title{RFformula - syntax to design random field models with trend
 or linear mixed models}

\description{It is described how to create a formula, which, for example, can be used as an argument of \command{\link{RFsimulate}} and
 \command{\link{RFfit}} to simulate and to fit data according to the
 model described by the formula.

 In general, the created formula serves two purposes:
 \itemize{
 \item to describe models in the \dQuote{Linear Mixed
 Models}-framework %including fixed and random effects
 \item to define models for random fields including trend surfaces from
 a geostatistical point of view.
 }

 Thereby, fixed effects and trend surfaces can be adressed via
 the expression \command{\link{RMfixed}} and the function
 \command{\link{RMshape}}. In simple cases, the trend can also
 be given in a very simple, see the examples below.
 The covariance structures of the zero-mean
 multivariate normally distributed %random effects and
 random field
 components are adressed by objects of class \code{\link[=RMmodel-class]{RMmodel}}, which
 allow for a very flexible covariance specification.


 Note: variables, \code{V} say, defined by the user, e.g. on the global
 environment, may not be used in formulae within \code{params} as
 the formulae are evaluated in a completely different environment,
 use local a variable,\code{W} say, instead:
\code{ V <- 5}\cr
\code{ params = list(X = ~ 5 + W)}\cr
\code{RFfit(params=params, W=V)} \cr
 
 Note that instead of \code{W} any other name can be used, even \code{V}
 itself, i.e.,
\code{params = list(X = ~ 5 + V)}\cr
\code{RFfit(params=params, V=V)}\cr

 
 See \link{RFformulaAdvanced} for rather complicated model definitions.
}

 
\note{
 Models given by formulae can be plotted only by \link{RFplot}.
}

\details{ 
 Define all formulae with \code{params}. Within the \code{RMmodel}, use
  only assignments such as \code{arg=abc}, where \code{abc} is defined
  in \code{params}. External variables can be used 

 Any definition with NAs needs an exact counterpart in the model, e.g.
 \code{params <- list(V = NA)}
 needs somewhere in the model a simple declaration (and exactly one of them)
 E.g., \code{RMmodel(var = V)} is correct; but
 \code{RMmodel(sc=V, var=V)},  \code{RMmodel(var=2 * V)}
 is wrong. Also not using V at all is wrong. However,
 \code{params <- list(V = NA, S = ~V);}
 \code{RMmodel1(sc=S, var=V)} is correct.
 If an introduced variable, X say, does not have a counterpart in the
 model, e.g. \code{params <- list(S=sin(X), V=X^2, X = NA);} together with
 \code{RMmodel(sc=V,  var=V)}, this variable must be \command{RMdeclare}d, e.g.
 \code{RMmodel(sc=V, var=V) + RMdeclare(X=X)}.
 

 The formula should be of the type
 \deqn{response ~ fixed effects %+ random effects
   + error term}
 or 
 \deqn{response ~ trend + zero-mean random field + nugget effect,}
 respectively.
 
 Thereby:
 
 \itemize{
 \item response\cr
 optional; name of response variable

 \item trend:\cr
 usual definition of trend.


 \item random effects: not programmed.
% optional, should be a sum (using \command{\link[=RMplus]{+}})%%check link
% of components either of the form \code{X@RMfixed(beta)} or
% \code{\link{RMshape}(...)} with \eqn{X} being a design matrix 
% and \eqn{\beta} being a vector of coefficients (see
% \command{\link{RMfixed}} and \command{\link{RMshape}}).\cr
% Note that a fixed effect of the form \eqn{X} is interpreted as
% \code{X@RMfixed(beta=NA)} by default (and \eqn{\beta} is estimated
% provided that the formula is used in \command{\link{RFfit}}).
 
% \item random effects/zero-mean random field:\cr
% optional, should be a sum (using \command{\link[=RMplus]{+}})%%check link
% of components of the form \code{Z@model}
%  where \eqn{Z} is a design matrix and \code{model} is an object of
%class \code{\link[=RMmodel-class]{RMmodel}}.\cr
%\code{Z@model} describes a vector of random effects which is
% normally distributed with zero mean and covariance matrix \eqn{Z
% \Sigma Z^T} where \eqn{Z^T} is the transpose of \eqn{Z} and
% \eqn{\Sigma} is the covariance matrix according to \code{model}.\cr
% Note that a random effect/random fluctuation of the form
% \code{model} is viewed as \code{I@model} where \eqn{I} is the
% identity matrix of corresponding dimension.
 
 \item error term/nugget effect\cr
 optional, should be of the form \code{\link{RMnugget}(...)}.
 \command{\link{RMnugget}} describes a vector
 of iid Gaussian random variables.

% Please note that the character \dQuote{@} in the RFformula-context can only
% be used to multiply design-matrices with corresponding vectors of
% fixed or random effects, whereas in the context of S4-classes \dQuote{@} is
% used to access slots of corresponding objects.
 }
 
}

\section{IMPORTANT}{
  Note that in formula constants within a canonic trend expression
  are interpreted as part of a linear
  model, i.e. the corresponding parameter has to be estimated
  (e.g. \code{~ 1 + ...}) whereas in models not given as formula this
  parameter must be given explicitly by \code{NA} to be estimated.

  Note also that \code{Aniso} may not be used as an argument if parts of
  \code{Aniso} are to be estimated. Instead, \code{anisoT} must be used.
}

\note{ %to do: make the following point cleares  
  (additional) argument names should always start with a capital
  letter. Small initial letters are reserved  for \command{\link{RFoptions}}.
}

\references{
 \itemize{
 \item Chiles, J.-P. and P. Delfiner (1999) \emph{Geostatistics. Modeling
 Spatial Uncertainty.} New York, Chichester: John Wiley & Sons.
 
 \item McCulloch, C. E., Searle, S. R. and Neuhaus, J. M. (2008)
 \emph{Generalized, linear, and mixed models.} Hoboken, NJ: John
 Wiley & Sons.
 
 \item Ruppert, D. and Wand, M. P. and Carroll, R. J. (2003)
 \emph{Semiparametric regression.} Cambridge: Cambridge University
 Press.
 }
}

\me

\seealso{
  \link{RFformulaAdvanced},
  \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}},
% \command{\link{RFplot}},
\code{\link[=RandomFields-package]{RandomFields}}.
}

\keyword{spatial}

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

RFoptions(modus_operandi="sloppy")

##############################################################
#
# Example  : Simulation and fitting of a two-dimensional
# Gaussian random field with exponential covariance function
#
###############################################################

V <- 10
s <- 0.3
M <- 3
model <- RMexp(var=V, scale=s) + M
x <- y <- seq(1, 3, 0.1)

simulated <- RFsimulate(model = model, x=x, y=y)
plot(simulated)


# an alternative code to the above code:
model <- ~ Mean + RMexp(var=Var, scale=Sc)
simulated2 <- RFsimulate(model = model,x=x, y=y, Var=V, Sc=s, Mean=M)
plot(simulated2)


# a third way of specifying the model using the argument 'param'
# the initials of the variables do not be captical letters
model <- ~ M + RMexp(var=var, scale=sc)
simulated3 <- RFsimulate(model = model,x=x, y=y,
                         param=list(var=V, sc=s, M=M))
plot(simulated3)


# Estimate parameters of underlying covariance function via
# maximum likelihood
model.na <- ~ NA + RMexp(var=NA, scale=NA)
fitted <- RFfit(model=model.na, data=simulated)

# compare sample mean of data with ML estimate, which is very similar:
mean(simulated@data[,1]) 
fitted




###########################################################
#
# Example 2: Simulate and fit a geostatistical model
#
###########################################################

# Simulate a Gaussian random field with trend m((x,y)) = 2 + 1.5 x - 3 y
# with vector of coordinates (x,y)
# and covariance function C(s,t) = 3*exp(-||(2*(s_1-t_1),s_2-t_2)||) +
# 1.5*exp(-||(2*(s_1-t_1),s_2-t_2)||^2)
# for s=(s_1,s_2) and t=(t_1,t_2)

C <- RMplus(RMexp(var=3), RMgauss(var=1.5), Aniso=matrix(c(2,0,0,1),nrow=2))
model <- ~ 2 + 1.5 * x - 3 * y + C ## tilde necessary
simulated <- RFsimulate(model=model, x=seq(0,2,0.1), y=seq(-1,3,0.2))


## equivalent model
C <- RMS(RMexp(var=3) + RMgauss(var=1.5), anisoT=matrix(c(2,0,0,1),nrow=2))
model <- RMshape(RMpolynome(1, value=c(2, 1.5, -3), dim=2)) + C # ~ optional
simulated2 <- RFsimulate(model=model, x=seq(0,2,0.1), y=seq(-1,3,0.2))
stopifnot(all.equal(simulated, simulated2))

simulated <- RFsimulate(model=model, x=seq(0,2,0.1), y=seq(-1,3,0.2))


# Estimate trend (polynomial of degree 1) and variance components such
# that var_exp = 2*var_gauss as in the true model used for simulation
Cest <- RMS(RMexp(var=2) + RMgauss(var=2),
               var=NA, anisoT=matrix(c(NA,0,0,NA),nrow=2))
model.na <- ~  RMpolynome(1, dim=2) + Cest

fit <- RFfit(model=model.na, data=simulated)


\dontshow{RFoptions(modus_operandi="normal")}

\dontshow{FinalizeExample()}}