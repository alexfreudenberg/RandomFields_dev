\name{RFplot}
\alias{RFplot}
\title{Plotting Models, Simulations, Empirical Statistics}
\description{
  \command{RFplot} plot various outputs of the package.
    It integrates the functions \command{\link{RFplotSimulation}}, 
    \command{\link{RFplotEmpVariogram}},
    \command{\link{RFplotModel}}.
    
  In many cases, \command{\link[base]{plot}} can be used instead of the
  command \command{RFplot}. The use of \command{RFplot} is mandatory,
    if \code{\link{RMmodels}} are defined by means of an
  \code{\link[=RFformula]{formula}}. 
}
\usage{
RFplot(x, \ldots)
}
\arguments{
  \item{x, ...}{see the commands \command{\link{RFplotSimulation}}, 
    \command{\link{RFplotEmpVariogram}},
    \command{\link{RFplotModel}}}
  for the significance of the arguments.
}
\details{
  If \command{RFplot} is called without argument \code{x},
  a list of possible specification of the arguments is given.
}

\value{
See \command{\link{RFplotSimulation}}, 
    \command{\link{RFplotEmpVariogram}},
    \command{\link{RFplotModel}}
}

\me

\seealso{
  \command{\link[base]{plot}},
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFplotEmpVariogram}},
}

\examples{\dontshow{StartExample(reduced=FALSE)}
model <- ~RMexp(s=S)
RFplot(model, 1:10, params=list(S=5))
\dontshow{FinalizeExample()}}
\keyword{spatial}






