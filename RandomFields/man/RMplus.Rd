\name{RMplus}
\alias{RMplus}
\alias{RMmodelplus}
\alias{RMtrendplus}
\alias{+}
\alias{RM_PLUS}
\title{Addition of Random Field Models}
\description{
 \command{\link{RMplus}} is an additive covariance model which depends on
 up to 10 submodels \eqn{C_0, C_1, ..., C_9}.
 In general, realizations of the created \command{\link{RMmodel}} are pointwise
 sums of independent realizations of the submodels.
 
 In particular, if all submodels are given through a covariance
 function, the resulting model is defined through its covariance
 function, which is the sum of the submodels' covariances.
 Analogously, if all submodels are given through a variogram.
}
\usage{
RMplus(C0, C1, C2, C3, C4, C5, C6, C7, C8, C9, trend,
       var, scale, Aniso, proj)
RMtrendplus(C0, C1, C2, C3, C4, C5, C6, C7, C8, C9, add.na=FALSE,
            var, scale, Aniso, proj)
RMmodelplus(C0, C1, C2, C3, C4, C5, C6, C7, C8, C9, 
            var, scale, Aniso, proj)
}
\arguments{
 \item{C0}{an \command{\link{RMmodel}}. See also  Details.}
 \item{C1,C2,C3,C4,C5,C6,C7,C8,C9}{optional; each an
   \command{\link{RMmodel}}.}
 \item{trend}{logical. If \code{TRUE} all summands are considered as
   trend (or shape) functions. If \code{FALSE} all summands are considered
   as part of the definition of the correlation structure. If not given,
   each summand is attributed separately by the package.
 }
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
   \command{\link{RMmodel}}. If not passed, the above
   model remains unmodified.}
 \item{add.na}{
   logical. 
   If \code{TRUE} the submodels
   are multiplied with \code{NAs} to estimate the parameters of the
   corresponding linear model. Note that
   \code{NAs} is only added if no multiplicative constants are involved
   in the respective submodel.
 }
}
\details{
 \command{Usually, \link{RMmodel}}s are summed up via the
 \code{+}-operator, e.g. C0 + C1, which is indeed only a short cut
 for using \command{RMplus}. \command{RMplus} shows a greater flexibility.
 In particular,
 the global arguments \code{var,scale,Aniso,proj} of \command{\link{RMplus}}
 applying to all submodels, e.g.\ usefull for a global variance parameter.

 \code{RMtrendplus} and  \code{RMmodelplus} have a two-fold purpose:
 (i) they can be used clarify whether a model defines a trend/shape or
 a definite function (ii) the arguments of 
 \code{RMtrendplus} and \code{RMmodelplus} are more flexible,
 allowing for both formulae and the most flexible list notation for models,
 see the Notes in \link{RMmodelsAdvanced}.

 Note, that \command{RMplus}, \code{RMtrendplus} and
 \code{RMmodelplus} also allow only one argument, in contrast to the
 binary operator \code{+}.
}
\value{
 \command{\link{RMplus}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.
}


\me
\seealso{
 \command{\link{RMmult}},
 \command{\link{RMmodel}},
 \command{\link{RMsum}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.
}
\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

model <- RMplus(RMgauss(), RMnugget(), var=0.5) ## common variance parameter
model2<- RMgauss() + RMnugget(var=0.1)
plot(model, "model.+"=model2, type=c("p", "l"), pch=20, xlim=c(0,3)) # the same

\dontshow{FinalizeExample()}}

\keyword{spatial}
\keyword{models}
