\name{RFfitoptimiser}
\alias{RFfitoptimiser}
\alias{RFfitOptimiser}
\title{Optimisers for
  fitting model parameters to spatial data
}
\description{
  See \command{\link{RFfit}} for a detailed description of
  the fitting procedure.
}

\details{
  Two parameters, see also \command{\link{RFoptions}}
  can be passed to \command{\link{RFfit}} that allow
  for choosing an optimiser different from \command{\link[stats]{optim}}:

  \code{optimiser} takes one of the values
  \code{"optim"}, \code{"optimx"}, \code{"soma"}, \code{"nloptr"},
  \code{"GenSA"}, \code{"minqa"}, \code{"pso"} or \code{"DEoptim"},
  see the corresponding packages for a description.

  If \code{optimiser="nloptr"}, then the additional parameter
  \code{algorithm} must be given which takes the values
  \code{"NLOPT_GN_DIRECT"},  
  \code{"NLOPT_GN_DIRECT_L"},
  \code{"NLOPT_GN_DIRECT_L_RAND"}, 
  \code{"NLOPT_GN_DIRECT_NOSCAL"},
  \code{"NLOPT_GN_DIRECT_L_NOSCAL"}, 
  \code{"NLOPT_GN_DIRECT_L_RAND_NOSCAL"},
  \code{"NLOPT_GN_ORIG_DIRECT"}, 
  \code{"NLOPT_GN_ORIG_DIRECT_L"},
  \code{"NLOPT_LN_PRAXIS"}, 
  \code{"NLOPT_GN_CRS2_LM"}, 
  \code{"NLOPT_LN_COBYLA"},
  \code{"NLOPT_LN_NELDERMEAD"}, 
  \code{"NLOPT_LN_SBPLX"},
  \code{"NLOPT_LN_BOBYQA"},
  \code{"NLOPT_GN_ISRES"},
  see \pkg{nloptr} for a description.
    
}

\me

\seealso{
  \command{\link{RFfit}}, \command{\link{RFoptions}}
}
\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

data(soil)
RFoptions(coordidx=1:2)

if (interactive()) {
   model <- moisture ~ 1 + RMwhittle(scale=NA, var=NA, nu=NA) +
                       RMnugget(var=NA)
} else model <- moisture ~ 1 + RMnugget(var=NA)

## standard optimiser 'optim'
print(system.time(fit <- RFfit(model, data=soil)))
print(fit)

opt <- "optimx" #  30 sec; better result
print(system.time(fit2 <- try(RFfit(model, data=soil, optimiser=opt))))%## OK
print(fit2)

\dontshow{\dontrun{ %in 6 sep 2017, soma has leak problems
opt <- "soma" #  450 sec % sehr schlecht
print(system.time(fit2 <- try(RFfit(model, data=soil, optimiser=opt))))%## OK
print(fit2)
}}

opt <- "minqa" # 330 sec %ok
print(system.time(fit2 <- try(RFfit(model, data=soil, optimiser=opt))))%## OK
print(fit2)


opt <- "nloptr"
algorithm <- RC_NLOPTR_NAMES
\dontshow{if(!interactive()) algorithm <- RC_NLOPTR_NAMES[1]}
for (i in 1:length(algorithm)) { % all algorithms are fast
  print(algorithm[i])
  print(system.time(fit2 <- try(RFfit(model, data=soil, optimiser=opt,%## OK
                                    algorithm=algorithm[i]))))
  print(fit2)
}


if (interactive()) {
## the following two optimisers are too slow to be run on CRAN.

opt <- "pso" # 600 sec
print(system.time(fit2 <- try(RFfit(model, data=soil, optimiser=opt))))%## OK
print(fit2)

opt <- "GenSA" #  10^4 sec
print(system.time(fit2 <- try(RFfit(model, data=soil, optimiser=opt))))%## OK
print(fit2)
}


\dontshow{RFoptions(modus_operandi="normal")}
\dontshow{FinalizeExample()}}

\keyword{spatial}
\keyword{optimize}

