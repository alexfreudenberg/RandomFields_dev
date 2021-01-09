
## Authors 
## Martin Schlather, schlather@math.uni-mannheim.de
##
##
## Copyright (C) 2012 -- 2017  Martin Schlather
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 3
## of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


## Version3
#RMstrokorbMono <- function(phi) stop("Please use 'RMm2r' instead")
#RMstrokorbBall <- function(phi) stop("Please use 'RMm3b' instead")
#RMstrokorbPoly <- function(phi) stop("Please use 'RMmps' instead")




PrepareModelOld <-  function(model, param, trend=NULL, 
                          nugget.remove=TRUE, method=NULL) {
  ## any of the users model definition (standard, nested, list) for the
  ## covariance function is transformed into a standard format, used
  ## especially in the c programs
  ##
  ## overwrites in some situation the simulation method for nugget.
  ## allows trend to be NA (or any other non finite value  -- is not checked!)
  ## trend has not been implemented yet!

  if (isS4(model))
    stop("S4 models cannot be combined with obsolete RandomFields functions")

  if (!is.null(method)) stop("to give method in suPrepareModelOld is obsolete")
  
  if (!is.null(trend))      
     if (!is.numeric(trend) || length(trend)!=1)
        stop("in the obsolete setting, only constant mean can used")
 
  if (is.list(model) && is.character(model[[1]]) &&
      (is.null(names(model)) || names(model)[[1]]=="")) {
     if (!missing(param) && !is.null(param))
        stop("param cannot be given in the extended definition")
     if (is.null(trend)) return(model)
     trend <- list(RM_TREND[2], mean=trend)
     if (model[[1]] %in% RM_PLUS) return(c(model, list(trend)))
     else return(list(SYMBOL_PLUS, model, trend))
  }
    
  printlevel <- RFoptions()$basic$printlevel
  STOP <- function(txt) {
    if (printlevel>=PL_ERRORS) {
      cat("model: ")
      if (!missing.model) Print(model) else cat(" missing.\n") #
      cat("param: ")
      if (!missing.param) Print(param) else cat(" missing.\n") #
      cat("trend: ")
      Print(trend) # 
     }
    stop("(in PrepareModelOld) ", txt, call.=FALSE)
  }
   
  transform <- function(model) {
    if (!is.list(model)) {
      STOP("some elements of the model definition are not lists")
    }
    m <- list(DOLLAR[1], var=model$v)
    lm <- length(model) - 3 # var, scale/aniso, name
    if (!is.null(model$a)) m$aniso <- model$a else m$scale <- model$scale
##    model <- c(model, if (!is.null(model$a))
##               list(aniso=model$a) else list(scale=model$s)) ## ???

    
    if (!is.na(p <- pmatch("meth", names(model), duplicates.ok=TRUE))) {
      if (printlevel>=PL_ERRORS)  Print(p, model) #
      stop("method cannot be given with the model anymore. It must be given as a parameter to the function. See 'RFoptions' and 'RFsimulate'")
     }
  
    if (!is.null(model$me))
      stop("'mean' seems to be given within the inner model definitions"); 
    if (!is.character(model$m)) {
       stop("'model' was not given extacly once each odd number of list entries or additional unused list elements are given.")
    }
    m1 <- list(model$m)
    if (!is.null(model$k)) {
      lm <- lm - 1
      if (length(model$k) != 0)
        for (i in 1:length(model$k)) {
          eval(parse(text=paste("m1$k", i, " <- model$k[", i, "]", sep="")))
      }
    }
    if (lm != 0) {
      if (printlevel>=PL_ERRORS) Print(lm, model) #
      stop("some parameters do not fit")
    }
    m <- c(m, list(m1))

    return(m)
    
  } # end transform

  op.list <- c(SYMBOL_PLUS, SYMBOL_MULT)  ## if others use complex list definition !
  missing.model <- missing(model)
  missing.param <- missing(param) || is.null(param)

  if (missing.param && is.null(model$param)) { ## full model
    if (RFoptions()$internal$warn_oldstyle)
      warning("the sequential list format is depreciated.")
    if (missing.model || (length(model)==0)) model <- list()
    else if (!is.list(model))
      STOP("if param is missing, model must be a list of lists (or a list in the extended notation)")
    if (is.null(trend) + is.null(model$mean) + is.null(model$trend)<2)
      STOP("trend/mean is given twice")
    if (!is.null(model$mean)) trend <- model$mean else
    if (!is.null(model$trend)) trend <- model$trend else trend <- NULL

    model$trend <- model$mean <- NULL
    ## the definition might be given at a deeper level as element
    ## $model of the list:
    if (is.list(model$model)) {
      if (!is.list(model$model[[1]]))
        STOP("if param is missing, the model$model must be a list of lists")
      model <- model$model
    }
    if (length(model)==0) { ## deterministic      
      return(if (is.null(trend)) NULL else list(RM_TREND[2], mean=trend))
    }
    if (length(model) %% 2 !=1) STOP("list for model definition should be odd")
    if (length(model)==1)
      return(if (is.null(trend) ||
                 is.numeric(trend) && length(trend)==1 && !is.na(trend)&&trend==0)
             transform(model[[1]]) 
             else list(SYMBOL_PLUS, transform(model[[1]]),
                       list(RM_TREND[2], mean=trend)));

    op <- pmatch(c(model[seq(2, length(model), 2)], recursive=TRUE),
                 op.list, duplicates.ok=TRUE) - 1
    if (!all(is.finite(op))) STOP("operators are not all allowed; see the extended list definition for extensions")
    model <- model[seq(1, length(model), 2)]

    plus <- which(op==0)
    if (length(plus) == 0) {
      m <- list("*", lapply(model, transform))
    } else {
      plus <- c(0, plus, length(op)+1)
      m <- list(SYMBOL_PLUS)
      for (i in 1:(length(plus) - 1)) {
        m[[i+1]] <-
          if (plus[i] + 1 == plus[i+1]) transform(model[[plus[i] + 1]])
          else list(SYMBOL_MULT,
                    lapply(model[(plus[i] + 1) : plus[i+1]], transform))
      }
    }
   model <- m
  } else { ## standard definition or nested model
    if (missing.param) { ## a simple list of the model and the
      ##                    parameters is also possible
      if (is.null(param <- model$p)) STOP("is.null(model$param)")
      stopifnot(is.null(trend) || is.null(model$trend))
      if (is.null(trend)) trend <- model$trend
      if (!is.null(model$mean)) {
        if (!is.null(trend)) STOP("mean and trend given twice")
        trend <- model$mean
      }
      model <- model$model
    }
    stopifnot(is.character(model), length(model)==1)
    if (is.matrix(param)) { ## nested
      if (nrow(param) == 1)
        return(PrepareModelOld(model=model, param=c(param[1], 0, param[-1]),
                            trend=trend))
      name <- model
      model <- list(SYMBOL_PLUS)#, method=method)
      for (i in 1:nrow(param)) {
        model <- c(model,
                   if (is.na(param[i, 2]) || param[i, 2] != 0)
                   list(list(DOLLAR[1], var=param[i, 1], scale=param[i, 2],
                             if (ncol(param) >2) list(name, k=param[i,-1:-2])
                             else list(name)))
                   else list(list(DOLLAR[1], var=param[i,1],
                                  list(RM_NUGGET[2]))))
      }
    } else if (is.vector(param)) {  ## standard, simple way
      ## falls trend gegeben, dann ist param um 1 Komponente gekuerzt
      if (is.null(trend)) {
        trend <- param[1]
        param <- param[-1]
      } else message("It is assumed that no mean is given so that the first component of param is the variance")
      if (model == RM_NUGGET[2]) {
        model <- transform(list(model=model, var=sum(param[1:2]), scale=1)) 
      } else {
        if  (length(param) > 3)
          model <- transform(list(model=model, var=param[1], scale=param[3],
                                  k=param[-1:-3]))
        else 
          model <- transform(list(model=model, var=param[1], scale=param[3]))
        if (is.na(param[2]) || param[2] != 0 || !nugget.remove) {# nugget
          model <- list(SYMBOL_PLUS,
                        model,
                        transform(list(model=RM_NUGGET[2], var=param[2], scale=1)))
        }
        ## if (!is.null(method)) model <- c(model, method=method) ## doppelt
      }
    } else stop("unknown format")  # end nested/standard definition
  }

  return(if (is.null(trend) ||
             is.numeric(trend) && length(trend)==1 &&  !is.na(trend) &&trend==0)
         return(model)
         else if (model[[1]] %in% RM_PLUS)
                 c(model, list(list(RM_TREND[2], mean=trend)))
         else list(SYMBOL_PLUS, model, list(RM_TREND[2], mean=trend)))
}



RFempiricalvariogram <- function(...) {
  obsolete <- as.character(list(match.call())[[1]])[1]
  new <- "RFvariogram"
  warning("'", obsolete, "' is an obsolete function.\n  Use '", new,
          "' which calculates both empirical and theoretical values.\n  (later versions might also be able to fit)")
  do.call(new, list(...))
}

RFempiricalmadogram <- function(...) {
  obsolete <- as.character(list(match.call())[[1]])[1]
  new <- "RFmadogram"
  warning("'", obsolete, "' is an obsolete function.\n  Use '", new,
          "' which calculates both empirical and theoretical values.\n  (later versions might also be able to fit)")
  do.call(new, list(...))
}


RFempiricalcovariance <- function(...) {
  obsolete <- as.character(list(match.call())[[1]])[1]
  new <- "RFcovariance"
  warning("'", obsolete, "' is an obsolete function.\n  Use '", new,
          "' which calculates both empirical and theoretical values.\n  (later versions might also be able to fit)")
  do.call(new, list(...))
}




  #### Version 2

RFoldstyle <- function(old=TRUE) {
  RFoptions(general.spConform = !old,
            internal.warn_newstyle = old,
            internal.warn_oldstyle = !old)
  invisible(NULL)
}

OldDistribution <- function(distribution) {
  if (missing(distribution) || is.null(distribution) || is.na(distribution))
    return(NULL)
  names <- c("Gauss", "MaxStable")
  Distr <- list("RMgauss", NULL) ## obsolete
  nr <- pmatch(distribution, names)
  if (!is.finite(nr)) stop("unknown method")
  return(Distr[[nr]])
}

OldMethod <- function(method) {
  if (missing(method) || is.null(method) || is.na(method)) return("RPgauss")
  names <-  c("circulant embed", 
              "cutoff CE",
              "intrinsic CE",
              "TBM2", 
              "TBM3",
              "spectral TBM", 
              "direct decomp.",
              "sequential",
              "Markov",
              "average",
              "nugget", 
              "coins",
              "hyperplanes",
              ##
              "max.MPP",
              "extremalGauss",
              "BrownResnick",
              "BooleanFunction")
  Meth <- c("circulant", "cutoff", "intrinsic", "tbm", "tbm",
            "spectral", "direct",  "sequential",
            "Markov does not work anymore",
            "average", "nugget", "coins", "hyperplane",
            "smith", "schlather", "brownresnick", "schlather")
  nr <- pmatch(method, names)
  if (!is.finite(nr)) stop("unknown method")
  return(paste("RP", Meth[nr], sep=""))
}


Variogram <- function(x,
                      model, param=NULL,
                      dim=NULL,
                      ## y=NULL,
                      Distances) {
  if (RFoptions()$internal$warn_oldstyle)
    warning("The function is obsolete. Use 'RFvariogram' instead.")
  if (is.null(dim) &&!missing(Distances))
    dim <- if (is.matrix(x)) ncol(x) else 1
  if (!missing(param) && is.na(param[1])) param[1] <- 0

##  Print(model, param)
  model <- PrepareModelOld(model, param)
  RFvariogram(x=x, model=model, dim=dim, distances=Distances)
}


Covariance <- CovarianceFct <-
  function(x, y=NULL,
           model, param=NULL,
           dim=NULL,
           Distances, fctncall=c("Covariance", "Variogram", "CovMatrix")) {
    if (RFoptions()$internal$warn_oldstyle)
      warning("The function is obsolete. Use 'RFcovariance' instead")
    if (is.null(dim) && !missing(Distances))
      dim <- if (is.matrix(x)) ncol(x) else 1
    if (!missing(param) && is.na(param[1])) param[1] <- 0
    model <- PrepareModelOld(model, param)
    fctncall <- match.arg(fctncall)
    rfeval(x=x, y=y, model=model, dim=dim, distances=Distances,fctncall=fctncall)
  }

CovMatrix <- function(x, y=NULL,
           model, param=NULL,
           dim=NULL,
           Distances) {
  if (RFoptions()$internal$warn_oldstyle)
    warning("The function is obsolete. Use 'RFcovmatrix' instead")
   if (is.null(dim) && !missing(Distances))
      dim <- if (is.matrix(x)) ncol(x) else 1
 if (!missing(param) && is.na(param[1])) param[1] <- 0
  model <- PrepareModelOld(model, param)
   RFcovmatrix(x=x, y=y, model=model, dim=dim, distances=Distances)
}

DoSimulateRF <- function (n = 1, register = 0, paired=FALSE, trend=NULL) {
  if (RFoptOld[[1]]$internal$warn_oldstyle)
    warning("The function is obsolete.\nSee documentation of 'RFSimulate' (advanced part) instead.")
 if (!is.null(trend))
    stop("parameter trend in DoSimulateRF cannot be used anymore")
  RFoptOld <-
    internal.rfoptions(register=register, gauss.paired=paired, spConform=FALSE)
  on.exit(RFoptions(LIST=RFoptOld[[1]]))


   RFsimulate(n=n)
}



InitSimulateRF <- function (x, y = NULL, z = NULL, T=NULL,
                            grid=!missing(gridtriple), model, param,
                            trend,  method = NULL,
                            register = 0, gridtriple, distribution=NA) {
  if (RFoptOld[[2]]$internal$warn_oldstyle)
    warning("This function is obsolete. Use 'RFsimulate' instead.")
 RFoptOld <-
    internal.rfoptions(register=register, #gauss.method=method,
                       spConform=FALSE)
  on.exit(RFoptions(LIST=RFoptOld[[1]]))

 
  OldDistribution(distribution) ## only check whether input is oK
  
  model <- c(list(OldMethod(method)),
             list(PrepareModelOld(model, param, trend)))

  if (!missing(gridtriple)) {
    if (gridtriple) {
      if (!is.null(x)) {
        if (is.matrix(x)) {
          x <- apply(x, 2,
                     function(s) {c(s[1], s[3], length(seq(s[1], s[2], s[3])))})
        } else {
          x <- seq(x[1], x[2], x[3])
          if (length(y)!=0)  y <- seq(y[1], y[2], y[3])
          if (length(z)!=0)  z <- seq(z[1], z[2], z[3])
        }
      }
    }
  }
  if (length(T)!=0)  T <- c(T[1], T[3], length(seq(T[1], T[2], T[3])))

  p <- list("Simulate", PrepareModel2(model, return_transform=FALSE)$model)
  warn.seed.not.na(RFoptOld, TRUE)
  rfInit(model=p, x=x, y =y, z = z, T=T, grid=grid, RFopt=RFoptOld[[2]])
}



InitGaussRF <- function(x, y = NULL, z = NULL, T=NULL,
                        grid = !missing(gridtriple),
                        model, param,
                        trend=NULL, method = NULL,
                        register = 0, gridtriple) {
  InitSimulateRF(x=x, y=y, z=z,  T=T, 
                        grid=grid, model=model, param=param, trend=trend,
                        method=method,
                        register=register,
                        gridtriple=gridtriple, distribution="Gauss")
}


GaussRF <- function (x, y = NULL, z = NULL, T=NULL,
          grid=!missing(gridtriple), model, param, trend=NULL, method = NULL, 
          n = 1, register = 0, gridtriple,
          paired=FALSE, PrintLevel=1, Storing=TRUE, ...) { #
  if (RFoptions()$internal$warn_oldstyle)
    warning("The function is obsolete. Use 'RFsimulate' instead.")


  model <- c(list(OldMethod(method)),
             list(PrepareModelOld(model, param, trend)))

#  Print(model)

 
  if (!missing(gridtriple)) {
    if (gridtriple) {
      if (!is.null(x)) {
        if (is.matrix(x)) {
          x <- apply(x, 2,
                     function(s) {c(s[1], s[3], length(seq(s[1], s[2], s[3])))})
        } else {
          x <- seq(x[1], x[2], x[3])
          if (length(y)!=0)  y <- seq(y[1], y[2], y[3])
          if (length(z)!=0)  z <- seq(z[1], z[2], z[3])
        }
      }
    }
  }
  if (length(T)!=0)  T <- c(T[1], T[3], length(seq(T[1], T[2], T[3])))

  RFoptOld <-
    internal.rfoptions(register=register, gauss.paired=paired,
                      # gauss.method=method,
                       spConform=FALSE)
  on.exit(RFoptions(LIST=RFoptOld[[1]]))

  
  RFsimulate(model=model, x=x, y=y, z=z, T=T, grid=grid, n=n,
             #gauss.method=method,
             register=register, gauss.paired=paired,
             printlevel = PrintLevel, storing=Storing,
             ...)
  
 ##   str(RFoptions())
 
}




## it does not make sense to me at the moment that a space-time model
## for extremes is defined.

InitMaxStableRF <- function(x, y = NULL, z = NULL, grid=NULL, model, param,
                            maxstable,
                            method = NULL,
                            register = 0, gridtriple = FALSE) {
  RFoptOld <-
    internal.rfoptions(register=register, #gauss.method=method,
                       spConform=FALSE)
  on.exit(RFoptions(LIST=RFoptOld[[1]]))
  if (RFoptOld[[2]]$internal$warn_oldstyle)
    warning("The function is obsolete. Use 'RFsimulate' instead.")
  meth <- if (!is.null(method)) method else maxstable
  if (is.null(meth)) stop("method not given")
  model <- c(list(OldMethod(meth)),
             list(PrepareModelOld(model, param)))

  if (!missing(gridtriple)) {
    if (gridtriple) {
      if (!is.null(x)) {
        if (is.matrix(x)) {
          x <- apply(x, 2,
                     function(s) {c(s[1], s[3], length(seq(s[1], s[2], s[3])))})
        } else {
          x <- seq(x[1], x[2], x[3])
          if (length(y)!=0)  y <- seq(y[1], y[2], y[3])
          if (length(z)!=0)  z <- seq(z[1], z[2], z[3])
         }
      }
    }
  }
  if (length(T)!=0)  T <- c(T[1], T[3], length(seq(T[1], T[2], T[3])))
  
  p <- list("Simulate", PrepareModel2(model, return_transform=FALSE)$model)
  warn.seed.not.na(RFoptOld, TRUE)
  rfInit(model=p, x=x, y=y, z=z, grid=grid, RFopt=RFoptOld[[2]])
}

  
MaxStableRF <- function (x, y = NULL, z = NULL, grid=NULL,
                         model, param, maxstable,
                         method = NULL,
                         n = 1, register = 0,
                         gridtriple = FALSE,
                         ...) {
  RFoptOld <-
   if (n>1) internal.rfoptions(..., register=register, #gauss.method=method,
                               spConform=FALSE, storing=TRUE)
   else internal.rfoptions(..., register=register, #gauss.method=method,
                           spConform=FALSE)
  on.exit(RFoptions(LIST=RFoptOld[[1]]))

  if (RFoptOld[[2]]$internal$warn_oldstyle)
    warning("The function is obsolete. Use 'RFsimulate' instead.")
  meth <- if (!is.null(method)) method else maxstable
  if (is.null(meth)) stop("method not given")
  model <- c(list(OldMethod(meth)),
             list(PrepareModelOld(model, param)))

  if (!missing(gridtriple)) {
    if (gridtriple) {
      if (!is.null(x)) {
        if (is.matrix(x)) {
          x <- apply(x, 2,
                     function(s) {c(s[1], s[3], length(seq(s[1], s[2], s[3])))})
        } else {
          x <- seq(x[1], x[2], x[3])
          if (length(y)!=0)  y <- seq(y[1], y[2], y[3])
          if (length(z)!=0)  z <- seq(z[1], z[2], z[3])
        }
      }
    }
  }
 
  return(RFsimulate(model=model, x=x, y=y, z=z, grid=grid, n=n))
}


EmpiricalVariogram <-
  function (x, y = NULL, z = NULL, T=NULL, data, grid=NULL, bins, gridtriple = FALSE,
            phi,  ## phi[1] erste richtung, phi[2] : anzahl der richtungen
            theta, ## aehnlich
            deltaT) {
    if (RFoptions()$internal$warn_oldstyle)
      warning("This function is obsolete. Use RFvariogram instead.")
 
  if (!missing(gridtriple)) {
    if (gridtriple) {
      if (!is.null(x)) {
        if (is.matrix(x)) {
          x <- apply(x, 2,
                     function(s) {c(s[1], s[3], length(seq(s[1], s[2], s[3])))})
        } else {
          x <- seq(x[1], x[2], x[3])
          if (length(y)!=0)  y <- seq(y[1], y[2], y[3])
          if (length(z)!=0)  z <- seq(z[1], z[2], z[3])
        }
      }
    }
  }
  if (length(T)!=0)  T <- c(T[1], T[3], length(seq(T[1], T[2], T[3])))
  
  rfempirical(x=x, y=y, z=z, T=T, data=data, grid=grid,
              bins=bins, nphi=phi, ntheta=theta, deltaT=deltaT, vdim=1)
}



Kriging <- function(krige.method, x, y=NULL, z=NULL, T=NULL,
                    grid=NULL, gridtriple=FALSE,
                    model, param, given, data, trend=NULL,            
                    pch=".", return.variance=FALSE,
                    allowdistanceZero = FALSE, cholesky=FALSE) {
                    
  RFoptOld <- internal.rfoptions(general.pch=pch, spConform=FALSE,
                                 return_variance=return.variance,
                                 allowdistanceZero=allowdistanceZero)
  on.exit(RFoptions(LIST=RFoptOld[[1]]))
  
   if (RFoptOld[[2]]$internal$warn_oldstyle)
    warning("The function is obsolete. Use 'RFinterpolate' instead.")
  if (!is.null(trend)) 
    stop("in the obsolete setting, Kriging may not be used with trend!\n")                  
                    
  model <- PrepareModelOld(model, param, trend)
  if (!missing(gridtriple)) {
    if (gridtriple) {
      if (!is.null(x)) {
        if (is.matrix(x)) {
          x <- apply(x, 2,
                     function(s) {c(s[1], s[3], length(seq(s[1], s[2], s[3])))})
        } else {
          x <- seq(x[1], x[2], x[3])
          if (length(y)!=0)  y <- seq(y[1], y[2], y[3])
          if (length(z)!=0)  z <- seq(z[1], z[2], z[3])
        }
      }
    }
  }
  if (length(T)!=0)  T <- c(T[1], T[3], length(seq(T[1], T[2], T[3])))

  data <- cbind(given, data)
  colnames(data) <- c(rep("", ncol(given)), "data")
  
  RFinterpolate(x=x, y=y, z=z, T=T, grid=grid, model=model, data=data)
}


CondSimu <- function(krige.method, x, y=NULL, z=NULL, T=NULL,
                     grid=NULL, gridtriple=FALSE,
                     model, param, method=NULL,
                     given, data, trend=NULL,
                     n=1, register=0, 
                     err.model=NULL, err.param=NULL, err.method=NULL,
                     err.register=1, 
                     tol=1E-5, pch=".", #
                     paired=FALSE,
                     na.rm=FALSE) {
  RFoptOld <- internal.rfoptions(register=register, gauss.paired=paired,
                                        # gauss.method=method,
                                 general.errregister=err.register,
                                 spConform=FALSE)
  on.exit(RFoptions(LIST=RFoptOld[[1]]))
  if (RFoptOld[[2]]$internal$warn_oldstyle)
     warning("The function is obsolete. Use 'RFsimulate' instead.")
  
  if (RFoptOld[[2]]$internal$warn_oldstyle)
    warning("The function is obsolete.\nUse 'RFsimulate' instead.")
  if (!is.null(err.method)) warning("err.method is ignored.")
  
  ##  Print(krige.method, is.character(krige.method), (is.na(pmatch(krige.method, c("S","O")))))
  
 if (is.character(krige.method) && is.na(pmatch(krige.method, c("S","O"))))
   stop("Sorry. The parameters of the function `CondSimu' as been redefined. Use `krige.method' instead of `method'.")

  model <- PrepareModelOld(model, param, trend)
  if (!is.null(trend)) stop("trend cannot be given anymore.")
  err.model <- PrepareModelOld(err.model, err.param)
  if (err.model[[1]]=="+" && length(err.model) == 3 &&
      err.model[[3]][[1]]=="trend" && err.model[[3]]$mean != 0) {
    warning("error model has non-zero mean. Set to zero.")
    err.model <- err.model[[2]]
  } 
  if (!missing(gridtriple)) {
    if (gridtriple) {
      if (!is.null(x)) {
        if (is.matrix(x)) {
          x <- apply(x, 2,
                     function(s) {c(s[1], s[3], length(seq(s[1], s[2], s[3])))})
        } else {
          x <- seq(x[1], x[2], x[3])
          if (length(y)!=0)  y <- seq(y[1], y[2], y[3])
          if (length(z)!=0)  z <- seq(z[1], z[2], z[3])
        }
      }
    }
  }
  if (length(T)!=0)  T <- c(T[1], T[3], length(seq(T[1], T[2], T[3])))
  
 
  data <- cbind(given, data)
  if (is.null(dim(given)))
    given <- matrix(given)
  colnames(data) <- c(rep("", ncol(given)), "data")
  return(rfCondGauss(model, x=x, y=y, z=z, T=T, grid=grid, n=n,
                     data=data, err.model=err.model)$simu)
  
}


  
hurst <- function(x, y = NULL, z = NULL, data,
                  gridtriple = FALSE, sort=TRUE,
                  block.sequ = unique(round(exp(seq(log(min(3000, dim[1] / 5)),
                    log(dim[1]), len=min(100, dim[1]))))),
                  fft.m = c(1, min(1000, (fft.len - 1) / 10)),
                  fft.max.length = Inf, ## longer ts are cut down
                  method=c("dfa", "fft", "var"),
                  mode=c("plot", "interactive"),
                   pch=16, cex=0.2, cex.main=0.85,
                  PrintLevel=RFoptions()$basic$printlevel,
                  height=3.5,
                  ...) {  
  if (RFoptions()$internal$warn_oldstyle)
    warning("'hurst' is obsolete. Use RFhurst instead.")
  fft.len <- min(dim[1], fft.max.length)
 
  if (missing(x)) {
    gridtriple <- TRUE
    ## stopifnot(grid) 
    x <- if (is.array(data)) rbind(1, dim(data), 1) else c(1, length(data), 1)
  }
  
  if (gridtriple) {
    if (!is.null(x)) {
      if (is.matrix(x)) {
        x <- apply(x, 2,
                   function(s) {c(s[1], s[3], length(seq(s[1], s[2], s[3])))})
      } else {
        x <- seq(x[1], x[2], x[3])
        if (length(y)!=0)  y <- seq(y[1], y[2], y[3])
        if (length(z)!=0)  z <- seq(z[1], z[2], z[3])
      }
    }
  }
  
  RFhurst(x=x, y=y, z=z, data=data, sort=sort,
          block.sequ=block.sequ, fft.m=fft.m, fft.max.length=fft.max.length,
          #gauss.method=method,
          #mode=mode, pch=pch, cex=cex, cex.main=cex.main,
          printlevel=PrintLevel, height=height, ...)
}



fractal.dim <-
  function(x, y = NULL, z = NULL, data,
           grid=TRUE, gridtriple = FALSE,
           bins,
           vario.n=5,
           sort=TRUE,
           #box.sequ=unique(round(exp(seq(log(1),
           #  log(min(dim - 1, 50)), len=100)))),
           #box.enlarge.y=1,
           #range.sequ=unique(round(exp(seq(log(1),
           #  log(min(dim - 1, 50)), len=100)))),
           fft.m = c(65, 86), ## in % of range of l.lambda
           fft.max.length=Inf,
           fft.max.regr=150000,
           fft.shift = 50, # in %; 50:WOSA; 100: no overlapping
           method=c("variogram", "fft"),# "box","range", not correctly implement.
           mode=c("plot", "interactive"),
           pch=16, cex=0.2, cex.main=0.85,
           PrintLevel = RFoptions()$basic$printlevel,
           height=3.5,
           ...) {

    if (RFoptions()$internal$warn_oldstyle)
      warning("`fractal.dim' is obsolete. Use RFfractaldim instead.")
    
    if (missing(x)) {
      gridtriple <- TRUE
      stopifnot(grid)
      x <- if (is.array(data)) rbind(1, dim(data), 1) else c(1, length(data), 1)
    }
    
    if (gridtriple) {
      if (!is.null(x)) {
        if (is.matrix(x)) {
          x <- apply(x, 2,
                     function(s) {c(s[1], s[3], length(seq(s[1], s[2], s[3])))})
        } else {
          x <- seq(x[1], x[2], x[3])
          if (length(y)!=0)  y <- seq(y[1], y[2], y[3])
          if (length(z)!=0)  z <- seq(z[1], z[2], z[3])
        }
      }
    }

    RFfractaldim(x=x, y=y, z=z, data=data, grid=grid, bins=bins, vario.n=vario.n,
                 sort=sort, fft.m=fft.m, fft.max.length=fft.max.length,
                 fft.max.regr=fft.max.regr, fft.shift=fft.shift,
                # gauss.method=method,
                 mode=mode, pch=pch, cex=cex, cex.main=cex.main,
                 printlevel=PrintLevel, heigth=height, ...)
    
  }




fitvario <-
  function(x, y=NULL, z=NULL, T=NULL, data, model, param,
           lower=NULL, upper=NULL, sill=NA, grid=!missing(gridtriple),
           gridtriple=FALSE,
            ...) {
 if (RFoptions()$internal$warn_oldstyle)
    warning("The function is obsolete. Use 'RFfit' instead.")
    fitvario.default(x=x, y=y, z=z, T=T, data=data, model=model, param=param,
                     lower=lower, upper=upper, sill=sill,
                     grid=grid, gridtriple,     
                     ...)
  }

fitvario.default <-
  function(x, y=NULL, z=NULL, T=NULL, data, model, param,
           grid=!missing(gridtriple), gridtriple=FALSE,     
           trend = NULL,
##         BoxCox ((x+c)^l - 1) / l; log(x+c); with c==1
##         BC.lambda : NA / value
##         BC.c: nur value
           BC.lambda, ## if missing then no BoxCox-Trafo
           BC.lambdaLB=-10, BC.lambdaUB=10,
           lower=NULL, upper=NULL, sill=NA,
           use.naturalscaling=FALSE,
           ## speed=FALSE, i.e. coordniates will be saved in GATTER
           PrintLevel=RFoptions()$basic$printlevel, optim.control=NULL,
           bins=20, nphi=1, ntheta=1, ntime=20,
           bin.dist.factor=0.5,
           upperbound.scale.factor=3, lowerbound.scale.factor=3, 
           lowerbound.scale.LS.factor=5,
           upperbound.var.factor=10, lowerbound.var.factor=100,     
           lowerbound.sill=1E-10, scale.max.relative.factor=1000,
           minbounddistance=0.001, minboundreldist=0.02,
           approximate.functioncalls=50, 
           pch=RFoptions()$general$pch, 
           transform=NULL,  standard.style=NA,
     ##      var.name="X", time.name="T",
           lsq.methods=c("self", "plain", "sqrt.nr", "sd.inv", "internal"),
           ## "internal" : name should not be changed; should always be last
           ##              method!
           mle.methods=c("ml"), # "reml", "rml1"),
           cross.methods=NULL,
  #       cross.methods=c("cross.sq", "cross.abs", "cross.ign", "cross.crps"),
           users.guess=NULL,  only.users = FALSE,
           distances=NULL, truedim,
           solvesigma = NA, 
           allowdistanceZero = FALSE,
           na.rm = TRUE) { ## do not use FALSE for mixed models !

    if (RFoptOld[[2]]$internal$warn_oldstyle)
      warning("fitvario is obsolete. Use RFfit instead.")

     
  RFoptOld <-
    internal.rfoptions(fit.boxcox_lb=BC.lambdaLB,
            fit.boxcox_ub=BC.lambdaUB,
     #       fit.sill=sill,
            fit.use_naturalscaling= use.naturalscaling,
            printlevel=PrintLevel, 
            empvario.bins=bins,
            empvario.nphi=nphi,
            empvario.ntheta=ntheta,
            empvario.deltaT=ntime,
            fit.bin_dist_factor = bin.dist.factor,
            fit.upperbound_scale_factor = upperbound.scale.factor,
            fit.lowerbound_scale_factor = lowerbound.scale.factor,
            fit.lowerbound_scale_ls_factor=lowerbound.scale.LS.factor,
            fit.upperbound_var_factor=upperbound.var.factor,
            fit.lowerbound_var_factor=lowerbound.var.factor,     
     #       fit.lowerbound_sill=lowerbound.sill,
            fit.scale_max_relative_factor=scale.max.relative.factor,
            fit.minbounddistance=minbounddistance,
            fit.minboundreldist=minboundreldist,
            fit.approximate_functioncalls=approximate.functioncalls,
            pch=pch, 
    #        fit.optim_var_elimination=
    #                   if(is.na(standard.style)) 'try' else if (standard.style)
    #                   'yes' else 'never',
            ## var.name="X", time.name="T",
            fit.only_users = only.users,
    #        fit.solvesigma = solvesigma , 
            allowdistanceZero = allowdistanceZero,
                       na_rm = na.rm)
  on.exit(RFoptions(LIST=RFoptOld[[1]]))


   model <- PrepareModelOld(model, param, trend)
    
  if (gridtriple) {
    if (!is.null(x)) {
      if (is.matrix(x)) {
        x <- apply(x, 2,
                   function(s) {c(s[1], s[3], length(seq(s[1], s[2], s[3])))})
      } else {
        x <- seq(x[1], x[2], x[3])
        if (length(y)!=0)  y <- seq(y[1], y[2], y[3])
        if (length(z)!=0)  z <- seq(z[1], z[2], z[3])
       }
    }
  }
  if (length(T)!=0)  T <- c(T[1], T[3], length(seq(T[1], T[2], T[3])))
 

  if (!missing(param)) {
    if (is.numeric(lower)) lower <- PrepareModelOld(model=model, param=lower)
    if (is.numeric(upper)) upper <- PrepareModelOld(model=model, param=upper)
  }
   
  RFfit(x=x, y=y, z=z, T=T, data=data, model=model,
        lower=lower, upper=upper, grid=grid,boxcox=BC.lambda,
        sub.methods=lsq.methods, methods=mle.methods,
        users.guess=users.guess, distances=distances,
        dim= truedim,
        optim.control=optim.control, transform=transform,
        spConform=FALSE)
}

DeleteAllRegisters <- function() {
  old <- RFoptions()$general$storing
  RFoptions(storing=FALSE, storing=old)
}

DeleteRegister <- function(nr=0){
 DeleteAllRegisters() 
}


RFparameters <- function(...) {
  if (RFoptions()$internal$warn_oldstyle)
    warning("The function is obsolete. Use 'RFoptions' instead.")
  l <- list(...)
  nl <- names(l)
  names <- c("Storing", "PrintLevel") #
  if (any(is.na(pmatch(nl, names))))
    warning("some options of 'RFparameters' in the package 'RandomFields' are not recognized anymore. Use 'RFoptions' instead")
  if ("Storing" %in% nl) RFoptions(storing=l$Storing)
  if ("PrintLevel" %in% nl) RFoptions(printlevel=l$PrintLevel)
  
}


plotWithCircles <- function(data, factor=1.0,
                            xlim=range(data[,1])+c(-maxr,maxr),
                            ylim=range(data[,2])+c(-maxr,maxr),
                            col=1, fill=0, ...) {
  ## marked point process: presents positive values of data as radii of circles
  CIRLLE.X <- cos(seq(0,2*pi,l=20))
  CIRLLE.Y <- sin(seq(0,2*pi,l=20))
  circle <- function(x,r) { polygon(x[1]+ r* CIRLLE.X,x[2]+ r* CIRLLE.Y,
                                    col=fill, border=col) }
  ##r <- x$NormedData - min(x$NormedData) +1
  ##r <- r/max(r)/nrow(x$coords) * diff(xlim) * diff(ylim) * 2.5;
  maxr <- max(data[,3])
  plot(Inf, Inf, xlim=xlim, ylim=ylim, xlab="", ylab="",...)
  for (i in 1:nrow(data)) { circle(data[i,c(1,2)], factor*data[i,3]) }
}




RMid <- function() {
  warning("'RMid' is an obsolete function. Use 'RMidcoord' instead.")
  RMidcoord()
}
