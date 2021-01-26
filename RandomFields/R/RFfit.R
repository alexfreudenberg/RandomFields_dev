
## Authors 
## Martin Schlather, schlather@math.uni-mannheim.de
##
##
## Copyright (C) 2015 -- 2017 Martin Schlather
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



### !!!!!!!!!!!! ACHTUNG !!!!!!!!!!!! TREND als cov-fct muss
### noch programmiert werden !!!

##   source("~/R/RF/RandomFields/R/MLES.R")

## PrintLevels
## 0 : no message
## 1 : important error messages
## 2 : warnings
## 3 : minium debugging information
## 5 : extended debugging information

## jetzt nur noch global naturalscaling (ja / nein)
## spaeter eine Funktion schreibbar, die den naturscaling umwandelt;
##   im prinzipt CMbuild, aber ruechwaers mit 1/newscale und eingefuegt
##   in eventuell schon vorhandene $ operatoren


#Beim paper lesen im Zug nach Muenchen heute morgen ist mir eine Referenz zu einem R Paket "mlegp: Maximum likelihood estimates of Gaussian processes" aufgefallen. Ist Dir aber sicher schon bekannt! 

#  stop("")
  # problem: natscale; im moment 2x implementiert, 1x mal ueber
  # scale/aniso (user) und einmal gedoppelt -- irgendwas muss raus

## LSQ variogram fuer trend = const.
## kann verbessert werden, insb. fuer fixed effects, aber auch eingeschraenkt
## fuer random effects -> BA/MA


## REML fehlt

## users.guess muss in eine List von meheren Vorschlaegen umgewandelt werden !!! Und dann muss RFfit recursiver call mit allen bisherigen Werden laufen !!


## zentrale C -Schnittstellen
##    .C(C_PutValuesAtNA, RegNr, param)

## bins bei Distances automatisch


## bei repet sind die Trends/fixed effects gleich, es muessen aber die
## random effects unterschiedlich sein.
## bei list(data) werden auch trend/fixed effects unterschiedlich geschaetzt.


## Erweiterungen: Emilio's Bi-MLE, Covarianz-Matrix-INversion per fft oder
## per INLA, grosse Datensaetze spalten in kleinere "unabhaengige".


###################################
## !!! Mixed Model Equations !!! ##
###################################

## accessing slots
accessByNameOrNumber <- function(x, i, j, drop=FALSE) {
  stopifnot(length(i)==1)
  if (is.numeric(i))  i <- slotNames(x)[i]
  return(accessSlotsByName(x=x, i=i, j=j, drop=drop))
}

setMethod("[", signature = CLASS_FITLIST, def=accessByNameOrNumber)


### to do : ask Paulo
#effects_RFfit <- function(OP, object, method) {
#  eff <- RFrandef(object=object, method=method, OP=OP)
#  linpart <- fitted_RFfit(OP=OP, object=object, method=method)
#  stop("unclear how these two results should be combined in the output")
#}
#effects_RMmodelFit <- function(...) stop("'effects' can only be used with the original and sp_conform output of 'R Ffit'.")
#setMethod(f="effects", signature='R Ffit',
#          definition=function(object, method="ml")
#          effects_RFfit("@", object=object, method=method))#
#setMethod(f="effects", signature='RMmodelFit',
#          definition=function(object, newdata=NULL) effects_RMmodelFit())#
#effects.RM_mod elFit <- function(object, ...) effects_RMmodelFit()
#effects.RF_ fit <- function(object, method="ml") effects_RMmodelFit()



RFhessian <- function(model) {
  method <- "ml"
  if (is(model, CLASS_FITLIST))
    return(if (isS4(model)) model[method]$hessian else model[[method]]@hessian)
  else stop("'model' is not an output of a fitting method")
}


boundary_values <- function(variab) {
  txt <- NULL
  if (length(variab) > 0) {
    upper.bound <- variab[4, , drop=FALSE]
    lower.bound <- variab[3, , drop=FALSE]
  # sd <- variab[2, ]
    variab <- variab[1, , drop=FALSE]
    lidx <- variab < lower.bound + 1e-8
    uidx <- variab > upper.bound - 1e-8
    nl <- sum(lidx, na.rm=TRUE)
    nu <- sum(uidx, na.rm=TRUE)
    if (nl + nu > 0) {
      lidx[is.na(lidx)] <- FALSE
      uidx[is.na(uidx)] <- FALSE    
    txt <-
      paste(sep="", "Note that the (possibly internal) fitted variable",
            if (nl > 0)
              paste(if (nl > 1) "s " else " ",
                    paste("'", colnames(variab)[lidx], "'", sep="", collapse=", "),
                    if (nl == 1)  " is " else " are ",                  
                    "close to or on the effective lower boundary", sep=""),
            if (nl > 0 && nu > 0) " and the variable",
            if (nu > 0)
              paste(if (nu > 1) "s " else " ", 
                    paste("'", colnames(variab)[uidx], "'",
                        sep="", collapse=", "),
                    if (nu == 1) "is" else "are",
                    "close to or on the effective upper boundary"),
            ".\nHence the gradient of the likelihood function might not be zero and none of the\nreported 'sd' values might be reliable.")
    }
  }
  return(txt)
}



summary.RMmodelFit <- function(object, ..., isna.param) {
#  print("objecyt")
#  str(object)
  if (isS4(object)) {
    OP <- "@"
    model <- object@modelAsList
  } else {
    OP <-"$"
    model <- object$model
  }

  for (i in c("covariat", "globalvariance", "param", "residuals", "variab",
              "likelihood", "AIC", "AICc", "BIC", "coordsystem", "params.list"))
    assign(i, do.call(OP, list(object, i)))
#  covariat <- do.call(OP, list(object, "covariat"))
#  globalvariance <- do.call(OP, list(object, "globalvariance"))
#  param <- do.call(OP, list(object, "param"))
#  residuals <- do.call(OP, list(object, "residuals"))
#  variab <- do.call(OP, list(object, "variab"))
#  do.call(OP, list(object, "likelihood")),
#  AIC = do.call(OP, list(object, "AIC")),
#  AICc= do.call(OP, list(object, "AICc")),
  ## BIC = do.call(OP, list(object, "BIC")),
  
  if (length(variab) + length(covariat) == 0) {
    l <- list(model=model)
  } else {  
    l <- list(model=model, loglikelihood=likelihood,
              AIC = AIC, AICc = AICc, BIC = BIC,
              residuals=if (length(residuals) == 1) residuals[[1]]
                        else residuals)
    if (missing(isna.param)) isna.param <- any(is.na(param)) 
    l$boundary <- boundary_values(variab)
    if (length(covariat) > 0)  covariat <- as.matrix(covariat)    
    
    if (length(param) > 0 && !any(is.na(param[1, ]))) {
      nr_p <- nrow(param)
      if (length(globalvariance) > 0)
        globalvariance <-
          c(globalvariance, rep(NA, nr_p -length(globalvariance)))
      l$param <- cbind(param, globalvariance,
                       if (length(covariat) > 0)
                         rbind(covariat, matrix(NA, ncol=ncol(covariat),
                                                nrow= nr_p - nrow(covariat))))
      l$params.list <- params.list
    }

    cartesian == all(coordsystem == "cartesisan")
    if (!cartesian) l$coordsystem <- if (all(coordsystem == coordsystem[1]))
                                     coordsystem[1] else coordsystem
    if (isna.param || !is.null(l$boundary) || length(l$coordsystem) > 1) {
      nr_v <- nrow(variab)
      if (length(globalvariance) > 0)
        globalvariance <-
          c(globalvariance, rep(NA, nr_v - length(globalvariance)))
      l$variab <- cbind(variab, globalvariance,
                        if (length(covariat) > 0)
                          rbind(covariat, matrix(NA, ncol=ncol(covariat),
                                                 nrow=nr_v - nrow(covariat)))
                      )
    }
  }
  class(l) <- "summary.RMmodelFit"
  l
}

setMethod(f="summary", signature=CLASS_SINGLEFIT, summary.RMmodelFit)#


print.summary.RMmodelFit <- function(x, ...) {

  printVariab <- function(x, coordsystem) {
    cat("Internal variables")
    if (length(coordsystem) == 1) cat(" (", coordsystem, ")", sep="")
    else if (length(coordsystem) > 1) cat(" (various coordinate systems used)")
    cat(":\n")
    v <- x$variab
    if (is.null(x$boundary)) v <- x$variab[1:2, , drop=FALSE]
    if (length(coordsystem) > 1) v <- rbind(v, coordsystem)
    print(v, ..., na.print="-", quote=FALSE) #
    cat("\n")
    return(ncol(x$variab))
  }

  printParam <- function(param, coordsystem) {
    cat("Model variables")
    if (length(coordsystem) == 1) cat(" (", coordsystem, ")", sep="")
    else if (length(coordsystem) > 1) cat(" (various coordinate systems used - see 'Internal variables')")
    cat(":\n")
    print(param, ..., na.print="-")#
    return(ncol(param))
  }

  printUsers <- function(param) {
    cat("\nUser's variables:\n")
    idx <- sapply(param, length) == 1
    print(unlist(param[idx]), ..., na.print="-")#
    if (!all(idx)) {
      p <- param[idx]
      n <- names(p)
      for (i in 1:length(p)){
        cat("user variable '", n[i], "':\n", sep="")
        print(p[[i]], ..., na.print="-")#
      }
    }
    return(ncol(param))
  }
  
  printRest <- function(...) {
    x <- unlist(list(...))
    stopifnot(length(x) == 3)
    names(x) <- c("#variab", "loglikelihood", "AIC")
    cat("\n")
    print(x) #
    cat("\n")
  }

  if (internalRFoptions(getoptions_="general")$detailed_output)
    str(x$model, no.list=TRUE) #
  cat("\n")
  np <- AIC <- ll <- nm <- NA
  cparam <- NULL
 
  if (length(x$param) > 0) {
    if (length(x$submodels) > 0) {
      cur_name <- ""
      len <- length(x$submodels)
      
      for (i in 1:len) {
        sm <- x$submodels[[i]]      
        n <- sm$report
        nnxt <- if (i==len) "" else x$submodels[[i+1]]      
        if (n != cur_name) {
          if (i > 1) {         
            if (!is.null(sm$param)) printParam(cparam, x$coordsystem)
             printRest(np, ll, AIC) #
            if (!is.null(sm$boundary)) cat(sm$boundary, "\n\n")
          }
          
          if (nnxt != n && length(sm$fixed) > 0) {
            
            nX <- paste(sep="", n, " (",
                        paste(c(if (length(sm$fixed$zero) > 0)
                                  paste(colnames(x$param)[sm$fixed$zero],"= 0"),
                                if (length(sm$fixed$one) > 0)
                                  paste(colnames(x$param)[sm$fixed$one],"= 1")),
                              sep=", "),
                        ")")
          } else nX <- n
          cat(if (!is.na(nm)) cat("\n"), nX, "\n",
              paste(rep("=", min(80, nchar(nX))), collapse=""),
              "\n", sep="")
          np <- 0
          AIC <- 0
          ll <- 0
          cparam <- NULL
          nm <- 1
        }
        if (!is.null(sm$variab)) {
          if (nm > 1 || (i<len && n==nnxt)) cat("model", nm, ", ")
          printVariab(sm, x$coordsystem)
        }
        if (!is.null(sm$param)) {
          param <- x$param * NA
          param[, sm$p.proj] <- sm$param
          fixed <- sm$fixed
          if (length(fixed) > 0) {
            param[1, fixed$zero] <- 0
            param[1, fixed$one] <- 1
          }
          
          ##  if (!is.null(cparam)) cparam <- rbind(cparam, NA)
          cparam <- rbind(cparam, param)
          
        }
        np <- np + length(sm$p.proj)
        ll <- ll + sm$loglikelihood
        AIC <- AIC + sm$AIC
        nm <- nm + 1;
        cur_name <- n
      }
      
      if (!is.null(sm$param)) printParam(param, x$coordsystem)
      printRest(np, ll, AIC) #
      if (!is.null(sm$boundary)) cat(sm$boundary, "\n\n")
      
      cat("\nuser's model\n", paste(rep("=", 12), collapse=""), "\n", sep="")
    }

    
    np <- NA
    np <- printVariab(x, x$coordsystem)
    if (length(x$param) > 0) np <- printParam(x$param, x$coordsystem)
    if (length(x$params.list) > 0) printUsers(x$params.list)
    
    printRest(np, x[c("loglikelihood", "AIC")])#
    if (!is.null(x$boundary)) cat(x$boundary, "\n\n")
    class(x) <- NULL
##    Print(x)
##    Print(x$params.list)
  } else {
    cat("The model does not have a parameter to fit.\n")
  }
   
  invisible(x)
}

print.RMmodelFit <- function(x, ...)
  print.summary.RMmodelFit(summary.RMmodelFit(x, ...))#

setMethod(f="show", signature=CLASS_SINGLEFIT,
          definition=function(object) print.RMmodelFit(object))#

getRFfitmethod <- function(object, method, all=FALSE) {
  IsS4 <- isS4(object)
  if (missing(method) || length(method) == 0) method <- METHOD_PREFLIST
  ok <- logical(length(method))
  for (m in 1:length(method)) {
#    Print(object)
    ob <- if (IsS4) object[method[m]] else object[[method[m]]]
    if (is.null(ob) && !all)
      ## if not !IsS4 and user has given non-existing method
       warning( paste("The following method does not exist: ", method[m]))
#    Print(ob)
    if ( (ok[m] <- length(if (isS4(ob)) ob@likelihood else ob$likelihood) > 0)
        && !all) return(method[m])
##     ok[m]  <- if (length(do.call(OP, list(newx$ml, "name")))>0)
  }
  if (all) return(method[ok])
  stop("no specified model found")
#} else stop("unknown class in getrffitmethod")
}
  
summary.RFfit <- function(object, ...,  method, full=FALSE) {
  IsS4 <- isS4(object)
  m <- getRFfitmethod(object, method)
  s <- summary.RMmodelFit(if (IsS4) object[m] else object[[m]])
  osm <- if (IsS4) object@submodels else object$submodels
  len <-  length(osm)
  if (full && len > 0) {
    submodels <- list()
    for (i in 1:len) {
      submodels[[i]] <- summary(osm[[i]][[method]],# 'summary'
                                isna.param=is.null(s$param))    # nicht
      submodels[[i]]$report <- osm[[i]]$report  # spezifizieren!
      submodels[[i]]$p.proj <- osm[[i]]$p.proj
      submodels[[i]]$fixed  <- osm[[i]]$fixed
    }     
    s$submodels <- submodels
  }
  s
}



print.RFfit <- function(x, ...,  method, full=FALSE) {
  print.summary.RMmodelFit(summary.RFfit(x, ..., method=method, full=full))
}

setMethod(f="show", signature=CLASS_FITLIST,
          definition=function(object) print.RFfit(object))#


print.AICRFfit<- function(x, ..., digits=3) {
  ## nur deshalb 
  fstcol <- 3
  sndcol <- 55
  trdcol <- 4
  forthcol<-9
  leer <- formatC("", width=fstcol)
  size <- max(abs(x[[2]]))
  size <- if (size>0) ceiling(log(size) / log(10)) else 1
  cat(leer, formatC("model", flag="-", width=sndcol), " ",
      formatC(names(x)[1], width=trdcol),
      formatC(names(x)[2], width=forthcol), "\n", sep="")
  names <- attr(x, "row.names")
  for (i in 1:length(names)) {
    cat(formatC(i, width=fstcol, flag="-"))
    if (nchar(xx <- names[i]) <= sndcol)
      cat(formatC(xx, width=sndcol, flag="-"))
    else {
      yy <- strsplit(xx, " \\* ")[[1]]
      for (j in 1:length(yy)) {
        ncyy <- nchar(yy[j])
        if (ncyy <= sndcol && j==length(yy))
          cat(format(yy[j], width=sndcol, flag="-"))
        else {
          if (ncyy <= sndcol - 2) {
            cat(yy[j])
          } else {
            zz <- strsplit(yy[j], ", ")[[1]]
            ncyy <-  0
            lenzz <- length(zz)
            for (k in 1:lenzz) {
              len <- nchar(zz[k])
              if (k > 1 && len > sndcol - 1) {
                cat("\n", leer, zz[k], sep="")
                if (k < lenzz)
                  cat(formatC(",", flag="-",  width=pmax(1, sndcol-len)))
              } else {
                if (ncyy + len > sndcol - 1) {
                  cat("\n", leer, sep="")
                  ncyy <- len
                } else {
                  ncyy <- ncyy + len
                }
                cat(zz[k])
                if (k < lenzz) {
                  cat(", ")
                  ncyy <- ncyy + 2
                }
              }
            } # for k 1:lenzz
          } # split according to commata
          if (j < length(yy)) cat(" *\n", leer, sep="")
          else if (ncyy < sndcol) cat(formatC("", width=sndcol-ncyy))
        }
      } # for 1:products
    } ## not be written in a single line
    cat("",
        formatC(x[[1]][i], width=trdcol),
        formatC(x[[2]][i], format="f", width=size + digits + 1,
                    digits=digits),"\n")
  }
}



fullAIC <- function(x, method="ml", AIC="AIC") {
  ats <- approx_test_single(x, method=method)$result
  values <- c("name", "df", AIC)
  model2 <- paste("model2.", values, sep="")
  ats2 <- ats[ !is.na(ats[, model2[2]]), model2]
  colnames(ats2) <- values
  if (ats2$df < 0) ats2 <- NULL
  ats <- ats[, paste("model1.", values, sep="")]
  colnames(ats) <- values
  if (ats$df < 0) ats <- NULL
  ats <- unique(rbind(ats, ats2))
  dimnames(ats) <- list(1:nrow(ats), colnames(ats))

  names  <- as.character(ats$name)
  ats <- ats[-1]
  attr(ats, "row.names") <- names  
  class(ats) <- "AICRFfit"
  ats
}

AIC.RFfit <- function(object, ..., k=2, method="ml", full=TRUE) {
  if (full) {
    fullAIC(object, method=method)
  } else {
    AIC <- if (isS4(object)) object[method]@AIC else object[[method]]$AIC
    names(AIC) <- "AIC"
    AIC
  }
}

AICc.RFfit <- function(object, ...,  method="ml", full=FALSE) {
  if (full) {
    stop("for 'AICc' the option 'full=TRUE' has not been programmed yet.")
    fullAIC(object, method=method)
  } else {
    AIC <- if (isS4(object)) object[method]@AICc else object[[method]]$AICc
    names(AIC) <- "AICc"
    AIC
  }
}

BIC.RFfit <- function(object, ..., method="ml", full=TRUE) {
  if (full) {
    fullAIC(object, method=method, AIC="BIC")
  } else {
    BIC <- if (isS4(object)) object[method]@BIC else object[[method]]$BIC
    names(BIC) <- "BIC"
    BIC
  }
}

setMethod(f="plot", signature(x=CLASS_FITLIST, y="missing"),
          function(x, y, ...) RFplotEmpVariogram(x, ...))
setMethod(f="plot", signature(x=CLASS_SINGLEFIT, y="missing"),
          function(x, y, ...) RFplotModel(x, ...))
setMethod(f="persp", signature(x=CLASS_FITLIST),
	  function(x, ...) RFplotEmpVariogram(x, ..., plotmethod="persp"))


contour.RFfit <- contour.RFempVariog <- 
  function(x,...) {
    stopifnot(!( (is(x, CLASS_FITLIST) && is.list(x@ev@centers))
                || (is(x, CLASS_EMPIR) && is.list(x@centers))
                ))
    RFplotEmpVariogram(x, ..., plotmethod="contour")
  }


ExpliciteGauss <- function(model) {
  if (model[[1]] != "RPgauss" && model[[1]] != "gauss.process") {
    boxcox <- getRFoptions("gauss")$boxcox
    if (any(is.na(boxcox)) || any(boxcox[c(TRUE, FALSE)] != Inf))
      return(list("RPgauss", boxcox=boxcox, model))
  }
  return(model)
}




RFfit <- function(model, x, y=NULL, z=NULL, T=NULL,  grid=NULL, data, 
           lower=NULL, upper=NULL, 
           methods, # "reml", "rml1"),
           sub.methods,
           ## "internal" : name should not be changed; should always be last
           ##              method!
           optim.control=NULL,
           users.guess=NULL,  
           distances=NULL, dim,
           params=NULL,
            ##type = c("Gauss", "BrownResnick", "Smith", "Schlather",
           ##             "Poisson"),
           ...)
{

  .C(C_NoCurrentRegister)

  RFopt <-internalRFoptions(xyz=length(y)!=0,..., fit.addNAlintrend = 2,
                            internal.examples_reduced = FALSE,
                            messages.warn_singlevariab=FALSE)
  if (!hasArg("COPY")) on.exit(optionsDelete(RFopt))
  
  if (length(params) > 0) {
    if ((!is.na(RFopt$fit$estimate_variance_globally) &&
         RFopt$fit$estimate_variance_globally) &&
        RFopt$basic$printlevel > 0)
      message("Value of option 'estimate_variance_globally' is ignored.")
    RFopt$fit$estimate_variance_globally <- FALSE
    setRFoptions(fit.estimate_variance_globally = FALSE)
  }
  fit <- RFopt$fit

  if (RFopt$general$vdim_close_together)
    stop("'vdim_close_together' must be FALSE")

  ## in U nifyData the further.models that contain only the parameter data
  ## are turned into genuine models 
  further.models <- list()
  models <- c("lower", "upper", "users.guess", "parscale")
  parscale <- optim.control$parscale ## could be given by RMmodel!!
  if (paramlist <- length(params) > 0)
    for (m in models) further.models[[m]] <- get(m)

##  Print(model, if(!missing(x)) x)

  Z <- UnifyData(model=model, x=x, y=y, z=z, T=T, grid=grid,
                 data=data, distances=distances, dim=dim,
                 RFopt=RFopt,
 		 further.models = further.models,
                 model.dependent.further.models = paramlist,
                 params=params, return_transform=TRUE,
                 ...)

#  Print(Z)
 
  Z[c("rangex", "mindist", "maxdist")] <-
    GetRanges(Z, RFopt$fit$smalldataset / 2)
  
  for (m in models)  ## 6.1.20 somewhere ReplaceC??
    if (!is.null(get(m)) && !is.numeric(get(m))) assign(m,Z$further.models[[m]])
  optim.control$parscale <- parscale ## transformed parscale

  new.model <- Z$model
  if (new.model[[1]] %in% c("RPpoisson", "poisson")) {
    res <- fit.poisson()
    type <- "poisson"
  } else if (new.model[[1]] %in% c("BRmixed", "BRshifted", "BRmixedIntern",
                               "RFbrownresnick")) {
    res <- fit.br()
    type <- "brownresnick"
 } else if (new.model[[1]] %in% c("RPschlather", "extremalgauss")) {
    res <- fit.extremal.gauss()
    type <- "schlather"
  } else if (new.model[[1]] %in% c("RPsmith", "smith")) {
    res <- fit.smith()
    type <- "smith"
  } else if (new.model[[1]] %in% c("RPbernoulli", "binaryprocess")) {
    res <- fit.bernoulli()    
    type <- "bernoulli"
  } else {   
    Z$model <- ExpliciteGauss(ReplaceC(Z$model))
    if (!is.null(Z$transform) && fit$estimate_variance_globally &&
        !hasArg("fit.estimate_variance_globally") &&
        !hasArg("estimate_variance_globally")) {
      stop("Since a transform of the parameters to be estimated is given, the parameter 'estimate_variance_globally=TRUE' might be faulty. On the user's full responsibility, the user might set 'estimate_variance_globally=TRUE' explicitely within the current function call.")
    }
    res <- do.call("rffit.gauss",
                   c(list(Z, lower=lower, upper=upper, users.guess=users.guess,
                          optim.control=optim.control,
                          recall = FALSE),
                     if (!missing(methods))  list(mle.methods = methods),
                     if (!missing(sub.methods)) list(lsq.methods=sub.methods)
                     ## "internal" : name should not be changed; should always
                     ## be last method!
                     ))
    type <- "gauss"
  }
  if (is.null(res)) return(res)
  if (RFopt$general$returncall)
    attr(res, "call") <- as.character(deparse(match.call()))

  attr(res, "type") <-type
    
  ## noch notwendig??? 
  attr(res, "coord_system") <- .Call(C_GetCoordSystem,
                                     as.integer(MODEL_MLE),
                                     RFopt$coords$coord_system,
                                     RFopt$coords$new_coord_system)
   return(res)
}
