
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

# RFsimulate:  Not implemented yet: If \code{model} is a formula or of class
#    \command{\dQuote{\link{RFformula}}},
#    the corresponding linear mixed model of the type
 #   \deqn{response = W*b + Z*u + e} is simulated

##   source("~/R/RF/RandomFields/R/MLES.R")

## Printlevels
## 0 : no message
## 1 : important error messages
## 2 : warnings
## 3 : minium debugging information
## 5 : extended debugging information 

## jetzt nur noch global naturalscaling (ja / nen)
## spaeter eine unktion schreibbar, die den naturscaling umwandelt;
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


## bins bei Distances automatisch


## bei repet sind die Trends/fixed effects gleich, es muessen aber die
## random effects unterschiedlich sein.
## bei list(data) werden auch trend/fixed effects unterschiedlich geschaetzt.


## Erweiterungen: Emilio's Bi-MLE, Covarianz-Matrix-INversion per fft oder
## per INLA, grosse Datensaetze spalten in kleinere "unabhaengige".


######################################################################
## kleine Helfer:
######################################################################

TRY <- function(...) {
  a <- try(...) 
##  print(a)
  a
}
##
TRY <- Try
#TRY <- function(A, ...) A

model2string <- function(model) {
##  options(warn=0)
  if (model[[1]] == RM_DECLARE) return("")
  ans <- paste(model[[1]], "(", sep="")
  if (length(model) > 1) {
    n <- names(model)
    j <- 0
    for (i in 2:length(model)) {
      sub <- if (isModel(model[[i]])) model2string(model[[i]]) else model[[i]]
      if (all(sub=="")) next ## === any
      j <- j + 1
      ans <- paste(ans, if (j>1) ",", n[i], "=", sub)
    }   
  }
  return(paste(ans, ")", sep=""))
}
   

######################################################################
##  Transform S3 fit result into S4 
#########
list2RFempVariog <- function(li) {
  return()
}



#############################################################

list2RMmodelFit <- function(x, RFsp.info=NULL, T) {
  ## final transformations ....

 # Print(RFsp.info)
#  str(x)
#  Print("T")
#  if (!missing(T)) str(T)
#  Print("start")
  
  stopifnot(is.list(x),
            all(c("model", "likelihood", "residuals") %in% names(x)))

  
  if (!is.null(RFsp.info)) {
    ## convert residuals to RFsp class
    err <- TRY({
        lres <- length(x$residuals)
        if (lres > 0) {
          for (i in 1:lres) {
            dta <- x$residuals[[i]]
            if (length(RFsp.info) < i)  {
              gT <- co <- NULL
              vdim <- n <- 1
            } else {
              gT <- RFsp.info[[i]]$gridTopology
              co <- RFsp.info[[i]]$coords
              n <- RFsp.info[[i]]$data.params$n
              vdim <- RFsp.info[[i]]$data.params$vdim
              if (!is.null(dta)) dim(dta) <- RFsp.info[[i]]$dimensions
#              Print(RFsp.info[[i]]$dimensions)              
            }
            if (!is.null(dta)) {
              x$residuals[[i]]<-
                conventional2RFspDataFrame(data=dta, coords=co,
                                           gridTopology=gT, n=n, vdim=vdim,
                                           T = T, vdim_close_together=FALSE)
            }
          }
        }      
    })
    
   if(is(err, CLASS_TRYERROR))
      warning(paste("residuals could not be coerced to class 'RFsp';",
                    err))
  } # !is.null(RFsp.info)

  if (!is.null(x$trend)) stop("'trend' as list currently does not work anymore")

  x$model <- list2RMmodel(x$model)
  return(do.call("new", c(list(CLASS_SINGLEFIT), x)))
}
  

addRMmatrix <- function(model, m) {
  if (length(model) == 0) return(model)
  if (model[[1]] %in% RM_PLUS)
    return (lapply(model, function(x) {
      if (!is.list(x)) x
      else if (x[[1]] %in% RM_PLUS) { addRMmatrix(x, m) }
      else if (x[[1]] %in% RM_MULT) {
        return( list(RM_MATRIX, M=m, x))        ## to do
        if (any(idx <- sapply(x, function(z) z[[1]] == R_C))) {
          idx <- which(idx)
          x[[idx[1]]] <- list(RM_MATRIX, M=m, x[[idx[1]]])                
          x
        }
        else list(RM_MATRIX, M=m, x)
      } else list(RM_MATRIX, M=m, x)
    }))
  else return(list(RM_MATRIX, M=m, model))
}
                                        

#RFlsq <- function() {}


######################################################################
##    'optimal' parscale for 'optim'
######################################################################
ParScale <- function(optim.control, current, lower, upper) {
  if (!is.null(parscale <- optim.control$parscale)) {
    if (!is.numeric(parscale)) {
      stop("non numeric parscale not allowed yet")
    } # else is.numeric: nothing to do
  } else parscale <- rep(NA, length(lower))

  if (!missing(current)) {
    idx <- is.finite(parscale)
    parscale[!idx] <- current[!idx]
    parscale[idx] <- sqrt(parscale[idx] * current[idx]) # geom mittel
  }
  
  if (any(idx <- !is.finite(parscale) | parscale == 0)) {
    parscale[idx] <- pmax(abs(lower[idx]), abs(upper[idx])) / 10 # siehe fit_scale_ratio unten
    idx <- idx & (lower * upper > 0)
    parscale[idx] <- sqrt(lower[idx] * upper[idx])
    stopifnot(all(is.finite(parscale)))
  }
  return(parscale)
}


GetValuesAtNA <- function(NAmodel, valuemodel, skipchecks, type="",
                          C_coords) {
  aux.reg <- MODEL_AUX
  
  info <- neu <- list()
  models <- list(NAmodel, ReplaceC(#P repareModel2( ## 13.7.19 geaendert
                              valuemodel)) #)
  
  for (m in 1:2) {
    info[[m]] <- .Call(C_SetAndGetModelLikelihood, aux.reg,
                       list("Covariance",ExpliciteGauss(models[[m]])),
                       C_coords, noConcerns,#raw ok, aber nicht gebraucht
                       original_model)
      neu[[m]] <- GetModel(register=aux.reg, modus=GETMODEL_DEL_MLE,
                           return.which.param=INCLUDENOTRETURN)
  }
  NAs <- sum(info[[1]]$NAs) - info[[1]]$NaNs
  ret <- TRY(.Call(C_Take2ndAtNaOf1st, aux.reg, list("Covariance", neu[[1]]),
                   list("Covariance", neu[[2]]), C_coords, NAs,
                   as.logical(skipchecks), noConcerns))
  if (!is.numeric(ret)) {
    model<-neu[[1]]
    "lower/upper/users.guess" <- neu[[2]]
    stop("'", type, "' does not match 'model'.\n  Better use the formula notation for the model given in ?RFformulaAdvanced.")
  }
  
  return(ret)
}



######################################################################
##   used for detecting splitting directions
######################################################################
vary.variables <- function(variab, lower, upper) {
  n.var <- length(lower)
  w.value <- runif(n.var, 3-0.5, 3+0.5) # weight
  n.modivar <- 5

   
  idx <- lower <= 0
  modilow <- modiupp <- numeric(n.var)
  
  modilow[idx] <- (lower[idx] + w.value[idx] * variab[idx]) / (w.value[idx] + 1)
  modilow[!idx] <-(lower[!idx]*variab[!idx]^w.value[!idx])^(1/(w.value[!idx]+1))
  modiupp[idx] <- (upper[idx] + w.value[idx] * variab[idx]) / (w.value[idx] + 1)
  modiupp[!idx] <-(upper[!idx]*variab[!idx]^w.value[!idx])^(1/(w.value[!idx]+1))
  
  modivar <- matrix(nrow=n.var, ncol=n.modivar)
  for (i in 1:n.var) {
    modivar[i,] <- seq(modilow[i], modiupp[i], length.out=n.modivar)
  }
  return(modivar)
}


C_vary.x <- function(rangex) { ## TO DO: verbessern?!
  ## x are in C-coding !
  ## letzter Eintrag ist 0
  npt <- 5
  totdim <- ncol(rangex)
  x <- matrix(0,ncol=npt+1, nrow=totdim)
  for (i in 1:totdim) {
    x[i, 1:npt] <- runif(npt, rangex[1, i], rangex[2, i])
  }
  x
}



ModelSplitXT <- function(splitReg, info.cov, trafo, variab,
                         lower, upper, rangex, modelinfo, model,
                         p.proj=1:length(variab), v.proj=1,
                         report.base) {
  
  vdim <- modelinfo$vdim
  tsdim <- modelinfo$tsdim
  ts.xdim <- modelinfo$ts.xdim
  xdimOZ <- modelinfo$xdimOZ
  dist.given <- modelinfo$dist.given

  if (ts.xdim == 1) return(list(p.proj=p.proj,
        x.proj = TRUE,
        v.proj=v.proj))
  lp <- length(variab[p.proj]) # not length(p.proj) sicherheitshalber
  statiso <- info.cov$trans.inv && info.cov$isotropic
   
  abs.tol <- 0

  truely.dep <- dep <- matrix(FALSE, nrow=lp, ncol=ts.xdim)
  varyx <- C_vary.x(rangex=rangex)
  varyy <- C_vary.x(rangex=rangex)

  modivar <- vary.variables(variab=variab, lower=lower, upper=upper)

  for (d in 1:ts.xdim) {
    x <- varyx
    x[-d, ] <- 0
    if (dist.given) y <- NULL
    else {
      y <- varyy
      y[-d, ] <- 0
    }

    S <- double(ncol(x) * vdim^2)
    for (i in 1:lp) {

      for (m0 in 1:ncol(modivar)) {
        var <- modivar[, m0]
      
        for (m1 in 1:ncol(modivar)) {
	  
	  ##         var[p.proj[i]] <- modivar[i, m1] ## olga: correction
	  var[p.proj[i]] <- modivar[p.proj[i], m1]	  

          .C(C_PutValuesAtNA, as.integer(splitReg), as.double(trafo(var))) 
          .Call(C_CovLocNonGrid, splitReg, x, y, S)
          
          base::dim(S) <- c(ncol(x), vdim, vdim)
         
          if (any(is.na(S)))
            stop("Model too complex to split it. Please set split=FALSE")

          ## Bei x_j=0, j!=i, gibt es Auswirkungen hinsichtlich
          ## des Parameters auf die Kov-Fkt?
          ## d.h. aendern sich die kov-Werte bei variablem Parameter,
          ## obwohl die x_j = 0 sind?
          if (m1==1) Sstore <- S[, v.proj, v.proj, drop=FALSE]
          else {
            difference <- S[,v.proj, v.proj, drop=FALSE] - Sstore
            if (any(abs(difference) > abs.tol)) {
              dep[i, d] <- TRUE

              ## ist die Kovarianzfunktion der Bauart C_1(x_1) + v C_2(x_2),
              ## und d=1, so bildet v eine Konstante bzg. C_1(x_1)        

              ## i.e.,
              ## note: s_1 C_1(x) + s_2 C_2(t)
              ## then s_2 dep(ends) on x, but not truely
              truely.dep[i, d] <- truely.dep[i, d] || 
                 !(all(apply(difference, 2:3, function(x) all(diff(x) == 0))))
            }
          }
        }
      }
    }

  }

  modelsplit <- NULL
  untreated_x <- rep(TRUE, ts.xdim)
  untrtd_p <- rep(TRUE, lp)
  
  for (d in 1:ts.xdim) {
    if (!untreated_x[d]) next
    rd <- dep[, d]
    if (!any(rd)) next
    same.class <- which(apply(dep == rd, 2, all))
    untreated_x[same.class] <- FALSE

    ## untreated_x wird false falls rd und der Rest keine Abhaengigkeiten
    ## mehr gemein haben
    untreated_x[d] <- !all(!rd | !dep[, -same.class, drop=FALSE])
 
    ## falls es Parameter gibt, die nur zur Koordinate d gehoeren,
    ## werden diese geschaetzt, zusammen mit allen anderen Paremetern,
    ## die zur Koordinate d (und zu anderen Koordinaten) gehoeren 
    trtd <- apply(rd & !dep[, -same.class, drop=FALSE], 1, all)
    if (any(trtd)) {
      untrtd_p <- untrtd_p & !trtd
      report <- "separable"
      modelsplit[[length(modelsplit) + 1]] <-
        list(p.proj = p.proj[dep[, d]], 
             v.proj = v.proj,
             x.proj=as.integer(same.class),
             report = (if (missing(report.base)) report else
		       paste(report.base, report, sep=" : "))
	     )
    }
  }

  # !truely.dep, die moegicherweise rausgeflogen waren
  nontruely <- apply(!truely.dep & dep, 1, any)
  if (any(nontruely)) {
    
    if (sum(nontruely) > 1) stop("more than 1 parameter is detected that is not truly dependent -- please use split=FALSE")
      ## Problem ist hier die nicht Eindeutigkeit bei
    ## C = v_1 C_1(x_1) + v_2 C_2(X_2) + v_2 C_3(x_3)
    ## da nun v_2 und v_3 als nontruely auftauchen und somit
    ## dass auf x_1 projezierte Modell nicht mehr identifizierbar ist.
    ## Letztendlich muessen die beiden (oder mehrere) zusammengefasst werden
    ## und in recursive.estimation, use.new.bounds adaequat beruecksichtigt
    ## werden 
    l <- length(modelsplit)
    modelsplit[[l + 1]] <- list(use.new.bounds=p.proj[nontruely])
    modelsplit[[l + 2]] <-
      list(p.proj = p.proj[nontruely], v.proj = v.proj,
           x.proj=which(apply(nontruely & dep, 2, any))
           ) ## oder vielleicht doch TRUE oder durchzaehlen fuer
    ##                     verschiedene Ebenen
    untrtd_p <- untrtd_p & !nontruely
  }

  ## Parameter die verwoben sind:
  if (any(untrtd_p & !nontruely)) {
    l <- length(modelsplit)
    report <- "simple space-time"
    modelsplit[[l + 1]] <- list(use.new.bounds = which(untrtd_p))
    modelsplit[[l + 2]] <-
      list(p.proj=which(untrtd_p), v.proj=v.proj,
           x.proj=which(apply(untrtd_p & dep, 2, any)),
           report = (if (missing(report.base)) report else
                     paste(report.base, report, sep=" : ")))
  }

  ## komplett, falls notwendig
  if ((any(nontruely) || any(untrtd_p)) && !all(untrtd_p) && !all(nontruely)) {
    l <- length(modelsplit)
    modelsplit[[l + 1]] <- list(use.new.bounds=p.proj)
    modelsplit[[l + 2]] <- list(p.proj=p.proj, x.proj=1:ts.xdim, v.proj=v.proj)
  }

  if (length(modelsplit) == 1) {
    modelsplit <- modelsplit[[1]]
    modelsplit$report <- NULL
  }

  return(modelsplit)
}


ModelSplit <- function(splitReg, info.cov, trafo, variab,
                       lower, upper, rangex, modelinfo, model) {
  vdim <- modelinfo$vdim
  tsdim <- modelinfo$tsdim
  ts.xdim <- modelinfo$ts.xdim
  dist.given <- modelinfo$dist.given
  xdimOZ <- modelinfo$xdimOZ
  refined <- modelinfo$refined

  abs.tol <- 0
  restrictive <- FALSE

  lp <- length(variab)
  statiso <- info.cov$trans.inv && info.cov$isotropic
 # if (xor(is.dist, statiso)) stop("mismatch in ModelSplit.", CONTACT)

  if (vdim == 1) {
    if (statiso) return(NULL)
    modelsplit <- ModelSplitXT(splitReg=splitReg, info.cov=info.cov,
                               trafo=trafo, variab=variab, lower=lower,
                               upper=upper, rangex=rangex,
                               modelinfo=modelinfo, model=model)

    if (is.list(modelsplit[[1]])) {
      l <- length(modelsplit)
      modelsplit[[l+1]] <- list(use.new.bounds=1:lp)
      modelsplit[[l+2]] <-
        list(p.proj=1:lp,
             x.proj=if (ts.xdim==1 && tsdim > 1) TRUE else 1:ts.xdim,
             v.proj=1:vdim)
      return(modelsplit)
    } else return(NULL)
  }
  
  
  overlap <- dep <- matrix(FALSE, nrow=lp, ncol=vdim)

  x <- C_vary.x(rangex=rangex)
  S <- double(ncol(x) * vdim^2)
  if (!dist.given) y <- C_vary.x(rangex=rangex)

  modivar <- vary.variables(variab=variab, lower=lower, upper=upper)

  for (i in 1:lp) {

    for (m0 in 1:ncol(modivar)) {
      var <- modivar[, m0]
      
      for (m1 in 1:ncol(modivar)) {
        var[i] <- modivar[i, m1]
       .C(C_PutValuesAtNA, splitReg, trafo(var))
       .Call(C_CovLocNonGrid, splitReg, x, if (dist.given) NULL else y, S) 
        base::dim(S) <- c(ncol(x), vdim, vdim)
        if (any(is.na(S)))
          stop("model too complex to split it. Please set split=FALSE")

         if (m1==1) Sstore <- S
        else {
          for (n in 1:vdim)
            dep[i, n] <-
              dep[i, n] | any(abs(S[,n,n] - Sstore[,n,n]) > abs.tol)
        }
      }
    }
  }

  modelsplit <- list()
  untreated_v <- logical(vdim)
  report <- "indep. multivariate"
  
  for (v in 1:vdim) {
    d1 <- dep[, v]
    if (!any(d1)) next
    overlap[, v] <- apply(dep[, -v, drop=FALSE] & d1, 1, any)
    untreated_v[v] <- if (restrictive) !any(overlap[, v])
             else (sum(overlap[, v]) < sum(d1))
    if (untreated_v[v]) {
      idx <- length(modelsplit) +1
       modelsplit[[idx]] <-
        ModelSplitXT(splitReg=splitReg, info.cov=info.cov, trafo=trafo, variab=variab,
                     lower=lower, upper=upper, rangex=rangex,
                     modelinfo=modelinfo, model=model,
                     p.proj=which(d1), v.proj = v, report.base=report)
      modelsplit[[idx]]$report <- report
    }
  }
 
  ## bislang oben nur auf univariat heruntergebrochen
  ## man koennte noch auf bivariate runterbrechen
  ## dies aber erst irgendwann;
  ## wuerde ich auch zu weiteren report-Ebenenen 2,3, etc. fuehren
  
  if (any(untreated_v)) {    
    overlapping <- apply(overlap[, untreated_v, drop=FALSE], 1, any)
    depending <- apply(dep[, untreated_v, drop=FALSE], 1, any)

    if (refined && any(depending) && any(overlapping)) {
      idx <- length(modelsplit)
      depend <- which(depending | overlapping)
      notok <- which(!depending & !overlapping)
      modelsplit[[idx + 1]] <- list(use.new.bounds=depend, fix=notok)
      report <- "simple multivariate"
      modelsplit[[idx + 2]] <-
        ModelSplitXT(splitReg=splitReg, info.cov=info.cov, trafo=trafo,
                     variab=variab,
                     lower=lower, upper=upper, rangex=rangex,
                     modelinfo=modelinfo, model=model,
                     p.proj = depend, v.proj = 1:vdim, report.base=report)
      modelsplit[[idx+2]]$report <- report
    } else notok <- which(!depending | overlapping)
 
    if (length(notok) > 0) {
      idx <- length(modelsplit) + 1
      modelsplit[[idx]] <-
        ModelSplitXT(splitReg=splitReg, info.cov=info.cov, trafo=trafo, variab=variab,
                     lower=lower, upper=upper, rangex=rangex,
                     modelinfo=modelinfo, model=model,
                     p.proj = notok, v.proj = 1:vdim)
    }
  }

  l <- length(modelsplit)
  modelsplit[[l+1]] <- list(use.new.bounds=1:lp)
  modelsplit[[l+2]] <-
    list(p.proj=1:lp,
         x.proj=if (ts.xdim==1 && tsdim > 1) TRUE else 1:ts.xdim,
         v.proj=1:vdim)

  modelsplit[[1]]$all.p <- dep

  return(modelsplit) ## return(NULL) # 11.9.12
}



recurs.estim <- function(split, level, splitReg, Z,
                         lower, upper,  guess,  
                         lsq.methods, mle.methods,
                         optim.control,
                         trafo, NaNs,
                         spConform,
                         practicalrange,
                         printlevel,
                         fit,
                         minmax,
                         sdvar,                         
			 origin,
                         IDX) {
##  Print(guess, upper, lower)
  M <- if (length(mle.methods) >= 1) mle.methods[length(mle.methods)]
  else lsq.methods[length(lsq.methods)]

  w <- 0.5
   sets <- length(Z$data)

  if (printlevel >= PL_FCTN_DETAILS) Print(split) #
  submodels <- NULL
  submodels_n <- 0
  fixed <- FALSE

  for (s in 1:length(split)) {
    
    sp <- split[[s]]

    p <- sp$use.new.bounds
    if (!is.null(p)) {
      if (printlevel >= PL_STRUCTURE) {
        cat("    calculating new lower and upper bounds for the parameters ",
            paste(p, collapse=", "), "\n", seq="")
      }
       
      if (!is.null(sp$fix)) {
        fix.zero <- (minmax[sp$fix, MINMAX_PMIN] <= 0) &
          (minmax[sp$fix, MINMAX_PMAX] >=0)
        fix.one <- (minmax[sp$fix, MINMAX_PMIN] <= 1) &
          (minmax[sp$fix, MINMAX_PMAX] >= 1)
        if (!all(fix.zero | fix.one))
          stop("Some parameters could not be fixed. Set 'split_refined = FALSE")
        fix.one <- sp$fix[fix.one & !fix.zero] ### zero has priority
        fix.zero <- sp$fix[fix.zero]
        
        guess[fix.one] <- 1
        guess[fix.zero] <- 0
        
        fixed <- TRUE ## info for the next split[[]] that some variables
        ##               have been fixed. The values must be reported in the
        ##               intermediate results for AIC and logratio test
      }
      next
    }

   
    if (is.list(sp[[1]])) {
     
      res <- recurs.estim(split=sp, level=level+1, splitReg=splitReg,
                          Z = Z,
                          lower=lower,
                          upper=upper,
                          guess= guess,  
                          lsq.methods=lsq.methods,
                          mle.methods=mle.methods,
                          optim.control=optim.control,
                          trafo=trafo, NaNs = NaNs,
                          spConform = spConform,
                          practicalrange = practicalrange,
                          printlevel=printlevel,
                          fit = fit,
                          minmax = minmax,
                          sdvar = sdvar,
			  origin = origin,
                          IDX = IDX)
      if (!is.null(sp$report)) {
        submodels_n <- submodels_n + 1
        if (fixed) res$fixed <- list(zero=fix.zero, one=fix.one)
        submodels[[submodels_n]] <- res
      }
      ##Print("Back from recursive")
    } else { # !is.list(sp[[1]])
      if (printlevel>=PL_RECURSIVE) {
        cat("splitting (",
            paste(rep(". ", level), collapse=""), format(s, width=2),
            ") : x-coords=", paste(sp$x, collapse=","),
            "; compon.=", paste(sp$v, collapse=","), sep="")
        if (printlevel>=PL_RECURSIVE)
          cat( "; parameters=", paste(sp$p, collapse=", "), sep="")
        cat(" ")
      }

      guess <- pmax(lower, pmin(upper, guess, na.rm=TRUE), na.rm=TRUE)

      .C(C_PutValuesAtNAnoInit, splitReg, trafo(guess), NAOK=TRUE) # ok
      old.model <- GetModel(register=splitReg, modus=GETMODEL_DEL_MLE ,
                            return.which.param=INCLUDENOTRETURN,
			    origin = origin)
      guess[sp$p.proj] <- NA
      .C(C_PutValuesAtNAnoInit, splitReg, trafo(guess), NAOK=TRUE) # ok
      
      new.model <- GetModel(register=splitReg, modus=GETMODEL_DEL_MLE,
 			    return.which.param=INCLUDENOTRETURN,
			    origin = origin)

      if (!all(is.na(lower))) {
        .C(C_PutValuesAtNAnoInit, splitReg, trafo(lower), NAOK=TRUE) # ok
        lower.model <- GetModel(register=splitReg, modus=GETMODEL_DEL_MLE,
				return.which.param=INCLUDENOTRETURN,
				origin = origin)
      } else lower.model <- NULL
      
      if (!all(is.na(upper))) {
        .C(C_PutValuesAtNAnoInit, splitReg, trafo(upper), NAOK=TRUE) # ok
        upper.model <- GetModel(register=splitReg, modus=GETMODEL_DEL_MLE,
                                return.which.param=INCLUDENOTRETURN,
                                origin = origin)
      } else upper.model <- NULL


      if ((new.vdim <- length(sp$v.proj)) < Z$vdim) {
        m <- matrix(0, nrow=new.vdim, ncol=Z$vdim)
        for (j in 1:new.vdim) m[j, sp$v.proj[j]] <- 1
        new.model   <- addRMmatrix(new.model, m)
        lower.model <- addRMmatrix(lower.model, m)
        upper.model <- addRMmatrix(upper.model, m)
        old.model   <- addRMmatrix(old.model, m)
        vproj <- rep(FALSE, Z$vdim)
        vproj[sp$v.proj] <- TRUE
        new.data <- lapply(Z$data, function(x) x[ , vproj, drop=FALSE])  
      } else {
        new.data <- Z$data
        vlen <- length(sp$v.proj)
      }

      for (i in 1:length(new.data)) {
        base::dim(new.data[[i]]) <-
          c(Z$coords[[i]]$totpts, vlen,
            length(new.data[[i]]) / (Z$coords[[i]]$totpts * vlen)) 
      }

      x.proj <- sp$x.proj
      ignore.T <- FALSE
      if (!Z$dist.given && !is.logical(sp$x.proj) &&
          length(sp$x.proj) < Z$tsdim) {
        ignored.x <- (1:Z$tsdim)[-x.proj]
        if (Z$coords[[1]]$grid) {
          ignore.T <- TRUE
          new.x <- Z$coords    
          for (i in 1:length(new.data)) {
            #Print(new.x[[i]]) 
            len <- c(new.x[[i]]$x[3, ], if (Z$dist.given) new.x[[i]]$T[3])
            d <- base::dim(new.data[[i]])
            new.d <-  c(len, d[-1]) ## (coord dim, vdim, repet)
            base::dim(new.data[[i]]) <- new.d
            new.data[[i]] <-
              aperm(new.data[[i]], c(x.proj, length(new.d) + (-1:0), ignored.x))
            repet <- d[3] * prod(len[ignored.x]) ## = repet * ignored
            base::dim(new.data[[i]]) <- c(prod(len[x.proj]), d[2], repet)
            new.x[[i]]$x[, ignored.x] <- c(0, 1, 1)
            new.x[[i]]$totpts <- prod(len)
            C_idx <- aperm(array(0:(Z$coords[[i]]$totpts-1), dim=len),
                           c(x.proj, ignored.x))
            new.x[[i]][[XLIST_RAWXIDX]] <-
              as.vector(C_idx)[1:prod(new.x[[i]]$x[3,])]
          }
        } else { ## not grid
          dummy <- new.data
          new.x <- new.data <- list()

          if (Z$has.time.comp) {
            ignore.T <- all(x.proj != Z$tsdim)
            x.proj <- x.proj[x.proj != Z$tsdim]
            ignored.x <- ignored.x[ignored.x != Z$tsdim]
            lenT <- sapply(Z$coords, function(x) x$T[3])
          } else {
            lenT <- rep(1, length(new.data))
          }

          xlen <- sapply(Z$coords, function(x) nrow(x$x))

          for (i in 1:length(dummy)) {
            d <- base::dim(dummy[[i]])
            base::dim(dummy[[i]]) <-
              c(xlen[i], lenT[i], length(dummy[[i]]) / (xlen[i] * lenT[i]))
            xyz <- Z$coords[[i]]$x
            C_idx <- 0:(Z$coords[[i]]$totpts - 1)
            
            ## problem, das hier geloest wird: wenn coordinaten auf 0 gesetzt
            ## werden, muessen die Punkte in diesen Koordinaten uebereinstimmen.
            ## deshalb werden untermengen gebildet, wo dies der Fall ist
            while(length(dummy[[i]]) > 0) {
              slot <- xyz[1, ignored.x]              
              xyz.ignored <- xyz[-1, ignored.x, drop=FALSE]
              if (length(xyz.ignored) > 0) {
                m <- colMeans(xyz.ignored)
                m[m == 0] <- 1 ## arbitrary value -- components are all 0 anyway
                idx <- colSums(abs(t(xyz.ignored) - slot) / m)
                idx <- c(TRUE, idx < min(idx) * 5) ## take also the first one
              } else idx <- TRUE
              last <- length(new.x) + 1
              new.x[[last]] <- Z$coords[[i]]
              new.x[[last]]$x <- xyz[idx, , drop=FALSE]
              new.x[[last]]$x[, ignored.x] <- 0
              new.x[[last]]$totpts <- nrow(new.x[[last]]$x)
              new.data[[last]] <- dummy[[i]][idx, , , drop=FALSE]
              dummy[[i]] <- dummy[[i]][!idx, , , drop=FALSE]
              xyz <- xyz[!idx, , drop=FALSE]
              new.x[[i]][[XLIST_RAWSET]] <- i - 1
              new.x[[i]][[XLIST_RAWXIDX]] <- C_idx[idx]
              C_idx <- C_idx[!idx]
            }
          }
        }
        
      } else {
        new.x <- Z$coords
      }

      
      if (!is.null(Z$transform)) {
        isna <- is.na(Z$transform[[2]](guess))
        idx <- Z$transform[[1]][isna]
        stopifnot(max(sp$p.proj) <= length(guess))
        
        f <- function(p) {
          q <- guess
          q[sp$p.proj] <- p
          z <- (Z$transform[[2]](q))[isna]
          z     
        }
        
        
        new.transform <- list(idx, f)
      } else new.transform <- NULL

      general_spConform <- if (s<length(split)) FALSE else spConform

      if (Z$has.time.comp && ignore.T) {        
        for (i in 1:length(new.x)) new.x[[i]]$T <- double(0)
      }

                                        #     if (length(dim(new.data[[1]])) > 2)
                                        #       dim(new.data[[1]]) <- c(dim(new.data[[1]])[1], prod(dim(new.data[[1]])[-1]))
                                        #       Print(new.model, new.x, new.data, split, s, sp,
                                        #             length(dim(new.data[[1]]))); 
                                        #      stopifnot(length(dim(new.data[[1]])) == 2)
                                        #     readline()
                                        #     stopifnot(s < 1)

      Z0 <- UnifyData(model=new.model, x=new.x, data=new.data,
                      RFopt=getRFoptions())
      Z2$transform <- new.transform

      #Print(Z0)
      
      res <-
        rffit.gauss(Z=Z0,
                    lower= if (!all(is.na(lower))) lower.model,
                    upper= if (!all(is.na(upper))) upper.model,
                    mle.methods=mle.methods,
                    lsq.methods=lsq.methods,                        
                    users.guess=old.model,  
                    optim.control=optim.control,
                    recall = TRUE,
                    fit.split = FALSE,
                    fit.factr = if (s<length(split)) fit$factr_recall
                                else fit$factr,
                    fit.pgtol = if (s<length(split)) fit$pgtol_recall
                                else fit$pgtol,
                    general.practicalrange = practicalrange,
                    general.spConform = general_spConform,
                    sdvar = if (length(sp$v.proj) > 1)
                              sdvar[sp$p.proj, sp$v.proj, drop=FALSE])

      if (!is.null(sp$report)) {
        submodels_n <- submodels_n + 1
        if (general_spConform) {
          res@p.proj <- as.integer(sp$p.proj)
          res@v.proj <- as.integer(sp$v.proj)
          res@x.proj <- sp$x.proj
          res@report <- sp$report
          res@true.tsdim <- as.integer(Z0$tsdim)
          res@true.vdim <- as.integer(Z0$vdim)
          if (fixed) res@fixed <- list(zero=fix.zero, one=fix.one)
        } else {
          res$p.proj <- as.integer(sp$p.proj)
          res$v.proj <- as.integer(sp$v.proj)
          res$x.proj <- sp$x.proj
          res$report <- sp$report
          res$true.tsdim <- as.integer(Z0$tsdim)
          res$true.vdim <- as.integer(Z0$vdim)
          if (fixed) res$fixed <- list(zero=fix.zero, one=fix.one)
        }
        submodels[[submodels_n]] <- res
      }

      table <- if (general_spConform) res@table else res$table  
      np <- length(sp$p.proj)
      
      lower[sp$p.proj] <- pmin(na.rm=TRUE, lower[sp$p.proj],
                               table[[M]][(np + 1) : (2 * np)])
      upper[sp$p.proj] <- pmax(na.rm=TRUE, upper[sp$p.proj],
                               table[[M]][(2 * np + 1) : (3 * np)])
      
      if (s==length(split)) {
        if (submodels_n >= 1) {
          if (general_spConform) res@submodels <- submodels
          else res$submodels <- submodels
        }
        return(res)
      }
    } # else, (is.list(sp[[1]])) 
    
    guess[sp$p.proj] <- res$table[[M]][1:length(sp$p.proj)]    
    fixed <- FALSE
    
  } # for s in split

  return(res)
} # fit$split




######################################################################
##  main function for fitting Gaussian processes
######################################################################

rffit.gauss <- function(Z, lower=NULL, upper=NULL,
                        mle.methods=MLMETHODS,
                        lsq.methods= LSQMETHODS,
                        ## "internal" : name should not be changed;
                        ## should always be last method!
                        users.guess=NULL,  
                        optim.control=NULL,
                        recall = FALSE,
                        sdvar = NULL,
                        ...) {

  ##  print(Z$model)
  ##  Print(Z);
  #print(Z$transform$fctn);  print(Z$transform$params.fctn); 

  ##  Print(lower, recall, upper, users.guess, Z, list(...))
  ##  Print("start", recall)
  ##  str(Z)

  ## ACHTUNG: durch rffit.gauss werden neu gesetzt:
  ##    practicalrange
  ##


  ## some definitions if RFoptions does not deliver enough information:
  Default.nbins <- 20
  Default.nangles <- 5
  Trace.n <- 5
  Trace.plots <- c(3, 5)
  Trace.mar <- c(1.95,1.95,0,0)
  Trace.cex=0.4
  Trace.cex.axis <- 1.5
  Trace.cex.lab <- 3
  Trace.pch <- 1
  Trace.col.emp <- "black"
  Trace.col <- rbind(plain=c("yellow", "orange"),
                     self=c("lightgray", "darkgray"),
                     sqrt.nr=c("tan1", "tan4"),
                     sd.inv=c("peachpuff", "peachpuff3"),
                     internal1=c("springgreen1", "springgreen4"),
                     internal2=c("chartreuse1", "chartreuse4"),
                     internal3=c("palegreen", "palegreen4"),
                     ml=c("red", "darkred")
                     )
  Trace.txt.cex <- 1.7
  Trace.txt.col <- "blue"
  Trace.line=-3

####################################################################
###                      Preludium                               ###
####################################################################

  transform <- Z$transform
  MLE_CONFORM <- if (is.null(transform)) mle_conform else original_model

  par.old <- par(no.readonly=TRUE)
  RFopt <- getRFoptions()
  if (!recall) {
    if (!exists(".Random.seed")) runif(1)
    old.seed <- .Random.seed

    set.seed(0)
    on.exit(.Random.seed <<- old.seed, add = TRUE)
  }

  if (RFopt$general$modus_operandi == MODE_NAMES[normal+1]) Help("normal_mode")
  
  show.error.message <- TRUE # options
  if (!show.error.message) warning("show.error.message = FALSE")

  
  save.options <- options()
  on.exit(options(save.options), add=TRUE)
  
  general <- RFopt$general
  pch <- general$pch
  basic <- RFopt$basic
  printlevel <- basic$printlevel
  if (printlevel < PL_IMPORTANT) pch <- ""
  
  fit <- RFopt$fit
  trace <- fit$trace
  emp_alpha <- fit$emp_alpha
  emp_alphaName <- if (emp_alpha == as.integer(emp_alpha))
                     FCTN_TYPE_NAMES[3 + emp_alpha]
                   else list("Pseudomadogram", alpha=emp_alpha)
  
  ##Print(emp_alpha, emp_alphaName , C_emp_alpha)
  
  if (general$practicalrange && fit$use_naturalscaling)
    stop("practicalrange must be FALSE if fit$use_naturalscaling=TRUE")
  if (fit$use_naturalscaling) setRFoptions(general.practicalrange = 3)
  
  if (printlevel>=PL_STRUCTURE) cat("\nfunction defintions...\n")
  ##all the following save.* are used for debugging only

  detailpch <- if (pch=="") "" else '+'

### definitions from 
  
  LiliReg <- MODEL_MLE ## for calculating likelihood; either it is
  ##                      not overwritten by 'split' or it is not used
  ##                      anymore later on; raw=noConcerns zwingend
  COVreg <- MODEL_LSQ  ## for calculating variogram and getting back models
  ## in particular in LStarget, similar to LiliReg
  ## also used for tracing
  splitReg <- MODEL_SPLIT ## for splitting,
  ## changing x coord, particularly CovLocNonGrid call


  
######################################################################
###                function definitions                            ###
######################################################################
  
  OPTIMIZER <- function(optimiser, max=TRUE) {
    fctn <-
      switch(optimiser,
             "optimx" =
               function(par, fn, lower, upper, control) {
                 TRY(optimx::optimx(par=par, fn=fn, lower=lower, upper=upper,
                                    control=control, method ="L-BFGS-B"))
               },
             "soma" =
               function(par, fn, lower, upper, control) {
                 TRY(soma::soma(if (max) function(...) - fn(...) else fn,
                                bounds=list(min=lower, max=upper),
                                options=list(), strategy="all2one"))
               },
             "nloptr" =
               function(par, fn, lower, upper, control) {
                 if (length(control$xtol_rel) == 0) control$xtol_rel <- 1e-4
                 opts <- control[-pmatch(c("parscale", "fnscale"),
                                         names(control))]
                 opts$algorithm <- fit$nloptr_algorithm
                 TRY(nloptr::nloptr(x0=par, eval_f=fn, lb=lower, ub=upper,
                                    opts= opts))
                 
               },
             "GenSA" =
               function(par, fn, lower, upper, control) {
                 control <- control[-pmatch(c("parscale", "fnscale", "pgtol",
                                              "factr"), names(control))]
                 TRY(GenSA::GenSA(par=par,
                                  if (max) fn=function(...) -fn(...) else fn,
                                  lower=lower, upper=upper, control=control))
               },
             "minqa" =
               function(par, fn, lower, upper, control) {
                 control <- control[-pmatch(c("parscale", "fnscale", "pgtol",
                                              "factr"), names(control))]
                 TRY(minqa::bobyqa(par=par,
                                   if (max) fn=function(...) -fn(...) else fn,
                                   lower=lower, upper=upper,
                                   control=control))
               },
             "pso" =
               function(par, fn, lower, upper, control) {
                 control <- control[-pmatch(c("fnscale", "parscale", "pgtol",
                                              "factr"), names(control))]
                 TRY(pso::psoptim(par=par, fn=fn, lower=lower, upper=upper,
                                  control=control))
               },
             "DEoptim" =
               function(par, fn, lower, upper, control) {
                 control <- control[-pmatch(c("parscale", "fnscale", "pgtol",
                                              "factr"), names(control))]
                 if (length(control)==0)
                   TRY(DEoptim::DEoptim(fn=fn, lower=lower,upper=upper))
                 else
                   TRY(DEoptim::DEoptim(fn=fn, lower=lower,upper=upper,
                                        control=control))
               },
             {
               optimiser <- "optim" # optim as default
               function(par, fn, lower, upper, control) {
                                        # Print(par)
                 TRY(optim(par=par, fn=fn, lower=lower, upper=upper,
                           control=control, method ="L-BFGS-B"))
               }
             }
             )
    ## Print(optimiser, fctn)
    if (optimiser != "optim") {
      if (!requireNamespace(optimiser, quietly = TRUE)) {
        stop("to use '", optimiser, "' its package must be installed")
      }
    }
    fctn
  }

  optim.control                     
  optimiser <- fit$optimiser 
  OPTIM <- OPTIMIZER(optimiser)
  SUBOPTIM <- OPTIMIZER(fit$sub_optimiser, max=FALSE)
  ##if (length(idx <- which("algorithm"==names(control))) > 0) control <- control[-idx]
  
  
  ## optim : standard
  ## optimx: slower, but better 
  ## soma  : extremely slow; sometimes better, sometimes worse
  ## nloptr: viele algorithmen, z.T. gut
  ## GenSA : extrem langsam
  ## minqa : gut, langsamer
  ## pso   : extrem langsam
  ## DEoptim: langsam, aber interessant; leider Dauerausgabe


  

######################################################################
  ##  lower, upper, users.model must match 'model'. This is checked here.
######################################################################
 


  SetUsers <- function(users, own=NULL, type) {
  #  Print(users, own, type)
    default <- switch (type,
                       "lower" = -Inf,
                       "upper" = Inf,
                       "users.guess" = NULL,
                       stop("unknown type"))
    if (!missing(users) && !is.null(users)) {
      if (is.numeric(users)) {
        if (length(users) != length(own))
          stop("number of parameters of '", type,
               "' does not match the model.",
               "Better use the model definition also for '", type,
               "', or the param definition in case of formula notation for the model.")
      } else {
        users <- GetValuesAtNA(NAmodel=Z$model, valuemodel=users,
                               skipchecks=!is.null(trafo), type=type,
                               C_coords = C_coords)
      }

      if (!is.null(own)) {
        idx <- !is.na(users)
        own[idx] <- users[idx]
        users[!idx] <- default
      }
    } else users <- rep(default, length(own))
    
    return(list(users=users, own=own))
  }


  nice_modelinfo <- function(minmax) {
    minmax <- minmax[, -MINMAX_NAN, drop=FALSE]
    minmax <- as.data.frame(minmax)
    minmax$type <- TYPEOF_PARAM_NAMES[minmax$type + 1]
    minmax
  }
  
  INVDIAGHESS <-
    if (!fit$return_hessian)
      function(variab, ...) 
        list(hessian=matrix(NA, ncol=length(variab), nrow=length(variab)),
             sd=rep(NA,length(variab)), invH <- NULL)
    else
      function(variab, control=NULL, MLELB, MLEUB, ...) {
        nvariab <- length(variab)
        if (nvariab == 0) return(list(hessian=NULL, sd=double(0)))
        
        ndeps <- rep(1e-3, nvariab)
        if (length(control) == 0) control <- list()
        control$ndeps <- ndeps
        ndeps <- ndeps * 2 # (1 + 1e-14)
        variabnew <- pmax(MLELB + ndeps, pmin(MLEUB - ndeps, variab))
        idx <- variabnew < MLELB | variabnew > MLEUB ## can happen if MLEUB-MLELB \approx 0
        variabnew[idx] <- variab[idx]
                                        #   idx <- variabnew != variab
        oH <- rawTry(optimHess(par=variabnew, fn=MLtarget, control=control))
        
        sd =rep(NA, nvariab)
        invH <- NULL
        if (is.numeric(oH)) {
          H <- oH
          zaehler <- 1
          while (zaehler <= 3) {
            invH <- rawTry(solve(H))
        if (is.numeric(invH)) {
          var <- -diag(invH)
          var[var < 0] <- Inf
          sd <- sqrt(var)
          break;
        }
            H <- H + rnorm(length(H), 0, 1e-7)
            zaehler <- zaehler + 1;
            invH <- NULL
      }
        } else  H <- matrix(NA, ncol=nvariab, nrow=nvariab)
        
        if (is.na(oH[1]) && !any(bayes))
          warning("The Hessian matrix is strange and not invertable")
        
        return(list(hessian=oH, sd=sd, invH=invH))
      }

  
  show <- function(nr, M, OPT, PARAM)
    cat("\n ", M, ", ", switch(nr, "start", "grid ", "re-do"), ": value=",
        format(OPT, dig=6), ", param=", format(PARAM, dig=2), sep="")

  
  Trace <- function(param, name, trace, opt) {
    ##Print(name)
    ## extern: emp_alpha, all Trace.*
    .C(C_PutValuesAtNA, COVreg, param)
    model.values <- .Call(C_MomentsIntern, COVreg, emp_alpha)
    if (length(Trace.model[[name]]) == 0) {
      Trace.model[[name]] <<- as.matrix(c(model.values))
      Trace.opt[[name]] <<- TRUE
      if (length(Trace.opt) > 1) { ## not the very first one
        prev <- length(Trace.opt) - 1
        idx <- which(Trace.opt[[prev]])
        t <- rep(FALSE, length(Trace.opt[[prev]]))
        t[if (length(idx) >= 1) idx[length(idx)] else length(t)] <- TRUE
        Trace.opt[[prev]] <<- t
        stopifnot(sum(Trace.opt[[prev]]) == 1)
      }
    } else {
      lastopts <- which(Trace.opt[[name]])
      del <- 1 + (trace==2 && length(lastopts) == 1 && lastopts == 1) +
        (length(Trace.opt[[name]]) < Trace.n) * Trace.n # delete only if full
      Trace.model[[name]] <<- cbind(Trace.model[[name]][ , -del], c(model.values))
      Trace.opt[[name]] <<- if (trace==1) c(rep(FALSE, ncol(Trace.model[[name]]) - 1), TRUE)
                            else c(Trace.opt[[name]][-del], opt)
    }

    n <- length(Trace.model)
    act.n <- length(Trace.opt[[name]])
    data <- matrix(nrow=nrow(Trace.model[[name]]), ncol=n-1)
    if (n > 1) for (i in 1:(n-1)) data[, i] <- Trace.model[[i]][, Trace.opt[[i]]]
    data <- cbind(data, Trace.ev, Trace.model[[name]])
    dim(data) <- c(Trace.dim, n + act.n)
    col <- c(if (n > 1) Trace.col[1:(n-1), 2], Trace.col.emp, Trace.col[n, 1 + Trace.opt[[n]]])
    type <- c(rep("l", n-1), "p", rep("l", length(Trace.opt)))
    for (i in 1:Trace.grdim[1]) {
      for (j in 1:Trace.grdim[2]) {
        matplot(ev$centers, data[ , i, j, ], type=type, col=col, pch=Trace.pch,
                cex.axis=Trace.cex.axis#, cex.label=Trace.cex.lab
                )
        if (i * j == 1)
          title(line=Trace.line, cex.main=Trace.txt.cex, col.main=Trace.txt.col,
                main=name, cex.sub=Trace.txt.cex, col.sub=Trace.txt.col,
                sub=if (prod(Trace.dim[2:3]) > 1) paste0(Trace.txt[2], " x ", Trace.txt[3]))
      }
    }
  }

  
  WarningMessage <- function (variab, LB, UB, txt) {
    cat("Note:", txt, ": forbidden values -- if there are too many warnings",
        "try narrower lower and upper bounds for the variables: (",
        paste(variab, collapse=","), ") not in [(",
        paste(LB, collapse=", "),  ") ; (",
        paste(UB, collapse=", "), ")]\n")
  }


  WarningMessage <- function (variab, LB, UB, txt) {
    cat("Note:", txt, ": forbidden values -- if there are too many warnings",
        "try narrower lower and upper bounds for the variables: (",
        paste(variab, collapse=","), ") not all differences to boundardies are non-negative",
        paste(variab - LB, collapse=", "),  ") ; (",
        paste(UB - variab, collapse=", "), ")]\n")
  }

  LSQsettings <- function(M) {    
    assign("LSQ.SELF.WEIGHING", M=="self", envir=ENVIR)
    if (!LSQ.SELF.WEIGHING) {
      assign("LSQ.WEIGHTS", weights[[M]], envir=ENVIR)
      if (globalvariance)
        assign("LSQ.BINNEDSQUARE",
               sum(binned.variogram^2 * LSQ.WEIGHTS, na.rm=TRUE),
               envir=ENVIR)
    }
  }

  
  LStarget <- function(variab) {
    variab <- variab + 0## unbedingt einfuegen, da bei R Fehler
    ##                     der Referenzierung !! 16.2.10
    if (printlevel>PL_FCTN_DETAILS) Print("LS", LSMIN, format(variab, dig=20))#

    ##    Print("begin", LSMIN, format(variab, dig=20))
    
    param <- as.double(trafo(variab))

    if (any((variab<LSQLB) | (variab>LSQUB))) {
      ## for safety -- should not happen, older versions of the optimiser
      ## did not stick precisely to the given bounds
      ## 13.12.03 still happens ...
      if (printlevel>=PL_STRUCTURE)
        WarningMessage(variab, LSQLB, LSQUB, "LSQ")      
      assign("BEYOND", BEYOND + 1, envir=ENVIR)
      variab0 <- pmax(LSQLB, pmin(LSQUB, variab))
      penalty <- sum(variab0 - variab)^2
      save <- list(LSMIN, LSPARAM, LSVARIAB)
      assign("LSMIN", +Inf, envir=ENVIR)
      res <- LStarget(variab0)
      res <- res + penalty * (1 + abs(res))
      if (res <= save[[1]]) {              
        assign("LSMIN", save[[1]], envir=ENVIR)
        assign("LSPARAM", save[[2]], envir=ENVIR)
        assign("LSVARIAB", save[[3]], envir=ENVIR)
      } else assign("LSMIN", res, envir=ENVIR)
      return(res)
    }

    .C(C_PutValuesAtNA, COVreg, param)
    model.values <- .Call(C_MomentsIntern, COVreg, emp_alpha)

    ##    Print(model.values, binned.variogram)
    stopifnot(length(model.values) == length(binned.variogram))

    ##    Print(model.values,(!all(is.finite(model.values))))
    
    if (!all(is.finite(model.values))) {
      
      if (printlevel>=PL_IMPORTANT) {
        message("LSQ missing values!")
      }
      return(1E300)
    }

    if (LSQ.SELF.WEIGHING) {
      ## weights = 1/ model.values^2
      
      gx <- binned.variogram / model.values     
      gx <- gx[is.finite(gx)]
      if (length(gx) == 0)  {
        res <- Inf
        stop("The specified model looks odd.")
      }
      if (globalvariance) {
        bgw <- sum(gx^2)
        g2w <- sum(gx)
        res <- bgw - g2w^2 / length(gx)
      } else {
        res <- sum((gx - 1)^2)
      }
    } else {
      if (globalvariance) {
        ## in this case the calculated variogram model is not the one
        ## to which we should compare, but:
        idx <- is.finite(binned.variogram)
        bgw <- sum(binned.variogram * model.values * LSQ.WEIGHTS, na.rm=TRUE)
        g2w <- sum( (model.values^2 * LSQ.WEIGHTS)[idx])

                                        #        Print(M, idx, bgw, g2w, LSQ.BINNEDSQUARE, bgw^2/g2w)
                                        #       print(binned.variogram)
                                        #      print(model.values)
        ##     print(       LSQ.WEIGHTS   )

        
        res <- LSQ.BINNEDSQUARE - bgw^2/g2w        
      } else {
        res <- sum((binned.variogram - model.values)^2 * LSQ.WEIGHTS,
                   na.rm=TRUE)
      }
    }

    ##    Print("Ende", LSMIN, format(variab, dig=20), res)
    ##    if (res < 0.00000001) q()

    if (res<=LSMIN) {
      assign("LSMIN", res, envir=ENVIR)
      assign("LSPARAM", param, envir=ENVIR)
      assign("LSVARIAB", variab, envir=ENVIR)
    }

    if (trace == 0) {
      if (printlevel>=PL_FCTN_DETAILS) Print(param, LSMIN)  #
    } else {
      if (abs(trace) >= 2 || res == LSMIN) {
        if (trace < 0) Print(M, res, variab, param) #
        else Trace(param=param, name=M, trace=trace, opt = res==LSMIN)
      }
    }
    
    return(res)
  }


  
  MLtarget <- function(variab) {
    ## new version based on C-Code starting wih 3.0.70    
    
    if (n.variab > 0) {
      variab <- variab + 0  ## unbedingt einfuegen, da bei R Fehler der Referenzierung !! 16.2.10
      
      if (printlevel>=PL_FCTN_DETAILS) {
        ##Print(format(variab, dig=20)) #
        if (printlevel>=PL_FCTN_SUBDETAILS) print(minmax) #
      }

      if (any((variab < MLELB) | (variab > MLEUB))) {
        ## for safety -- should not happen, older Uversions of the optimiser
        ## did not stick precisely to the given bounds
        ## 23.12.03 : still happens
        if (printlevel>=PL_STRUCTURE)
          WarningMessage(variab, MLELB, MLEUB, "MLE")   
        assign("BEYOND", BEYOND + 1, envir=ENVIR)
        penalty <- variab
        variab <- pmax(MLELB, pmin(MLEUB, variab)) 
        penalty <-  - sum(variab - penalty)^2 ## not the best ....
        save <- list(MLEMAX, MLECOVAR, MLEPARAM, MLEVARIAB)
        ##    Print(save, variab)
        assign("MLEMAX", -Inf, envir=ENVIR)
        res <- MLtarget(variab)
        res <- res + penalty * (1 + abs(res))       
        if (res < save[[1]]) {
          assign("MLEMAX", save[[1]], envir=ENVIR)
          assign("MLECOVAR", save[[2]], envir=ENVIR)      
          assign("MLEPARAM", save[[3]], envir=ENVIR)
          assign("MLEVARIAB", save[[4]], envir=ENVIR)        
        } else assign("MLEMAX", res, envir=ENVIR)

        return(res)
      }
      
      param <- as.double(trafo(variab))

      ##Print(param, variab)
      .C(C_PutValuesAtNA, LiliReg, param)
      options(show.error.messages = show.error.message)
    } else param <- NULL

    if (printlevel > PL_FCTN_SUBDETAILS) {
      cat("\n\nAufruf von MLtarget\n===================\n")
      Print(variab, param, MLELB, MLEUB) ##
    }

    ans <- TRY(.Call(C_EvaluateModel, double(0), integer(0), LiliReg))
    
    ## e.g. in case of illegal parameter values
    
    if (is(ans, CLASS_TRYERROR) || is.na(ans[1])) {

      if (printlevel > PL_DETAILS) Print(ans) ##
      
      assign("ML_failures", ML_failures + 1)
      if (printlevel > PL_IMPORTANT) ("model evalation has failed")
      if (MLEMAX == -Inf) assign("MLECOVAR", NA, envir=ENVIR)
      return(-1e300)
    }
    res <- ans[1]
    
    ##    Print(ans, MLEMAX, variab, param, MLELB, MLEUB);
    ##print(minmax)
    ##  stopifnot(all(is.finite(MLEUB)))

    if (res >= MLEMAX) {## >= und nicht > da -Inf, -Inf auftreten kann
      ## z.B. bei autostart
      assign("MLEMAX", res, envir=ENVIR)
      assign("MLECOVAR", ans[-1], envir=ENVIR)      
      assign("MLEPARAM", param, envir=ENVIR)
      assign("MLEVARIAB", variab, envir=ENVIR)
    }
    
    if (trace == 0) {
      if (printlevel>=PL_FCTN_DETAILS) Print(ans, MLEMAX)  #
    } else {
      if (abs(trace) >= 2 || res == MLEMAX) {
        if (trace < 0) Print(M, res, variab, param, MLECOVAR) #
        else Trace(param=param, name=M, trace=trace, opt = res==MLEMAX)
      }
    }
    
    return(res)
  } # mltarget


  get.residuals <- function(LiliReg) {
    return( .Call(C_get_logli_residuals, as.integer(LiliReg)))
  }

  get.var.message <- TRUE
  get.var.covariat <- function(Variab) {
    max <- MLEMAX
    covar <- MLECOVAR
    param <- MLEPARAM
    variab <- MLEVARIAB

    assign("MLEMAX", -Inf, envir=ENVIR)

    old.trace <- trace
    trace <<- 0
    MLtarget(Variab)
    result <- MLECOVAR
    trace <<- old.trace
    
    assign("MLEMAX", max, envir=ENVIR)
    assign("MLECOVAR", covar, envir=ENVIR)
    assign("MLEPARAM", param, envir=ENVIR)
    assign("MLEVARIAB", variab, envir=ENVIR)
    if (any(!is.finite(result))) {
      if (printlevel >= PL_IMPORTANT && get.var.message) {
	assign("get.var.message", FALSE, envir=ENVIR)
	message(paste("Variance and/or covariates could not be calculated for some results. They are all set to 1e-8.", if (printlevel <= PL_DETAILSUSER) " Set 'printlevel' at least to 4 to get more information about the reasons."))
      } else cat("%")
      return(1e-8)
    } else return(result)
  }

  
  ## to avoid warning on "no visible binding" we define the following
  ## variable that are used in the local functions:
  ENVIR <- environment()
  LSQ.SELF.WEIGHING <- LSQ.WEIGHTS <- LSQ.BINNEDSQUARE <- 
    DO.REML <- DO.REML1 <- RML.A <- RML.data <-
      REML.CORRECTION <- DO.RML1 <- 
        ML.RESIDUALS <- MLEMAX <- MLEINF <- MLECOVAR <- MLEPARAM <- 
          CROSS.DIST <- CROSS.KRIGE <- CROSS.VAR <- CROSSMODEL <-
            LOGDET <- NULL
  BEYOND <- 0
  ML_failures <- 0

  
######################################################################
###              End of definitions of local functions             ###
######################################################################


######################################################################
###    Initial settings, coords part I (without model info)         ###
######################################################################

  if (printlevel>=PL_STRUCTURE) cat("\ninitial settings...\n")

  dist.given <- Z$dist.given
  spatialdim <- Z$spatialdim
  time <-  Z$has.time.comp
  xdimOZ <- Z$xdimOZ
  vdim <- Z$vdim
  if (vdim > 1 && emp_alpha < 0 && emp_alpha != VARIOGRAM)
    stop("madograms cannot be used in the multivariate case")
  
  coords <- Z$coords 
  C_coords <- Z$C_coords
  tsdim <- Z$tsdim
  mindistance <- Z$mindist
  maxdistance <- Z$maxdist

################    analyses of orginal model        ###############
##### variables needed for analysis of trend, upper and lower input --
##### user cannot know what the internal represenatation is

###
##  Print(Z)

  if (printlevel>=PL_STRUCTURE) cat("\nfirst analysis of model  ...\n")


  info.cov <- .Call(C_SetAndGetModelLikelihood, LiliReg,
                    list("RFloglikelihood", data = Z$data, Z$model,
                         standardize=fit$standardizedL),
                    C_coords, noConcerns, MLE_CONFORM) ## LiliReg wieder verwendet !!!
  
  
  trans.inv <- info.cov$trans.inv ## note: only with respect to the
  ##              coordinates, mixed effect koennen andere wirkung haben
  isotropic <- info.cov$isotropic
  minmax<- info.cov$minmax # 4 Spalten: 1:min, 2:max, 3:type, 4:is.nan, not na
  NaNs <- as.logical(minmax[, "NAN"])
  coordsystem <- minmax[, "iso"][!NaNs]
  minmax.names <- attr(minmax, "dimnames")[[1]]
  variabnames <- minmax.names[!NaNs]
  n.param <- nrow(minmax)
  n.variab <- n.param - info.cov$NaNs
  Z$effect <- effect <- info.cov$effect # length == number of summands in "+"s
  if (!any(effect >= RandomEffect))
    stop("there must be an error component in the model")

  if (is.null(likeli.info <- .Call(C_get_likeliinfo, LiliReg)))
    stop("bug in likelihood. Please inform author.")
  n.covariat <- likeli.info$betas
  betanames <- likeli.info$betanames
  globalvariance <- likeli.info$estimate_variance
  sum.not.isna.data <- likeli.info$sum_not_isna_data


  ##  Print(n.param, info.cov)
  ##  print(minmax)
  
  stopifnot(n.param == sum(info.cov$NAs)) ## NAs per model
  
  
  ##

  ## hier zum ersten mal model verwendet wichtig,
  ## da interne Darstellung abweichen kann. Z.B. dass ein optionaler Parameter
  ## auf einen Standardwert gesetzt wird
  ## -- wichtig fuer u.a. GetValuesAtNA

### ACHTUNG: Z hat nicht VAR immer noch auf NA (falls globalvariance!!) ?????

  
  if (any(sapply(coords, function(x) x$dist.given)) &&
      (!trans.inv || !isotropic))
    stop("only domain and isotropic models go along with distances")
  
  
  coordnames <- c("cartesian", "sphere", "earth")
  coordsystem <- coordnames[1 + (coordsystem >= FIRST_SPHERICAL) +
                            (coordsystem >= FIRST_EARTH)]


######################################################################
  ## trafo
######################################################################
  if (printlevel>=PL_STRUCTURE) cat("\ntrafo...")
  
  params.fctn <- trafoidx <- trafo <- NULL
  if (is.null(transform)) {
    if (any(minmax[, MINMAX_NAN]==1)) ## nan, not na
      stop("NaN only allowed if transform is given.")
    trafo <- function(x) x;
  } else {
    if (!recall) message("\nNote: if 'transform' or 'params' is given and off-diagonal elements of an anisotropy matrix are estimated, then 'RMS' should be given explicitely and 'anisoT' instead of 'Aniso' should be used within 'RMS'.")

    if (trafo <- is.list(transform) && length(transform)==3) {
      trafoidx <- transform$isNA
      params.fctn <- transform$params.fctn

      if (n.variab == sum(trafoidx) && n.param == length(trafoidx))
        trafo <- transform$fctn
    }

    if (printlevel > PL_IMPORTANT) print(minmax[!NaNs, ]) #
    
    if (!is.function(trafo)) {
      values <- rep(NaN, n.param)
      values[!NaNs] <- as.double((1:n.variab) * 1.111)

      .C(C_PutValuesAtNA, LiliReg, values, NAOK=TRUE)
      
      model_with_NAs_replaced_by_ranks <-
        GetModel(register=LiliReg, modus=GETMODEL_DEL_MLE,
                 which.submodels="user.but.once+jump", origin = original_model)
      cat("'Transform must be a list of two elements:\n* A vector V of logicals whose length equals the number of identified NAs in the model (see below).\n* A function taking a vector whose length equals sum(V). The output\n  is the original model where the NAs are replaced by real numbers.\n\nThe currently identified", n.variab,  "NAs are (positions are given by n.nnn) :\n")
      str(model_with_NAs_replaced_by_ranks, vec.len=20, digits=4) #
      L <- length(trafoidx)
      cat("\nHowever, the fitting method was called by",
          if (L < n.variab) paste("only", L)
          else if (L > n.variab) paste("as many as", L)
          else paste("the following ", L),
          "variables",
          ":\n")
      print(transform$isNA, vec.len=10) #      
      cat("('TRUE' means 'to be estimated'.)\n")
      if (L > n.variab)
        cat("\nNote that the parameters of the trend are optimised analytically, hence they may not be considered, here.\n")
      else if (L < n.variab)
        cat( sep="","\nLikely too many variables have been declared",
            " (maybe by '", RM_DECLARE, "').\n", 
            "Else: if you use argument 'transform' explicitely, check it.\n",
            "Else: inform the maintainer about this error.\n")
      return(NULL)
    }

    if (any(minmax[!trafoidx, MINMAX_NAN] != 1) ||
        any(minmax[trafoidx, MINMAX_NAN] != 0))
      stop("NaNs do not match logical vector of transform")
    minmax <- minmax[trafoidx, , drop=FALSE]    
  }
  if (printlevel >= PL_SUBIMPORTANT + recall) print(minmax) #
  
  ptype <- minmax[, MINMAX_TYPE]
  diag.idx <- which(ptype == DIAGPARAM)
  if (length(diag.idx)>0) minmax[diag.idx[1], MINMAX_PMIN] <- fit$min_diag

  if (Z$matrix.indep.of.x.assumed && !info.cov$matrix.indep.of.x)
    stop("x-coordinates are neither given by 'x' nor by 'distances' nor by 'data',\n  but the model seem to require them")
  

  ts.xdim <- as.integer(xdimOZ + time)
  sets <- length(Z$data)
  repet <- Z$repetitions   
  N <- S <- Sq <- 0 
  for (i in 1:sets) {
    N <- N + rowSums(matrix(colSums(!is.na(Z$data[[i]])), ncol=repet[i]))
    S  <- S  + rowSums(matrix(colSums(Z$data[[i]], na.rm=TRUE),
                              ncol=repet[i]), na.rm=TRUE)    
    Sq <- Sq + rowSums(matrix(colSums(Z$data[[i]]^2, na.rm=TRUE),
                              ncol=repet[i]), na.rm=TRUE)
  }
  if (vdim == 1) {
    N <- sum(N)
    S <- sum(S)
    Sq <- sum(Sq)
  }
  mean.data <- S / N
  var.data <- Sq / N - mean.data^2
  



#######################   upper,lower,user     ########################
  if (printlevel>=PL_STRUCTURE) cat("\nlower and upper ...\n")
  
  SU <- SetUsers(lower, minmax[, MINMAX_PMIN], "lower")
  users.lower <- SU$users
  lower <- SU$own
  
  SU <- SetUsers(upper, minmax[, MINMAX_PMAX], "upper")
  users.upper <- SU$users
  upper <- SU$own

  if (any(idx <- users.upper <= users.lower))
    stop("Lower bounds must be strictly smaller than the upper ones, ",
         "what is not true for ",
         paste("'", variabnames[idx], "'",  sep="", collapse=", "), ".")

  ## nur vollstaendige "guesses" erlaubt
  users.guess <- SetUsers(users.guess, NULL, "users.guess")$users

  
###########################  upper, lower, transform     #######################
  ## either given bu users.transform + users.min, users.max
  ## DIESER TEIL MUSS IMMER HINTER SetUsers STEHEN
  if (printlevel>=PL_STRUCTURE) cat("\ntransform ...\n")
  
  ##  delete.idx <- rep(FALSE, length(lower))
  for (lu in c("lower", "upper")) {
    z <- TRY(trafo(get(lu)))
    ##    print(trafo); Print(lu, get(lu), z, minmax)
    if (!is.numeric(z) || !all(is.finite(z)) || !is.vector(z)) {
      stop("A '", lu, "' bound is either not explicitely set or is not sound",
           if (is(z, CLASS_TRYERROR)) paste(":", z$message) else ".")
    }
                                        # Or the dummy variables have not been defined by '", RM_DECLARE, "'.")
    if (length(z) != n.param)
      stop("\n'transform' returns a vector of length ",
           length(z), ", but one of length ", n.param, " is expected.",
           if (length(z) > n.param) " Note that the parameters of the trend are optimised analytically, hence they may not be considered, here.",
           " Call 'RFfit' with `transform=list()' to get more information on the parameters of the model.\n\n")
  }

  
  ## Achtung which, da upper,lower etc um box-cox-Variable verlaengert
  ## werden koennten !
  SDVAR.IDX <- ptype == SDPARAM | ptype == VARPARAM | ptype == NUGGETVAR
  SIGN.VAR.IDX <- ptype == SIGNEDVARPARAM
  SIGN.SD.IDX <- ptype == SIGNEDSDPARAM
  ALL.SDVAR <- SDVAR.IDX | SIGN.VAR.IDX | SIGN.SD.IDX

  
  if (is.null(sdvar)) {
    sdvar <- matrix(FALSE, nrow=n.variab, ncol=vdim)
    for (i in 1:vdim) {
      sdvar[ , i] <- (minmax[, MINMAX_COLS] == i |
                      minmax[, MINMAX_ROWS] == i) & SDVAR.IDX
      ##                     minmax[, MINMAX_ROWS] == i) & ANY.SDVAR
      ##Print(i, MINMAX_COLS, MINMAX_ROWS, minmax[, MINMAX_COLS],  minmax[, MINMAX_ROWS] == i, SDVAR.IDX)
    }    
  } else stopifnot(all(rowSums(sdvar[SDVAR.IDX, ]) >= 1))
  ##  } else stopifnot(all(rowSums(sdvar[ANY.SDVAR, ]) >= 1)) {

  SCALE.IDX <- ptype == SCALEPARAM  ## large capitals 
  var.idx <- which(ptype == VARPARAM)
  sd.idx <- which(ptype == SDPARAM)
  nugget.idx <- which(ptype == NUGGETVAR)

  if (vdim ==1) {
    varmin <- varmax <- rep(var.data, n.variab)
  } else {
    ## sdvar : matrix of indices
    varmax <- apply(sdvar, 1, function(x) if (any(x)) max(var.data[x]) else NA)
    varmin <- apply(sdvar, 1, function(x) if (any(x)) min(var.data[x]) else NA)
    if (printlevel >= PL_IMPORTANT && !recall) {
      mean.var <- mean(var.data)
      if (mean.var == 0) message("all variables are constants")
      else if (mean(abs(mean.data)) != 0 &&
               any(log(sd(mean.data) / mean(abs(mean.data))) > 1.5))
        message("Are the average values of the components rather different? If so, it might be\n worth thinking of standardising the values before calling RFfit.\n")
      else if (any(abs(log(var.data / mean.var)) > 2.0))
        message("The standard deviations of the components are rather different. It might be\n better to standardise the components of the data before calling RFfit.\n")
    }
  } # vdim > 1
  
#################################################################
##############     prepare constants in S, X,etc      ###########
#################################################################
  if (printlevel>=PL_STRUCTURE) cat("\ndistances and data...")

##############         distances              #################
  ## to do: distances auf C berechnen falls vorteilhaft!


  
##############         Coordinates & data    #################
  ## note: the direct C call needs matrix where points are given column-wise
  ##       whereas the R function CovarianceFct need them row-wise,
  ##                   except for fctncall==CovarianceMatrix
  

  
  ## to do: missing values -- should be distinguished between
  ## lots of missing and not that many?
  ## balanced <- rep(TRUE, sets)
  

  if (vdim>1 && printlevel>=PL_IMPORTANT && !recall)
    message("Due to the covariance model a ", vdim,
            "-variate random field is expected. Therefore, \nthe data matrix",
            " is assumed to consist of ", repet,
            " independent measurements for\neach point.",
            " Each realisation is given as the entries of ", vdim,
            " consecutive \ncolumns.")

  
##############      find upper and lower bounds      #################
  if (printlevel>=PL_STRUCTURE) cat("\nbounds...")
  
  txt <- "lower and upper are both lists or vectors of the same length or NULL"
  lengthmismatch <- "lengths of bound vectors do not match model"
  structuremismatch <- "structures of bounds do not match the model"

  
  ## autostart will give the starting values for LSQ
  ## appears if trafo is given. Then better do not search for
  ## automatic bounds

  autostart <- numeric(length(lower))
  neg <- lower <= 0
  autostart[neg] <-  0.5 * (lower[neg] + upper[neg])
  autostart[!neg] <- sqrt(lower[!neg]*upper[!neg])

  idx <- which(SDVAR.IDX)

  if (length(idx) > 0) {
    ## lower bound of first model is treated differently!
    ## so the "main" model should be given first!             !!!!!
    
    
    ## lower[idx] <- 0
    ## first.idx <- nugget.idx
    ## if (is.null(first.idx)) first.idx <- var.idx
    ## if (is.null(first.idx)) first.idx <- sd.idx
    lower[idx] <- varmin[idx] / fit$lowerbound_var_factor / length(idx)
                                        #  lower[idx] <- 0 ## ??
    
    if (fit$lowerbound_var_factor == Inf && length(idx)>1) {
      idx2 <- which(users.guess[idx] == max(users.guess[idx], na.rm=TRUE))
      if (length(idx2) == 0) idx2 <- 1
      idx2 <- idx[idx2[1]]
      lower[idx2] <- varmin[idx2] / 1e8
    }
    
    upper[idx] <- varmax[idx] * fit$upperbound_var_factor
    autostart[idx] <- sqrt(varmin[idx] * varmax[idx]) /length(idx)
    
    if (length(sd.idx) > 0) {
      lower[sd.idx] <- sqrt(lower[sd.idx])
      upper[sd.idx] <- sqrt(upper[sd.idx])
      autostart[sd.idx] <- sqrt(autostart[sd.idx])
    }       
  }
  
  SIGN.VAR.IDX <- SIGN.VAR.IDX & !is.na(varmax)
  if (any(SIGN.VAR.IDX)) {
    lower[SIGN.VAR.IDX] <-
      -(upper[SIGN.VAR.IDX] <- varmax[SIGN.VAR.IDX]*fit$upperbound_var_factor);
    autostart[SIGN.VAR.IDX] <- 0 
  }
  
  SIGN.SD.IDX <- SIGN.SD.IDX & !is.na(varmax)
  if (any(SIGN.SD.IDX)) {
    lower[SIGN.SD.IDX] <-
      -(upper[SIGN.SD.IDX]<-sqrt(varmax[SIGN.SD.IDX]*fit$upperbound_var_factor))
    autostart[SIGN.SD.IDX] <- 0 
  }
  
  lb.s.ls.f <- fit$lowerbound_scale_ls_factor
  up.s.f <- fit$upperbound_scale_factor
  if (RFopt$coords$earth_coord_names[1] %in% coords[[1]]$coordunits) {
    if ("km" %in% coords[[1]]$new_coordunits) {
      lb.s.ls.f <- lb.s.ls.f / 40 # 40 = approx 7000 / 180
      up.s.f <- up.s.f * 40
    } else if ("miles" %in% coords[[1]]$new_coordunits) {
      lb.s.ls.f <- lb.s.ls.f / 20 # 20 = approx 7000 / 180
      up.s.f <- up.s.f * 20      
    } else stop("unknown coordinate transformation")
  }
  if (any(idx <- ptype == DIAGPARAM)) {
    lower[idx] <- 1 / (up.s.f * maxdistance)
    upper[idx] <- lb.s.ls.f / mindistance
    autostart[idx] <- 8 / (maxdistance + 7 * mindistance)
  }

  
  if (any(idx <- ptype == ANISOPARAM)) {
    if (is.null(trafo))
      warning("The algorithms RandomFields transpose the matrix Aniso to aniso -- this may cause problems when applying transform to the anisotropy parameters. To be safe, use only the parameter anisoT in RMfit.")
    lower[idx] <- -lb.s.ls.f / mindistance
    autostart[idx] <- 0
  }

  if (any(SCALE.IDX)) {
    idx <- which(SCALE.IDX)
#    Print(lower, idx, mindistance, lb.s.ls.f)
    lower[idx] <- mindistance / lb.s.ls.f
    upper[idx] <- maxdistance * up.s.f
    autostart[idx] <- (maxdistance + 7 * mindistance) / 8      
  }

  
###########################        split       #######################
  if (printlevel>=PL_STRUCTURE) cat("\nsplit...")
  
  if (FALSE && fit$split > 0 && length(autostart)>=fit$split) {
    stopifnot(fit$split > 1)    
    
    new.param <- if (is.null(users.guess)) autostart else users.guess


### Achtung!! delete.idx darf davor nur fuer trafo gesetzt werden!!

    stopifnot(spatialdim == tsdim - time)
    if (length(C_coords[[1]]$y)>0) stop("x/y mismatch.", CONTACT)
    splitxy <- matrix(as.double(1:(spatialdim^2)), ncol = spatialdim)
    split_l <- spatialdim
    if (dist.given) {
      if (xdimOZ == 1) {
        splitxy <- t(as.vector(dist(splitxy)))
        split_l = length(splitxy)
      } else stop("vector-valued distances currently not allowed")
    }
     
    .Call(C_SetAndGetModelLikelihood, splitReg, list("Covariance", Z$model), 
          C_coords, neverGlobalXT, MLE_CONFORM)
    .Call(C_LocNonGrid, COVreg,
          list(x=splitxy,                 #0
               y=if (!dist.given) splitxy else double(0),#1
               T=C_coords[[1]]$T,          #2
               grid = FALSE,              #3
               spatialdim = spatialdim,   #4
               has.time.comp = C_coords[[1]]$has.time.comp,  #5
               dist.given = dist.given,   #6 ok
               totpts = nrow(C_coords[[1]]$x),     #7
               l = spatialdim,            #8
               coordunits = C_coords[[1]]$coordunits,
               new_coordunits = C_coords[[1]]$new_coordunits))

    
       
    stopifnot(ncol(Z$rangex) == ts.xdim)
     
    split <- TRY(ModelSplit(splitReg=splitReg, info.cov=info.cov,trafo=trafo,
                            variab=new.param,
                            lower=lower, upper=upper,
                            rangex = Z$rangex,
                              ## ts.xdim != tsdim falls distances (x has.time.compy)
                            modelinfo=list(ts.xdim=ts.xdim, tsdim=tsdim,
                                           xdimOZ = xdimOZ, vdim=vdim,
                                           dist.given=dist.given,
                                           refined = fit$split_refined),
                            model=Z$model))
    
      
    if (is(split, CLASS_TRYERROR)) {
      message("Splitting failed (", split[[1]], "). \nSo, standard optimization is tried")
    } else {
      if (printlevel>=PL_STRUCTURE) cat("\nsplitted (", length(split), ") ...")

      if (length(split) == 1)
        stop("split does not work in this case. Use split=FALSE")
      if (length(split) > 1) {
        if (printlevel >= PL_RECURSIVE) {
          if (printlevel>=PL_STRUCTURE) cat("\n")
          Print(split) #
        } 
        
# 	Print("spitting. call start", C_Set AndGetModelLikelihood, splitReg, list("Covariance", Z$model),  C_coords)
        .Call(C_SetAndGetModelLikelihood, splitReg,
              list("Covariance", Z$model),
              C_coords, noConcerns, MLE_CONFORM) ## not in the previous version


#	print("spitting. call end")
        
        
        return(recurs.estim(split=split, level=0,  splitReg=splitReg,
                            Z = Z,
                            lower= if (is.null(transform)) rep(NA, n.variab)
                                   else lower,
                            upper= if (is.null(transform)) rep(NA, n.variab)
                                   else upper,
                            guess=new.param, # setzt default werte
                            lsq.methods = LSQMETHODS,
                            mle.methods=mle.methods,
                            optim.control=optim.control,
                            trafo=trafo,
                            NaNs=NaNs,
                            spConform = general$spConform,
                            practicalrange = general$practicalrange,
                            printlevel=printlevel,
                            minmax=minmax,
                            fit = RFopt$fit,
                            sdvar = apply(split[[1]]$all.p, 2,
                                function(x) {x[!ALL.SDVAR] <- FALSE; x}),
                            ## sdvar in spaltenrichtung vdim, in zeilenrichtung
			   ## die parameter
                            origin = MLE_CONFORM,
                            IDX = list("lower"=(n.variab + 1) : (2 * n.variab),
                                       "upper"=(2*n.variab+1) : (3 * n.variab))
			   ))
      }
    }
  }

###  delete.idx <- which(delete.idx) ## !!
  
 
######################################################################
###                                                                ###
###   check which parameters are NA -- only old style is allowed   ###
###                                                                ###
###     certain combinations of NA allow for faster algorithms     ###
###                                                                ###
###     !is.na(sill) needs special treatment, hence must be        ###
###     identified --- currently deleted; see before 3.0.70        ###
###                                                                ###
###                                                                ###
###     scaling method must be identified                          ###
###                                                                ###
###     some autostart values are calculated                       ###
###                                                                ###
###                                                                ###
######################################################################

  if (printlevel>=PL_STRUCTURE) cat("\nauto...")  
  
  ssm <-  nonugget <- novariance <- FALSE
  var.global <- var.idx

  
  idx <- which(users.lower > -Inf)
  probably.wrong <- users.lower[idx] > upper[idx]
  for (bound in c("lower", "upper")) {
    if (any(probably.wrong)) {
      probably.wrong <- idx[probably.wrong]
      b <- if (bound == "lower") users.upper else users.lower
      txt <- paste("The user given", bound, "bound of",
                   paste("'", variabnames[probably.wrong], "'", sep="",
                         collapse = ", "),
                   "is", if (bound == "lower") "larger than the upper"
                         else "smaller than the lower",
                   "bound calculated internally.\n")
      if (all(is.finite(b[probably.wrong])))
        message(txt, "Although considered as wrong the user bound are taken.")
      else
        stop(txt,"In this case the upper bounds must also be given and finite.")
      upper[probably.wrong] <- users.upper[probably.wrong] # notwendig trotz
      ## SetUsers, da z.T.us daten ermittelt unabhaengig vom Vorwert
      lower[probably.wrong] <- users.lower[probably.wrong]
    }
    ## loop incremental part:
    idx <- which(users.upper < Inf)
    probably.wrong <- users.upper[idx] < lower[idx]
  }

  if (!fit$suggesting_bounds) { ## i.e. users bound are always taken over
    idx <- which(users.lower > -Inf)
    lower[idx] <- users.lower[idx]
    idx <- which(users.upper < Inf)
    upper[idx] <- users.upper[idx]
  }


  bounds <- minmax[, c(MINMAX_MIN, MINMAX_MAX), drop=FALSE]
  bayes <- as.logical(minmax[, MINMAX_BAYES])
   if (any(bayes)) {
    lower[bayes] <- pmax(lower[bayes], minmax[bayes, MINMAX_PMIN])
    upper[bayes] <- pmin(upper[bayes], minmax[bayes, MINMAX_PMAX])    
    bounds[bayes, 1] <- pmax(bounds[bayes, 1], minmax[bayes, MINMAX_PMIN])
    bounds[bayes, 2] <- pmin(bounds[bayes, 2], minmax[bayes, MINMAX_PMAX])
  }

  bounds <- apply(abs(bounds), 1, max)

## aelter als 21.11.20 ausgeblendet:    
##  if (length(delete.idx)>0) {
##    upper <- upper[-delete.idx]
##    lower <- lower[-delete.idx]
## #   users.lower <- users.lower[-delete.idx]
## #   users.upper <- users.upper[-delete.idx]
##    autostart <-autostart[-delete.idx]
##    variabnames <- variabnames[-delete.idx]
##    SCALE.IDX <- SCALE.IDX[-delete.idx]
##    SDVAR.IDX <- SDVAR.IDX[-delete.idx]
##    ptype <- ptype[-delete.idx]
##    bounds <- bounds[-delete.idx]    
  ##  }

  if (any(autostart<lower) || any(autostart>upper)) {
    if (printlevel >= PL_ERRORS) Print(cbind(lower, autostart, upper)) #
    autostart <- pmin(upper, pmax(lower, autostart))
  }

  
  if (n.variab == 0) { ## diese Abfrage erst spaet wegen globalvariance
    if (globalvariance) Help("onlyvar")
    else Note("no_fit", n.covariat)
  }

  if (any(idx <-lower >= upper)) {
    lu <- cbind(lower=lower, upper=upper, idx)  
    stop(paste("Some lower bounds for the parameters are greater than ",
               "or equal to the upper bound\n",
               paste(collapse="\n ", dimnames(lu)[[1]], ":",
                     apply(lu, 1, function(x)
                           paste("\tlower=", signif(x[1]),
                                 ",\tupper=", signif(x[2]),
                                 if (x[3]) "  \tnot ok!", sep=""))
                     )
               ))
  }


  fill.in <-  trafo(autostart)
#  .C(C_PutValuesAtNA, R e g, trafo(autostart))   

  
######################################################################
######################################################################
###                     Estimation part itself                     ###
######################################################################
######################################################################

  ## check optim.control 
  ## parscale will give the magnitude of the parameters to be eliminated
  ##     passed to optim/optimise so that the optimiser eliminates
  ##     values around 1 ##to do: irgendwo in paper die wichtigkeit beschreiben?

  parscale <- ParScale(optim.control, lower=lower, upper=upper) 
  fit.fnscale <- optim.control$fnscale

  
  if (length(optim.control)>0) {
    opt.control <- optim.control
    stopifnot(is.list(opt.control))
    forbidden.param <- c("parscale", "fnscale", "algorithm")
    ## fnscale=-1 turns the problem into a maximisation problem, see below
    forbidden <- which(!is.na(pmatch(names(opt.control), forbidden.param)))
    forbidden.opt <- opt.control[forbidden]  
    if (length(forbidden) > 0)  opt.control <- opt.control[-forbidden]
   } else {
    opt.control <- list()
  }

  if (length(fit$algorithm) > 0 && fit$algorithm != "")
      opt.control$algorithm <- fit$algorithm
  if (length(optim.control)>0) {
    if (length(forbidden.opt$algorithm) > 0)
      opt.control$algorithm <- forbidden.opt$algorithm
  }

  if (fit$optimiser=="optim") {
    if (length(opt.control$pgtol)==0) #not given by user
      opt.control$pgtol <- fit$pgtol
   if (length(opt.control$factr)==0) #not given by user
      opt.control$factr <- fit$factr
  }


###################  preparation  ################
  if (printlevel>=PL_STRUCTURE) cat("\npreparing fitting...")
  ## methods
  formals <- formals()
              
  if (length(MLMETHODS) != 1 || MLMETHODS != "ml")
    stop("reml currently not programmed")

  
  allmethods <- c(PRIMMETHODS, LSQMETHODS, MLMETHODS, CROSSMETHODS)

  ## how preceding methods have been considered ?
  ## note cm is used again at the very end when error checking
  cm <- cumsum(c(0, length(PRIMMETHODS), length(LSQMETHODS),
                     length(MLMETHODS), length(CROSSMETHODS)))
  cm <- cbind(cm[-length(cm)] + 1, cm[-1])
  cm <- apply(cm, 1, function(x) x[1] : x[2])
  names(cm) <- c("prim", "lsq", "mle", "cross")


   methodprevto <-
    if (fit$only_users) {
      list(lsq="users.guess",mle="users.guess",cross="users.guess")
    } else list(lsq=c(cm$prim),
		mle=c(cm$prim, cm$lsq),
		cross=c(cm$prim, cm$lsq, cm$cross)
		)

  ## index (start, end) to the various categories of
  ## information to be stored
  IDX <- function(name) { idx <- tblidx[[name]]; idx[1]:idx[2]}

  tblidx <- cumsum(c(0,
                     n.variab, # variables used in algorithm
                     n.variab, # their lower bounds
                     n.variab, # ... and upper bounds
                     n.variab,  # sd of variabs
                     n.param,# param values to be eliminated
                     rep(1, length(allmethods) - length(PRIMMETHODS)),#method
                     ##                                                 score
                     1, 1, 1, ## AIC, AICc, BIC,
                     as.integer(globalvariance),
                     n.covariat, #
                     # whether there has been a global variance to estimated
                     # coeff to eliminated for covariates, i.e.
                     ##           trend parameters
                     n.param ## sd of params
                    ))

  if (printlevel>=PL_STRUCTURE) cat("\npreparing fitting (part 2)...")

  tblidx <- rbind(tblidx[-length(tblidx)] + 1, tblidx[-1])
  idx <- tblidx[1, ] > tblidx[2, ]
  tblidx[, idx] <- 0 
   if (printlevel>=PL_STRUCTURE) cat("\npreparing fitting (part 3)...")

   
  dimnames(tblidx) <- list(c("start", "end"),                           
                           c("variab", "lower", "upper", "sdvariab",
                             "param", 
                             allmethods[-1:-length(PRIMMETHODS)],
                             "AIC", "AICc", "BIC",
                             "glbl.var",
                             "covariat",
                             "sdparam"
                             ##,  "lowbeta", "upbeta", only used for
                             ## cross-validation
                             ))
  if (tblidx[2, "covariat"] != 0) tblidx[2, "glbl.var"] <- tblidx[2, "covariat"]
  if (tblidx[1, "glbl.var"] == 0) tblidx[1, "glbl.var"] <- tblidx[1, "covariat"]

  if (printlevel>=PL_STRUCTURE) cat("\npreparing fitting (part 4)...")

  maxtblidx <- max(tblidx)
  tblidx <- data.frame(tblidx)

  if (printlevel>=PL_STRUCTURE) cat("\npreparing fitting (part 5)...")

  ## table of all information; col:various method; row:information to method
  tablenames <- if (n.variab == 0) NULL
                else c(paste("v", variabnames, sep=":"),        
                       paste("lb", variabnames, sep=":"),
                       paste("ub", variabnames, sep=":"),
                       paste("sdv", variabnames, sep=":"))
    
  tablenames <-
    c(tablenames,
      minmax.names,
      ##       if (n.param > 0) paste("p", 1:n.param, sep=":")
      ##                       else minmax.names,
      allmethods[-1:-length(PRIMMETHODS)],
      "AIC", "AICc", "BIC",
      if (globalvariance) "glbl.var",  
      betanames,
      ## do not try to join the next two lines, since both
      ## variabnames and betanames may contain nonsense if
      ## n.variab==0 and n.covariat==0, respectively
      if (n.variab > 0)  paste("sd", minmax.names, sep=":")
      )

  param.table <- data.frame(matrix(NA, nrow=maxtblidx, ncol=length(allmethods),
                                   dimnames=list(tablenames, allmethods)))

  fit.fnscale <- if (is.null(fit.fnscale)) rep(NA, length(allmethods)) else
                  -abs(fit.fnscale)


#############################################################
## end preparation; remaining part is elimination  ###########
#############################################################

##################################################
###############    PRIMITIVE METHODS   ###########
##################################################
 
  MLELB <- LSQLB <- lower
  MLEUB <- LSQUB <- upper


 
  ##****************    autostart    *****************
  if (printlevel>=PL_STRUCTURE) cat("\nautostart...")
  M <- "autostart"
  primMethods <- M
  default.param <- param.table[[M]][IDX("variab")] <- autostart

  param.table[[M]][IDX("param")] <- trafo(autostart) 
  MLEVARIAB <- autostart

  param.table[[M]][IDX("glbl.var")] <- get.var.covariat(autostart)

  ## ****************    user's guess    *****************
  if (!is.null(users.guess)) {
    M <- "users.guess"
    primMethods <- c(primMethods, M)
#    if (length(delete.idx) > 0) users.guess <- users.guess[-delete.idx]
    idx <- users.guess < lower | users.guess > upper
    if (any(idx)) {
      if (recall) users.guess <- NULL
      else if (general$modus_operandi ==  MODE_NAMES[careless + 1]) {
        lower <- pmin(lower, users.guess)
        upper <- pmax(upper, users.guess)
        MLELB <- LSQLB <- lower
        MLEUB <- LSQUB <- upper
     } else {
        m <- cbind(lower, users.guess, upper, idx)
        dimnames(m) <- list(rep("", length(lower)),
                            c("lower", "user", "upper", "outside"))
        cat("\n")
        print(m) ## nicht loeschen!
        
        stop("not all users.guesses within bounds\n change values of `lower' and `upper' or \nthose of the `lowerbound*.factor's and `upperbound*.factor's")
      }
    }

    if (length(users.guess) > 0) {
      param.table[[M]][IDX("variab")] <- users.guess
      param.table[[M]][IDX("param")] <- trafo(users.guess)
      MLEVARIAB <- users.guess
      param.table[[M]][IDX("glbl.var")] <- get.var.covariat(users.guess)
    }
  }


##################################################
################### Empirical Variogram   ########################

 
##################################################
################### Empirical Variogram   ########################
 if (printlevel>=PL_STRUCTURE) cat("\nempirical variogram...")

  ## see above for the trafo definitions
  ##
  ## zuerst regression fit fuer variogram,
  ## dann schaetzung der Parameter, dann berechnung der covariates
  ## T.o.D.o.: kann verbessert werden durch einschluss der Trendschaetzung
  ## Xges auf RFsimulate basieren

  
  lsqMethods <- NULL
  ev <- list()
  if (length(lsq.methods) == 0) {  # not trans.inv
    message("nNo submethod is not allowed.")
    if (trace != 0) message("'trace' only works if a submethod is allowed")
    trace <- -abs(trace)
  } else { # trans.inv
    sd <- vector("list", vdim)
    index.bv <- NULL

    ## to do: non-translationinvariant case!!!
    
    if (vdim == 1 && !dist.given) {
      residuals <- .Call(C_simple_residuals, LiliReg) 
    } else {
      residuals <- list()
      for (i in 1:sets) {
        residuals[[i]] <- Z$data[[i]]
        base::dim(residuals[[i]]) <-
          c(coords[[i]]$totpts, vdim, repet[i])
      }
    }

    bins <- RFopt$empvario$bins
    if (length(bins) == 0) bins <- Default.nbins
    if (length(bins) == 1 && bins > 1) {
      bins <-c(-1, if (ts.xdim==1 && Z$coords[[1]]$grid && length(bins)==1)
             seq(Z$coords[[1]]$x[2] / 2, bins * Z$coords[[1]]$x[2], len=bins+1)
                   else seq(0, fit$bin_dist_factor * maxdistance, len=bins+1))
    }
    if ((nphi <- RFopt$empvario$nphi) <= 1 && spatialdim >= 2 && !isotropic)
      nphi <- Default.nangles
    if ((ntheta<-RFopt$empvario$ntheta) <= 1 && spatialdim>=3 && !isotropic)
      ntheta <- Default.nangles
    if (any((deltaT <- RFopt$empvario$deltaT) == 0) &&
        length(coords[[1]]$T) > 0) {
      deltaT <- as.integer(min(coords[[1]]$T[3] / 3, 100))
    }
  
    if (!dist.given) { ## todo
      ev <- rfempirical(x=coords, data=residuals, bins=bins,
                        nphi=nphi,
                        ntheta=ntheta,
			deltaT=deltaT, spConform=FALSE,
                        alpha=emp_alpha,
			vdim=vdim)
      n.bin <- ev$n.bin
      sd <- ev$sd
      binned.variogram <- ev$empirical
    } else {
      stop("das ist doch auch in C programmiert?!")
      
      binned.variogram <- NULL ## ja nicht loeschen !
      len <- sapply(Z$coords, function(x) x$totpts)      
      for (j in 1:vdim) {
        if (Z$xdimOZ != 1) stop("Distance vectors are not allowed.")
        n.bin <- vario2 <- vario <- rep(0, length(bins))
        for (i in 1:sets) {
          dimW <- dim(residuals[[i]])
          W <- residuals[[i]][, j, ] 
          dim(W) <- dimW[c(1,3)] ## necesary! as W could have dropped dimension
          lc <- len[i]
          rep <- 2 * repet[i]
          k <- 1
          for (g in 1:(lc-1)) {
            for (h in (g+1):lc) {
              idx <- sum(coords[[i]]$x[k] > bins) ## distances
              dW <- W[g, ] - W[h, ]
              n.bin[idx] <- n.bin[idx] + sum(!is.na(dW))
              vario[idx] <- vario[idx] + sum(dW^2, na.rm=TRUE)
              vario2[idx] <- vario2[idx] + sum(dW^4, na.rm=TRUE)
              k <- k + 1
            }
          }
        }
        n.bin <- n.bin * 2
        vario <- vario / n.bin ## Achtung! n.bin ist bereits gedoppelt        
        sdvario <-
          sqrt(pmax(0, vario2 / n.bin / 2.0 - vario^2)) ## numerische Fehler
        sdvario[1] <- vario[1] <- 0
        n.bin[1] <- sum(len)
        centers <- 0.5 * (bins[-1] + bins[-length(bins)])
        centers[1] <- 0
        empirical <- vario[-length(n.bin)]
        sdv <- sdvario[-length(n.bin)]
        tmp <- n.bin[-length(n.bin)]
        dims <-  c(length(empirical), rep(1, 5))
        dim(empirical) <- dim(sdv) <- dim(tmp) <- dims
        ev <- list(centers=centers, empirical=empirical, sd=sdv, n.bin=tmp,
                  stop(" Fehtl* noch * viel") ## auch alpha fehlt hier
                   )
	sd[[j]] <- ev$sd # j:vdim; ev contains the sets
	if (is.null(binned.variogram))
	  binned.variogram <- array(dim=c(length(ev$empirical), vdim, vdim))
	binned.variogram[, j, j] <- ev$empirical
      } # j in 1:vdim
    } ## dist.given

    ##   Print(ev)
    ## print(binned.variogram)

    if (!(sum(binned.variogram, na.rm=TRUE) > 0))
      stop("all value of the empirical variogram are NA; check values of bins and bin_dist_factor")
    

    if (!any(abs(binned.variogram) != 0))
      stop("all value of the empirical variogram are NA; check values of bins and bin_dist_factor")
   
    bin.centers <- as.matrix(ev$centers)

    if (length(ev$phi) > 0 || length(ev$theta)>0) {
      if (length(ev$phi) > 0) {
        if (spatialdim<2) stop("x dimension is less than two, but phi is given")
        bin.centers <- cbind(as.vector(outer(bin.centers, cos(ev$phi))),
                             as.vector(outer(bin.centers, sin(ev$phi))))
      }
      if (length(ev$theta) > 0) {
        if (spatialdim<3)
          stop("x dimension is less than three, but theta is given") 
        if (ncol(bin.centers)==1) bin.centers <- cbind(bin.centers, 0)
        bin.centers <- cbind(as.vector(outer(bin.centers[, 1],
                                             cos(ev$theta))),
                             as.vector(outer(bin.centers[, 2],
                                             cos(ev$theta))),
                             rep(sin(ev$theta), each=nrow(bin.centers)))
      }
     } else {
      
      ##  warning("must be ncol()")      
      if (ncol(bin.centers) < spatialdim) { # dimension of bincenter vector
        ##                       smaller than dimension of location space      
        bin.centers <- 
          cbind(bin.centers, matrix(0, nrow=nrow(bin.centers),
                                    ncol=spatialdim - ncol(bin.centers)
                                    ))
      }
    }


    ## es muessen beim direkten C-aufruf die componenten der Punkte
    ## hintereinander kommen (siehe auch variable coords, Xdistance). Deshalb t()
    ## aus den vdim sd's muss ein Gewicht gemacht werden
    
    if (length(sd) > 0 &&
	length(sd[[1]])>0 ## not satisfied if variogram is estimated by fft
	) { 
       if (is.list(sd)) {
        if (length(sd[[1]]) > 0)  {
          evsd <- sapply(sd, function(x) x^2) ## ALWAYS vector or matrix        
          if (is.matrix(evsd)) evsd <- rowMeans(evsd, na.rm=TRUE)
        } else evsd <- 1
      } else if (is.matrix(sd)) {
        evsd <- as.vector(apply(sd, 1, function(x) diag(x)))
      } else evsd <- sd
      evsd[!is.finite(evsd)] <- 10 * sum(evsd, na.rm=TRUE) ## == "infinity"
      evsd <- as.double(evsd)
    } else evsd <- 1

    L.n.bin             <- length(n.bin)
    binned.n         <- as.integer(n.bin)
    inv.evsd <- 1.0 / evsd
    inv.evsd[evsd == 0] <- 0

    weights <- cbind(self=NA,                      # self
                     plain=rep(1, L.n.bin),            # plain 
                     sqrt.nr=sqrt(binned.n),          # sqrt(#)
                     sd.inv=inv.evsd,                # sd^-1
                     internal1 = sqrt(L.n.bin:1 * as.double(binned.n)), # internal1 # kann sonst
                     ##                   fehler verursachen, da integer overflow
                     internal2 = L.n.bin:1,                  # internal2
                     internal3 = sqrt(L.n.bin:1)             # internal3
                     )
    stopifnot(ncol(weights)==length(LSQMETHODS))
    weights <- data.frame(weights)

    ##################################################
    #######################  LSQ  ####################
    ##***********   elimination part itself   **********     
    ## find a good initial value for MLE using weighted least squares
    ## and binned variogram
    ##
    ## background: if the number of observations (and the observation
    ## field) tends to infinity then any least square algorithm should
    ## yield the same result as MLE
    ## so the hope is that for a finite number of points the least squares
    ## find an acceptable initial values

    ## advantage of the following way is that the for-loop is run through
    ## in an ordered sense -- this might be useful in case partial results
    ## are reused
    if (printlevel>=PL_STRUCTURE) cat("\nelimination part...")

    lsqMethods <- LSQMETHODS[pmatch(lsq.methods, LSQMETHODS)]

 #   Print(Z, emp_alphaName, C_coords)
    
    .Call(C_SetAndGetModelLikelihood, COVreg, list(emp_alphaName, Z$model),
          C_coords, neverGlobalXT, MLE_CONFORM)

  #  Print(Z, C_coords, C_UnifyXT(x = bin.centers, T = ev$T))
    .Call(C_LocNonGrid, COVreg, C_UnifyXT(x = bin.centers, T = ev$T))
 
    firstoptim <- TRUE
    LSMIN <- Inf
    if (length(lsqMethods) > 0 &&
        any(is.na(lsqMethods))) stop("not all lsq.methods could be matched")

    if (trace > 0) {
      txt <- c("phi", "theta", "dT")
      d <- dim(binned.variogram) # c(n.bins, n.phi, n.theta, n.delta, vdim,vdim)
      if (length(d) != 6 || d[6] != vdim)
        stop("Programming error in 'Trace'.", CONTACT)
      xdim <- d[1]
      d <- d[-c(1,5,6)]
       if (vdim == 1) {
        w <- which(d != 1)
        if (length(w) == 3) {## all set
          idx <- 1:(1 + (d[1] <= d[3]))
          idx2 <- (1:3)[-idx]
          d <- c(prod(d[idx]), prod(d[idx2]))
          txt <- c(paste0(txt[idx], collapse="/"),
                   paste0(txt[idx2], collapse="/"))
        } else if (length(w)==2) {
          d <- d[w]
          txt <- txt[w]
        } else if (length(w) == 1){
          ## last one could be still pretty large. Make a rectangle out of it
          nc <- as.integer(sqrt(d[w]))
          d <- c(nc, d[w] / nc)
          txt <- paste0(txt[w], "[", c("1:", ". %% "), nc, "]")
        } else {
          d <- c(1,1)
        }
      } else {
        if (prod(d)== 1) { ## split vdim
          d <- rep(vdim, 2)
          txt <- rep("vdim", 2)
        } else { ## split between vdim and rest
          txt<- paste(txt[d > 1], collapse="/")
          d <- c(prod(d), vdim^2)
        }
      }
      Trace.grdim <- pmin(d, Trace.plots)
      Trace.dim <- c(xdim, d)
      Trace.txt <- c("distance", txt)      
      Trace.ev <- c(binned.variogram)
      Trace.model <- list()
      Trace.opt <- list()

      on.exit(par(par.old), add=TRUE)
      par(mar=Trace.mar, mfcol=Trace.grdim, cex=Trace.cex)   
    }

    for (M in lsqMethods) {
     if (printlevel>=PL_STRUCTURE) cat("\n", M) else cat(pch)

      param.table[[M]][IDX("variab")] <- default.param
      
      LSQsettings(M)
     
      LSMIN <- Inf ## must be before next "if (n.variab==0)"
      LSPARAM <- LSVARIAB <- NA 
      ##      if (n.variab == 0) {
                                        #       warning("trivial case may cause problems")
 #     } else {
      param.table[[M]][IDX("lower")] <- LSQLB
      param.table[[M]][IDX("upper")] <- LSQUB
      options(show.error.messages = show.error.message) ##

      if (n.variab == 0) {
        LStarget(param.table[IDX("variab"), methodprevto$lsq[1]])
      } else {
        min <- Inf
        min.variab <- NULL
        for (i in methodprevto$lsq) { ## ! -- the parts that change if
          ##                               this part is copied for other method
          if (!any(is.na(variab <- param.table[IDX("variab"), i]))) {
	    value <- LStarget(variab) ##
	    if (is.finite(value)) {
              param.table[tblidx[[M]][1], i] <- value
              if (value < min) {
                min.variab <- variab
                min <- value
              } else {
                param.table[tblidx[[M]][1], i] <- NaN
                next
              }
            }
          }
        }

        stopifnot(min==LSMIN) ## check
        if (any(min.variab != LSVARIAB) && printlevel > PL_SUBIMPORTANT) {
          Print("optimal variables, directly returned and by optim, differ",#
                min.variab, LSVARIAB)
        }
          
        
        fnscale <- if (length(fit.fnscale) == 0 || is.na(fit.fnscale[M]))
          min else fit.fnscale[M]

        lsq.optim.control <-
          c(opt.control, list(parscale=parscale, fnscale=fnscale))

        ##  lsq.optim.control <- c(opt.control, list(parscale=c(1,1), fnscale=-1))
        ##        Print(M, LSVARIAB, lower, upper, autostart, info.cov, parscale, lsq.optim.control,optimiser);
#        print(info.cov)

        SUBOPTIM(LSVARIAB, LStarget, lower = LSQLB, upper = LSQUB,
              control=lsq.optim.control)

      } # n.variab > 0
      options(show.error.messages = show.error.message)  
      ## side effect: minimum so far is in LSMIN and LSPARAM
      ## even if the algorithm finally fails


      if (is.finite(LSMIN)) {
        param.table[[M]][tblidx[[M]][1]] <- LSMIN
        param.table[[M]][IDX("variab")] <- LSVARIAB
        param.table[[M]][IDX("param")] <- LSPARAM

        ps <- abs(LSVARIAB)
        zero <- ps == 0
        parscale[!zero] <- ps[!zero]
        
      } else {
        param.table[[M]] <- if (n.variab==0) NA else NaN
      }

      param.table[[M]][IDX("glbl.var")] <- get.var.covariat(LSVARIAB)

      firstoptim <- FALSE
    } # for M   
  } # trans.inv


  if (length(LSQMETHODS) > 0) {
    ps <- matrix(NA, nrow=n.variab, ncol=length(LSQMETHODS))
    for (iM in 1:length(LSQMETHODS)) {
      M <- LSQMETHODS[iM]

      if (!is.na(param.table[[M]][tblidx[[M]][1]])) {
        ps[ , iM] <- abs(param.table[[M]][IDX("variab")])
      }
    }
    ps <- apply(ps, 1, median, na.rm=TRUE)
    parscale <- ParScale(optim.control, ps, lower, upper)
  }


##################################################
### optional parameter grid for MLE and CROSS  ###
  if (printlevel>=PL_STRUCTURE) cat("\nmle param...")
  
  
  idx <- IDX("variab")
  glblvar.idx <- IDX("glbl.var")

  
  gridmax <- as.matrix(param.table[idx, cm$lsq])
  if (!any(is.finite(gridmax))) gridmax <- param.table[idx, , drop=FALSE]
   
  gridmin <- apply(gridmax, 1, min, na.rm=TRUE)
  gridmax <- apply(gridmax, 1, max, na.rm=TRUE)

  gridbound <- lower
  gridbound[!is.finite(gridbound)] <- NA
  idx <- !is.na(gridbound)
  abase <- 0.25
  a <- is.na(gridmin[idx]) * (1-abase) + abase
  ## maybe there have not been any lsq elimination; then a=1
  gridmin[idx] <- (1-a) * gridmin[idx] + a * gridbound[idx]
  stopifnot(all(gridmin >= lower))
  gridbound <- upper
  gridbound[!is.finite(gridbound)] <- NA
  idx <- !is.na(gridbound)
  a <- is.na(gridmax[idx]) * (1-abase) + abase
  gridmax[idx] <- (1-a) * gridmax[idx] + a * gridbound[idx]
  stopifnot(all(gridmax <= upper))

 

##################################################
###################   MLE    #####################


  if (printlevel>=PL_STRUCTURE) cat("\nMLE XXX...")
  mleMethods <- (if (length(mle.methods)==0) NULL else
                 MLMETHODS[pmatch(mle.methods, MLMETHODS)])

##    Print(mleMethods)
  
  if ("reml" %in% mleMethods && n.covariat == 0)
    mleMethods <- c("ml")# to do, "reml", "rml")

  ## lowerbound_scale_ls_factor <  lowerbound_scale_factor, usually
  ## LS optimisation should not run to a boundary (what often happens
  ## for the scale) since a boundary value is usually a bad initial
  ## value for MLE (heuristic statement). Therefore a small
  ## lowerbound_scale_ls_factor is used for LS optimisation.
  ## For MLE elimination we should include the true value of the scale;
  ## so the bounds must be larger. Here lower[SCALE] is corrected
  ## to be suitable for MLE elimination
  
  if (any(MLELB > MLEUB))
    stop("the users lower and upper bounds are too restricitve")

  ## fnscale <- -1 : maximisation
  for (M in mleMethods) {
    if (printlevel>=PL_STRUCTURE) cat("\n", M) else cat(pch)
    param.table[[M]][IDX("variab")] <- default.param    
    
    if (M!="ml" ##&& !anyFixedEffect
        ) { ## also reml nicht nochmal rechnen...
      param.table[[M]] <- param.table[["ml"]]
      ## ML-value is now REML value:
      param.table[[M]][tblidx[[M]][1]] <- param.table[[M]][tblidx[["ml"]][1]]
      next
    }
     
    MLEMAX <- -Inf ## must be before next "if (nMLEINDEX==0)"
    MLEVARIAB <- NULL ## nachfolgende for-schleife setzt MLEVARIAB
     MLEPARAM <- NA
    onborderline <- FALSE
    if (n.variab == 0) MLtarget(NULL)
    else {
      param.table[[M]][IDX("lower")] <- MLELB
      param.table[[M]][IDX("upper")] <- MLEUB
      options(show.error.messages = show.error.message) ##
      max <- -Inf

      for (i in methodprevto$mle) { ## ! -- the parts that change if
        ##                             this part is copied for other methods
        ## should mle be included when M=reml?
        ## same for lsq methods as well: should previous result be included?

        if (!any(is.na(variab <- param.table[IDX("variab"), i]))) {

          value <- MLtarget(variab) ## !
#          Print(value)

          if (is.finite(value)) {
            param.table[tblidx[[M]][1], i] <- value
            if (value > max) {
              max.variab <- variab
              max <- value
            }
          } else {
            param.table[tblidx[[M]][1], i] <- NaN
            next
          }
        }
        
      }

     
      fnscale <-
        if (length(fit.fnscale)==0 || is.na(fit.fnscale[M]))
          -max(abs(max), 0.1) else fit.fnscale[M]

      mle.optim.control <-
        c(opt.control, list(parscale=parscale, fnscale=fnscale))

                                        ##
#      Print(parscale, MLEVARIAB, info.cov, methodprevto$mle, globalvariance)
#      print(minmax)

      if (length(parscale) > 0 && length(parscale) != length(MLEVARIAB))
        stop(#"length of 'parscale' (", length(parscale), ") differs from the length of the variables (", length(MLEVARIAB), "). ",
            if (length(MLEVARIAB) != 0) stop(CONTACT)
            else "Likely, there is a problem with the model.")
             
      
      MLEINF <- FALSE

      if (fit$critical < 2) {
        OPTIM(MLEVARIAB, MLtarget, lower = MLELB, upper=MLEUB,
              control=mle.optim.control)

        if (ML_failures > 20 && printlevel >= PL_IMPORTANT)
          message("There are many failures when trying to evaluate the model. Is the model OK?")
        
         if (MLEINF) {
          if (printlevel>=PL_STRUCTURE)
            Print("MLEINF", MLEVARIAB, MLEMAX) else cat("#") #
          OPTIM(MLEVARIAB, MLtarget, lower = MLELB, upper=MLEUB,
                control=mle.optim.control)
          if (printlevel>=PL_STRUCTURE)
            Print("MLEINF new", MLEVARIAB, MLEMAX) #
        }
              
        options(show.error.messages = TRUE) ##
        mindist <-
          pmax(fit$minbounddistance, fit$minboundreldist * abs(MLEVARIAB))
        
        onborderline <- 
          (abs(MLEVARIAB - MLELB) <
           pmax(mindist,  ## absolute difference
                fit$minboundreldist * abs(MLELB) ## relative difference
                )) |
        (abs(MLEVARIAB - MLEUB) <
         pmax(mindist, fit$minboundreldist * abs(MLEUB)))
      }
    } # length(MLELB) > 0
  
    if (printlevel>=PL_STRUCTURE)
      Print("mle first round", MLEVARIAB, MLEPARAM, MLEMAX) #
    
    if (!is.finite(MLEMAX)) {
      if (printlevel>=PL_IMPORTANT) message(M, ": MLtarget I failed.")
      param.table[[M]] <- MLEPARAM <- NaN
      variab <- MLELB ## to call for onborderline
      ml.residuals <- NA
    } else {
      param.table[[M]][tblidx[[M]][1]] <- MLEMAX
      param.table[[M]][IDX("variab")] <- MLEVARIAB

      stopifnot(all(MLEVARIAB >= MLELB & MLEVARIAB <= MLEUB))
      
      param.table[[M]][IDX("param")] <- MLEPARAM
      param.table[[M]][IDX("glbl.var")] <- get.var.covariat(MLEVARIAB)     
      ml.residuals <- ML.RESIDUALS


      if (FALSE) Print(recall, fit$critical>=2 || any(onborderline),#
             !fit$only_users , fit$critical >= 0)
      
      if (fit$critical>=2  || (any(onborderline)
        && !fit$only_users && fit$critical >= 0)) {
        ## if the MLE result is close to the border, it usually means that
        ## the algorithm has failed, especially because of a bad starting
        ## value (least squares do not always give a good starting point,helas)
        ## so the brutal method:
        ## calculate the MLE values on a grid and start the optimization with
        ## the best grid point. Again, there is the believe that the
        ## least square give at least a hint what a good grid is
        
        if (fit$critical == 0) {
          MLEgridmin <- gridmin
          MLEgridmax <- gridmax
          
          if (any(is.na(MLEgridmin)) || any(is.na(MLEgridmax))) {
            if (printlevel >= PL_SUBIMPORTANT) {
              Print(cbind(MLELB, variab, MLEUB, onborderline), #
                    MLEgridmin, MLEgridmax)
            }
            warning(paste(M, "converged to a boundary value -- ",
                          "better performance might be obtained",
                          "when allowing for more lsq.methods"))
          } else {
            if (printlevel>=PL_FCTN_SUBDETAILS)
              show(1, M, MLEMAX, MLEVARIAB) else cat(detailpch)
            MLEgridlength <-
              max(3, round(fit$approximate_functioncalls^(1/n.variab)))
            ## grid is given by the extremes of the LS results
            ## so, therefore we should examine above at least 4 different sets
            ## of weights wichtig: gridmin/max basiert auf den reduzierten
            ## Variablen 
            step <- (MLEgridmax - MLEgridmin) / (MLEgridlength-2) # grid starts
                                        # bit outside
            MLEgridmin <- pmax(MLEgridmin - runif(length(MLEgridmin)) * step/2,
                               MLELB)   # the extremes of LS
            MLEgridmax <- pmin(MLEgridmax + runif(length(MLEgridmax)) * step/2,
                               MLEUB)
            step <- (MLEgridmax - MLEgridmin) / (MLEgridlength-1)
            
            startingvalues <- vector("list", length(step))
            for (i in 1:length(step)) {
              startingvalues[[i]] <-
                MLEgridmin[i] + step[i] * 0:(MLEgridlength-1)
            }
            
            startingvalues <- do.call(base::expand.grid, startingvalues)
            
            limit <- 10 * fit$approximate_functioncalls
            if ((rn <- nrow(startingvalues)) > limit) {              
              if (printlevel>=PL_STRUCTURE)
                cat("using only a random subset of the", rn, "grid points")
              rand <- runif(rn)
              startingvalues <-
                startingvalues[rand < quantile(rand, limit / rn), , drop=FALSE]
              gc()
            }
            
            ## MLEMAX <- -Inf
            
            apply(startingvalues, 1, ## side effect: MLEVARIAB
                  function(x) {
##                    TRY(MLtarget(x))
##                    stopifnot( MLEMAX == MLtarget(MLEVARIAB))
              ## try-result:
                  ##     if (!is.numeric(z) && abs(z)<1e10) cat(paste(x, sep=","), "\n") else cat(".\n")
                    ##  if (MLEINF) stop("stop")
                    TRY(MLtarget(x))
                })
            
            if (printlevel>=PL_FCTN_SUBDETAILS)
              Print("mle grid search", MLEVARIAB, MLEPARAM, MLEMAX) #
            
            ## side effect:Maximum is in MLEMAX!
            ##                             and optimal parameter is in MLEVARIAB
            if (printlevel>=PL_FCTN_SUBDETAILS)
              show(2, M, MLEMAX, MLEVARIAB)
            
            cat(detailpch)
            options(show.error.messages = show.error.message) ##

            OPTIM(MLEVARIAB, MLtarget, lower = MLELB, upper = MLEUB,
                  control=mle.optim.control)
            options(show.error.messages = TRUE) ##
            if (!is.finite(MLEMAX) &&(printlevel>=PL_IMPORTANT))
              message("MLtarget II failed.\n")
            ## do not check anymore whether there had been convergence or not.
            ## just take the best of the two strategies (initial value given by
            ## LS, initial value given by a grid), and be happy.
            if (printlevel>=PL_FCTN_SUBDETAILS) show(3, M, MLEMAX, MLEVARIAB)
	    else if (printlevel>=PL_STRUCTURE)
              Print("mle second round", MLEVARIAB, MLEPARAM, MLEMAX) #
            
            if (is.finite(MLEMAX) && MLEMAX >=param.table[[M]][tblidx[[M]][1]]){
              param.table[[M]][tblidx[[M]][1]] <- MLEMAX
              param.table[[M]][IDX("variab")] <- MLEVARIAB
              stopifnot(all(MLEVARIAB >= MLELB & MLEVARIAB <= MLEUB))
              param.table[[M]][IDX("param")] <- MLEPARAM
              param.table[[M]][IDX("glbl.var")] <- get.var.covariat(MLEVARIAB)
              ml.residuals <- ML.RESIDUALS
                     
            }
          } # (is.na(MLEgridmin[1]))
        } else { # fit$critical > 0  

          critical <- ptype == CRITICALPARAM;
          if (fit$critical>=3 || !any(critical)) critical <- rep(TRUE, n.variab)
          ncrit <- as.integer(fit$n_crit^(1/sum(critical)))
          if (!is.finite(ncrit)) ncrit <- 2
          if (ncrit^sum(critical) > 100 && printlevel>=PL_IMPORTANT)
            message("The optimisation may last a pretty long time!")
          
          if (ncrit > 1 || fit$critical>=2) {
            if (!is.null(transform)) {
              stop("if 'transform' is given, 'critical' must be '0'")
            }
            w <- 0.5
            lowlist <- upplist <- list()
            for (i in 1:n.variab) {
              if (critical[i]) {               
                newparam <- seq(if (SDVAR.IDX[i]) 0 else MLELB[i], MLEUB[i],
                                len=ncrit+1)
                lowlist[[i]] <- newparam[-length(newparam)]
                upplist[[i]] <- newparam[-1]
              } else {
                lowlist[[i]] <- MLELB[i]
                upplist[[i]] <- MLEUB[i]
              }
            }
            lowlist <- as.matrix(do.call(base::expand.grid, lowlist))
            upplist <- as.matrix(do.call(base::expand.grid, upplist))
            
            orig.MLEVARIAB <- MLEVARIAB
            b.idx <- is.finite(bounds)

            stopifnot(length(lowlist) > 1)
          
            for (i in 1:nrow(lowlist)) {
              cat(detailpch)
              new.lower.vector <- trafo(lowlist[i, ])
              new.upper.vector <- trafo(upplist[i, ])
              new.bounds <- TRUE
                        
              new.parscale <- guess <- fnscale <- NULL   
              first.passage <- TRUE
              while (TRUE) {
                if (new.bounds) {
                  new.bounds <- FALSE
                  .C(C_PutValuesAtNA, LiliReg, as.double(new.lower.vector))
                  new.lower <- GetModel(register=LiliReg,modus=GETMODEL_DEL_MLE,
                                        which.submodels = "user.but.once+jump",
					return.which.param=INCLUDENOTRETURN,
					origin = MLE_CONFORM)

                  .C(C_PutValuesAtNA, LiliReg, as.double(new.upper.vector))
                  new.upper <- GetModel(register=LiliReg,modus=GETMODEL_DEL_MLE,
                                        which.submodels="user.but.once+jump",
 					return.which.param=INCLUDENOTRETURN,
					origin = MLE_CONFORM)
                }
               
                if (printlevel>=PL_STRUCTURE) cat("\nrecall rffit...\n")
               
                res <-
                  rffit.gauss(Z= Z,
                              lower=new.lower,
                              upper=new.upper,
                              mle.methods="ml",
                              lsq.methods=lsq.methods,
                              users.guess = guess,
                              optim.control=c(opt.control,
                                fnscale=list(fnscale),
                                parscale=list(new.parscale)),
                              recall = TRUE,
                              general.pch = if (pch == "") "" else ":",
                              general.practicalrange = general$practicalrange,
                              general.spConform = FALSE,
                              fit.critical = -1,
                              fit.split =FALSE)
                
 
                guess <- res$table[[M]][IDX("param")]
                ## scale_ratio <- abs(log(abs(guess[-delete.idx]/new.parscale)))

                stopifnot(length(users.lower) # + length(delete.idx)
                          == length(new.lower.vector))
                
                if (!all(guess >= new.lower.vector & guess <= new.upper.vector)
                    || !all(new.lower.vector > users.lower &
                            new.upper.vector < users.upper)) {
                  new.lower.vector <- pmin(new.lower.vector, guess)
                  new.upper.vector <- pmax(new.upper.vector, guess)
                  new.lower.vector <- pmax(new.lower.vector, users.lower)
                  new.upper.vector <- pmin(new.upper.vector, users.upper)
                  guess <- pmax(new.lower.vector, pmin(new.upper.vector, guess))
                  new.bounds <- TRUE
                }

                
                likelihood <- res$table[[M]][tblidx[[M]][1]]

                if (is.finite(likelihood)) {
                  if (likelihood > param.table[[M]][tblidx[[M]][1]]) {
                    if (printlevel > PL_RECURSIVE && !first.passage)
                      cat("parscale: table improved by ",
                          likelihood -  param.table[[M]][tblidx[[M]][1]], "\n") 
                    param.table[[M]][tblidx[[M]][1]] <- MLEMAX <- likelihood
                    for (j in c("variab", "param", "glbl.var"))
                      param.table[[M]][IDX(j)] <- res$table[[M]][IDX(j)]
                    ml.residuals <- res$ml$residuals
                  }
                }
                
                if (first.passage) {
                  old.likelihood <- likelihood
                } else {
                  if (printlevel > PL_RECURSIVE) {
                    if (likelihood > old.likelihood) {                      
                      cat("parscale: mle improved by", likelihood-old.likelihood,
                          "\n")
                    }
                  }
                  break;
                }
                            
                ## urspruengliches Abbruchkriterium, das nicht gut fkt:
                ##value.ratio <- abs(log (abs(likelihood) / abs(MLEMAX)))
                ##outside <- (scale_ratio > fit$scale_ratio)  & MLEVARIAB != 0
                ##outside <- outside & (!b.idx | new.parscale >
                ##                      exp(fit$scale_ratio) * 
                ##                      (w * abs(guess) + (1-w) * bounds))
                ## if (!any(outside) && value.ratio <= fit$scale_ratio) {
               ##   break
               ## }

                
                zero <- new.parscale == 0
                new.parscale[zero] <-
                  pmax(abs(lower[zero]), abs(upper[zero])) / 10                
                new.parscale[b.idx] <-
                  w * new.parscale[b.idx] + (1-w) * bounds[b.idx]
                #new.parscale <- abs(guess[-delete.idx])

                restable <- as.matrix(res$table)
                used <- which(!apply(is.na(restable), 2, all)[-1:-2]) # ohne auto,user
                names.tbl <- names(tblidx)
                start <- which(names.tbl %in% LSQMETHODS)              
                start <- if (length(start) == 0) 0 else min(start) - 1
                if (start>0) names.tbl <- names.tbl[-1:-start]
                 
                fnscale <- numeric(length(used))                
                for (j in 1:length(used)) {
                  fnscale[j] <- abs(restable[tblidx[[start + used[j]]][1],
                                              2 + used[j]])
                  if (names.tbl[used[j]] == "ml")
                    fnscale[j] <- -max(0.1, fnscale[j])
                }

                .C(C_PutValuesAtNA, LiliReg, guess, NAOK=TRUE) # ok
                guess <- GetModel(register=LiliReg, modus=GETMODEL_DEL_MLE,
                                  which.submodels = "user.but.once+jump",
 				  return.which.param=INCLUDENOTRETURN,
				  origin = MLE_CONFORM)
                first.passage <- FALSE
              } # while true
            } # for 1:nrow(lowlist)            
          } # ncrit  > 1          
        } # fit$critical > 0
      }
      if (printlevel>=PL_STRUCTURE) cat("\nend MLE critical part ...\n")
 
      if (!fit$only_users && (fit$reoptimise || any(SDVAR.IDX)) &&
          fit$critical >= 0 && fit$critical <= 1){
        if (pch!="") cat("$")
        new.parscale <- abs(revariab <- param.table[[M]][IDX("variab")])
        new.parscale[new.parscale == 0] <- 1e-3
        fnscale <- -max(0.1, abs(MLEMAX <- param.table[[M]][tblidx[[M]][1]]))
      
        old.control<- if (exists("mle.optim.control")) mle.optim.control else NA
        mle.optim.control <- c(opt.control,
                               list(parscale=new.parscale, fnscale=fnscale))


        if (any(SDVAR.IDX))
          MLELB[SDVAR.IDX] <- pmax(MLELB[SDVAR.IDX]/100,
                                   users.lower[SDVAR.IDX])
       

        OPTIM(revariab, MLtarget, lower = MLELB, upper=MLEUB, control=mle.optim.control)

        if (MLEMAX >= param.table[[M]][tblidx[[M]][1]]) {   
          if (printlevel > PL_SUBIMPORTANT)
            Print(old.control, fnscale, revariab, #
                  param.table[[M]][IDX("param")],
                  MLEMAX, MLEVARIAB,  MLEPARAM, MLELB, MLEUB)
          param.table[[M]][tblidx[[M]][1]] <- MLEMAX
          param.table[[M]][IDX("variab")] <- MLEVARIAB         
          param.table[[M]][IDX("param")] <- MLEPARAM
          param.table[[M]][IDX("glbl.var")] <- get.var.covariat(MLEVARIAB)
        }
      }
    } # is.finite(MLEMAX)     
  } ## M(le)



######################################################################
###                     error checks                               ###
######################################################################
  ## if rather close to nugget and nugget==fixed, do exception handling.
  ## This part is active only if
  ## scale_max_relative.factor < lowerbound_scale_ls_factor
  ## By default it is not, i.e., the following if-condition
  ## will "always" be FALSE.
  if (printlevel>=PL_STRUCTURE) cat("\nerror checks ...\n")
  trace <- 0
  
  if (globalvariance && length(nugget.idx) == 0 && any(SCALE.IDX)) {
    idx <- IDX("variab")
    alllsqscales <- param.table[idx, cm$lsq][SCALE.IDX, ]
   
    if (any(alllsqscales < mindistance/fit$scale_max_relative_factor,
            na.rm=TRUE))
      warning(paste(sep="", "Chosen model seems to be inappropriate!\n \
                    Probably a (larger) nugget effect should be considered"))
  }


######################################################################
###                   prepare lower and upper                      ###
######################################################################

  ## if the covariance functions use natural scaling, just
  ## correct the final output by GNS$natscale
  ## (GNS$natscale==1 if no rescaling was used)
  ##

  lower <- trafo(lower)
  upper <- trafo(upper)
  idx <- lower == upper
  lower[idx] = upper[idx] = NA

  
  .C(C_PutValuesAtNA, LiliReg, as.double(lower), NAOK=TRUE) # ok
  lower <- GetModel(register=LiliReg, modus=GETMODEL_SOLVE_MLE,
                    which.submodels = "user.but.once+jump",
                    C_conform = FALSE, origin = MLE_CONFORM)
  .C(C_PutValuesAtNA, LiliReg, as.double(upper), NAOK=TRUE) # ok
  upper <- GetModel(register=LiliReg, modus=GETMODEL_SOLVE_MLE,
                    which.submodels = "user.but.once+jump",
		    C_conform = FALSE, origin = MLE_CONFORM)

######################################################################
###                   prepare models for returning                 ###
######################################################################

  if (printlevel>=PL_STRUCTURE) cat("preparing for returning ...\n");

  idxCovar <- IDX("covariat")
  idxpar <- IDX("param") # henceforth
  if (globalvariance) idxvar <- IDX("glbl.var")
  res <- values.res <- list()
  totalparam <- as.integer(n.variab + n.covariat) ## not n.param
  AICconst <-
    2 * totalparam + 2 * totalparam * (totalparam + 1) / (sum.not.isna.data - totalparam - 1) 

  allMethods <- c(primMethods, lsqMethods, mleMethods)

  DataNames <- Z$DataNames
  Hessians <- invH <- list()

  for (i in OneTo(length(allmethods))) {
    ##    Print(i)
    Meth_i <- allmethods[i]
    if (!(Meth_i %in% allMethods)) next
    
    if (is.na(param.table[[Meth_i]][1])## && !is.nan(param.table[1,i]) ?? unklar
        && (length(idxpar)!=1 || idxpar != 0)) next ## result with NAs
    
    ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ##     calculate all target values for all optimal parameters     +++
    ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
    if (printlevel>= PL_STRUCTURE) cat("calculating method", i, "... ")
    v <- param.table[[Meth_i]][IDX("variab")] 
    for (M in LSQMETHODS) {
      if (M %in% lsqMethods) {          
        LSQsettings(M)
        param.table[[Meth_i]][tblidx[[M]][1]] <- LStarget(v)
      }
    } 

    for (M in MLMETHODS) {
      cur <- param.table[[Meth_i]][tblidx[[M]][1]]
      if (is.na(cur[1]) && !is.nan(cur[1]) && M %in% mleMethods) {
        param.table[[Meth_i]][tblidx[[M]][1]] <- MLtarget(v)
      }

      if (Meth_i== M) {
        H <- INVDIAGHESS(param.table[[M]][IDX("variab")], MLELB=MLELB,
                      MLEUB=MLEUB, control=mle.optim.control)
        param.table[[M]][IDX("sdvariab")] <- H$sd
        Hessians[[Meth_i]] <- H$hessian
        invH[[Meth_i]] <- invH
      }
    }


    for (M in CROSSMETHODS) {
      cur <- param.table[[Meth_i]][tblidx[[M]][1]]
      if (is.na(cur) && !is.nan(cur) && M %in% NULL) { # crossMethods) {
        stop("not programmed ")
          ##  crosssettings(M) ## uncomment
        variab <- v
        if (n.covariat > 0) {
          variab <- c(variab, param.table[[Meth_i]][idxCovar])
        }
        param.table[[Meth_i]][tblidx[[M]][1]] <- stop("") # crosstarget(variab)
      }
    }

    ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ## die nachfolgenden Zeilen aendern die Tabellenwerte -- somit
    ## muessen diese Zeilen unmittelbar vor dem return stehen !!!
    ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if (fit$use_naturalscaling && any(SCALE.IDX)) {
      param <- as.double(param.table[[Meth_i]][idxpar]  + 0.0) ## + 0.0 MUSS STEHEN
      ## register muss mit cov initialisiert sein, sonst
      ## wird laueft er in cov->key rein!
      .C(C_PutValuesAtNA, LiliReg, param)
      .C(C_expliciteDollarMLE, LiliReg, param)
      param.table[[Meth_i]][idxpar]  <- param
    }


    ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ##  models with NAs filled
    ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if (printlevel>= PL_STRUCTURE)cat("with nas filled ... ")

    .C(C_PutValuesAtNA, LiliReg, as.double(param.table[[Meth_i]][idxpar] ))


    if (globalvariance) .C(C_PutGlblVar, as.integer(LiliReg),
                           as.double(param.table[[Meth_i]][glblvar.idx]))

    modelres <- GetModel(register = LiliReg, modus = GETMODEL_SOLVE_MLE,
                         C_conform = FALSE,
			solve_random = TRUE,
			which.submodels = "user.but.once+jump",
			origin = original_model)

    ## print("modelres")
    ## print(modelres)
     
    if (!is.null(DataNames)) {
      DataNames <- if (length(DataNames) == 0) ""
      else paste("c(", paste(DataNames, collapse=","), ")")
      txt <- paste(DataNames, "~", model2string(modelres))
      formel <- as.formula(txt)
    } else formel <- NULL
    
    ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ##   AIC etc for MLE
    ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if (printlevel>= PL_STRUCTURE) cat("AIC...\n")
    M <- allmethods[i]
 #   if (M == "ml")  {
      ## jetzt muss RFlikelihood-initialisierung sein!

    
    if ("ml" %in% mleMethods) {
      likelihood <- param.table[[Meth_i]][tblidx[["ml"]][1]]
      lilihood <- TRY(.Call(C_EvaluateModel, double(0), integer(0), LiliReg))
      if (!is(lilihood, CLASS_TRYERROR)) {
        residu <- get.residuals(LiliReg)       
        large.diff <- abs(lilihood[1] - likelihood) > 1e-7 * (abs(likelihood)+1)
        if (is.na(large.diff))
          stop("The likelihood could not be calculated",
               maintainer.info("for '", allmethods[i], "'"),
               ". Maybe the model definition is incorrect or the parameters ",
               "have no or incorrect bounds.")
        if (large.diff) {
          ##stop("The two ways of calculating the likelihood do not match: ",
          ##      lilihood[1], " != ", likelihood)
        }
      } else residu <- lapply(Z$data, function(x) x * NA)
      AIC <- 2 * totalparam  - 2 * likelihood
      AICc <- AICconst - 2 * likelihood
      BIC <- log(sum.not.isna.data) * totalparam - 2 * likelihood      
      param.table[[Meth_i]][tblidx[["AIC"]][1]] <- AIC
      param.table[[Meth_i]][tblidx[["AICc"]][1]] <- AICc
      param.table[[Meth_i]][tblidx[["BIC"]][1]] <- BIC
    } else {
      likelihood <- AIC <- AICc <- BIC <- NA_real_
      residu <- lapply(Z$data, function(x) x * NA)
    }
    
    
    if (n.covariat > 0) {
      c.table <- rbind(param.table[[Meth_i]][IDX("covariat")], NA)
      dimnames(c.table) <- list(NULL, betanames) 
    } else c.table <- NULL
    v.table <- rbind(param.table[[Meth_i]][IDX("variab"), drop=FALSE],
                     param.table[[Meth_i]][IDX("sdvariab"), drop=FALSE],
                     param.table[[Meth_i]][IDX("lower"), drop=FALSE],
                     param.table[[Meth_i]][IDX("upper"), drop=FALSE])

     
    dimnames(v.table) <- list(c("value", "sd", "lower", "upper"), variabnames)
    p.table <- if (n.variab > 0) rbind(param.table[[Meth_i]][IDX("param")], NA)
               else matrix(ncol=0, nrow=2)
    dimnames(p.table) <- list(c("value", "sd"), minmax.names)
    
    
    res[[M]] <-
      list(model=modelres,
           formel = formel,
           variab = v.table,
           coordsystem = coordsystem,
           param = p.table,
           covariat = c.table,
           globalvariance=
             if (globalvariance) param.table[[Meth_i]][IDX("glbl.var")][1],
           hessian = NULL,
           likelihood = likelihood,
           AIC = AIC, AICc = AICc, BIC = BIC,
           residuals = residu,
           params.list = if (is.null(params.fctn)) list() else params.fctn(v)
           )


    ##Print(res[[M]])
    
    class(res[[M]]) <- CLASS_SINGLEFIT
  }

  if (pch!="" && !recall) cat("\n")
                     
  ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ## Hess-Matrix for the parameters (for the variables see above)
  ## and further error checking
  ##
  ## next lines must be the very last call as standards are overwritten
  ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
   if (printlevel>= PL_STRUCTURE) cat("hessian ...\n")

  ## prepare MLtarget for getting 'param' directly
  if (n.param > 0) {
    if (is.null(transform)) {
      for (i in OneTo(length(mleMethods))) {
        M <- mleMethods[i]
        res[[M]]$param[2, 1:n.param] <- param.table[[M]][IDX("sdparam")] <-
           param.table[[M]][IDX("sdvariab")]
        res[[M]]$hessian <- Hessians[[M]]
      }
    } else {
      for (i in OneTo(length(mleMethods))) {
        M <- mleMethods[i]
        variab <- param.table[[M]][IDX("variab")] ## before orig.trafo
        param <- param.table[[M]][IDX("param")] ##  after orig.trafo
        eps.max <- 1e-5
        eps <- pmin(pmax(variab - MLELB, MLEUB - variab), eps.max)        
        eps <- diag(nrow=n.variab, eps * (2 * (variab <= MLEUB - eps) - 1))
        dtrafo <- apply(eps, 1, ## zeile: param, spalte: variab
                        function(h) (trafo(variab + h) - param) / sum(h))
        s <- svd(dtrafo, nv=Inf, nu=Inf)
        pseudoinv <- s$v %*% diag(ncol=n.param, nrow=n.variab, c(1/s$d)) %*%
          t(s$u)
     
        ## Print(dtrafo, n.param, n.variab, pseudoinv, pseudoinv %*% dtrafo, dtrafo %*% pseudoinv)
      
        ## TODO: not clear, whether next line makes sense of not
        res[[M]]$hessian <- t(pseudoinv) %*% Hessians[[M]] %*% pseudoinv
        ## Print(res[[M]]$hessian)
        if (length(invH[[M]]) == 0) res[[M]]$param[2, 1:n.param] <- NA
        else {
          res[[M]]$param[2, 1:n.param] <-
            sqrt(-diag(dtrafo %*% invH[[M]] %*% t(dtrafo)))
        }
      } #for methods 
    } ## not identity

    ## at the very end a check:
    L <- trafo(MLELB)
    U <- trafo(MLEUB)
    trafo <- function(x) x
    for (i in OneTo(length(mleMethods))) {
      M <- mleMethods[i]
      param <- param.table[[M]][IDX("param")] ##  after orig.trafo
      MLELB <- pmin(L, U, param, na.rm=TRUE)
      MLEUB <- pmax(L, U, param, na.rm=TRUE)
      cur <- param.table[[M]][tblidx[["ml"]][1]]
      
      if (!is.na(cur)) { ## kann bei autostart + Bayes auftreten
        mle <- MLtarget(param)
        difference <- abs(cur - mle)
        if ( (difference > 1e-10 && printlevel > PL_IMPORTANT) ||
             (difference > 1e-4 && printlevel >= PL_IMPORTANT) )
          message("likelihoods are calculated in two different ways leading",
                  " to the values ", cur, " and ", mle,
                  ", which are not the same.", CONTACT)
      }
    }
  
    
  } ## n.param > 0
        
   
  if (printlevel >= PL_IMPORTANT && BEYOND > 100)
    cat("\nNote: There are ",
        if (BEYOND > 1000)"very strong" else
        if (BEYOND > 250) "strong" else "some",
        " indications that the model might be overparametrised",
        "\nor that the bounds for the variables are too wide. ",
        "Try narrower lower and\nupper bounds for the variables in the ",
        "latter case. One of the critical\totalparameters is ",
        "'lowerbound_var_factor' whose value (", fit$lowerbound_var_factor,
        ") might be reduced.\n", sep="")
  
  if (printlevel>= PL_STRUCTURE) cat("S3 / S4 ...\n")

  L <- list(ev = ev,
            table=param.table,
            n.variab = as.integer(n.variab),
            n.param = as.integer(n.param),
            n.covariates = as.integer(n.covariat),
            lowerbounds=lower,
            upperbounds=upper,            
            number.of.data= as.integer(sum.not.isna.data),
            number.of.parameters = totalparam,
            modelinfo = nice_modelinfo(minmax),               
            p.proj = integer(0),
            v.proj = integer(0),
            x.proj = TRUE,
            fixed = NULL,
            true.tsdim = as.integer(tsdim),
            true.vdim = as.integer(vdim),
            trafo = trafo,
            report = "",
            submodels = NULL)

  if (!general$spConform) {
    V <- "'$vario' is defunctioned. Use '$ml' instead of '$value$ml'!"
    Res <- c(L,
	     list(vario = V), ## must be the very last of the parameter list!
             res)
    class(Res) <- CLASS_FITLIST
    return(Res)
  }
  
  res2 <- vector("list", length(res))
  nm <- names(res)

 ## str(res, max=1)
  for (i in 1:length(res)) { ## laeuft ueber Methoden, nicht ueber data.sets!
    ##    Print(i, res[[i]], coords, coords[[i]])
    res2[[i]] <- list2RMmodelFit(res[[i]], RFsp.info=Z$RFsp.info,
                                 T=coords[[1]]$T) ## nicht i !!
    attr(res2, "type") <- "gauss"
    ## wie uebergebe ich hier multiple Datenssaetze???
  }
  names(res2) <- nm

  return(rffit.gauss2sp(res2, L, Z))
}

rffit.gauss2sp <- function(res2, L, Z) {
  if (length(L$ev) > 0) L$ev <- do.call("new", c(list(CLASS_EMPIR), L$ev))
  L$lowerbounds <- list2RMmodel(L$lower)
  L$upperbounds <- list2RMmodel(L$upper)
  do.call.par <- c(list(Class = CLASS_FITLIST,
                        Z=Z,
                        coordunits=Z$coordunits,
                        varunits= if (length(L$ev) > 0) L$ev@varunits
                        else ""
                       ),
                   L,
                   res2)			
##  str(do.call.par, max.level=2)
  do.call(methods::new, args=do.call.par)
}

