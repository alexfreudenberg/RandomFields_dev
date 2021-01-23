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

# source("rf.R")
# see getNset.R for the variable .methods


prepare4RFspDataFrame <- function(x=NULL, RFopt, coordnames=NULL) {
  if (length(coordnames)==0) coordnames <- SystemCoordNames(x, opt=RFopt$coord)
  
  ## coords or GridTopology 
  if (x$grid) {
    coords <- NULL
    xgr <- cbind(if (length(x$xgr) == 0) x$x else x$gr, x$T)
##     Print(xgr, coordnames)
   colnames(xgr) <- coordnames
    xgr[is.na(xgr)] <- 0
    gridTopology <- sp::GridTopology(xgr[1, ], xgr[2, ], xgr[3, ])
  } else { ## grid == FALSE
    gridTopology <- NULL
    
    ## cbind of locations from x-matrix and T (if given)
    coords <- if (length(x$spatialtotpts) > 0) ## i.e. from KeyInfo.cc
                t(x$x) else x$x
    if (x$has.time.comp) {
      coords <- as.matrix(apply(coords, 2, rep, times=x$T[3]))
      coords <- cbind(coords, rep(seq(x$T[1], by=x$T[2], len=x$T[3]),
                                  each=x$totpts / x$T[3]))
    }
 ##   Print(coords, coordnames)
    if (is.matrix(coords)) colnames(coords) <- coordnames
  }
  Help("newstyle")
  return(list(coords=coords, gridTopology=gridTopology))
}


RFboxcox <- function(data, boxcox, vdim=1, inverse=FALSE, ignore.na=FALSE) {
  ## Fctn darf intern nicht aufgerufen werden!
  internalRFoptions(COPY=TRUE)
  boxcoxIntern(data = data, boxcox=boxcox, vdim=vdim, inverse=inverse,
                 ignore.na=ignore.na)
}
  
boxcoxIntern <- function(data, boxcox, vdim=1, inverse=FALSE, ignore.na=FALSE) {
  if (missing(boxcox)) boxcox <- .Call(C_get_boxcox, NULL)
  if (any(is.na(boxcox)) && !ignore.na)
    stop("non-finte values in Box-Cox transformation")
  if (!all(is.finite(boxcox))) return(data)
  if (is.list(data)) {
    for (i in 1:length(data))
      data[[i]] <- boxcoxIntern(data[[i]], boxcox=boxcox, vdim=vdim,
                                inverse=inverse, ignore.na=ignore.na)
    return(data)
  }
  Data <- data + 0
  .Call(C_BoxCox_trafo, as.double(boxcox), as.double(Data), as.integer(vdim),
        as.logical(inverse));
  return(Data)
}

RFlinearpart <- function(model, x, y = NULL, z = NULL, T=NULL, grid=NULL,
                         params,  distances, dim, set=0, ...) {
  Reg <- MODEL_USER  
  RFopt <- internalRFoptions(...)
  if (!hasArg("COPY")) on.exit(optionsDelete(RFopt, register=Reg))

  new <- UnifyXT(x, y, z, T, grid=grid, distances=distances, dim=dim)

  if (isS4(model) || !is(model, CLASS_CLIST)) {
    pm <- PrepareModel2(if (hasArg("trend")) list(...)$trend else model,
                        params=params, x=new, ...)
    new <- pm$C_coords
    model <- pm$model
  }
  if (model[[1]] != "linearpart") model <- list("linearpart", model)
 
  rfInit(model=model, x=new, reg=Reg, RFopt=RFopt)
  .Call(C_get_linearpart, Reg, as.integer(set))
}



setvector <- function(model, preceding, len, factor) {
  if (model[[1]] == SYMBOL_PLUS) {
    return(c(SYMBOL_PLUS,
             lapply(model[-1], setvector, preceding=preceding, len=len,
                    factor=factor)))
  }
  if (found <- model[[1]] == R_C) {
    if (length(model) != len + 1) stop("bug")
    model <- c(model[1], rep(0, preceding), model[-1])
    if (!missing(factor)) model <- list(SYMBOL_MULT, model, factor)
  } else if (found <- model[[1]] == R_CONST) {
    if (length(model[[1]]) != len) stop("bug")
    model[[2]] <- c(rep(0, preceding), model[[2]])
    if (!missing(factor)) model <- list(SYMBOL_MULT, model, factor)
  } else if (model[[1]] == SYMBOL_MULT) {
    for (i in 2:length(model)) {
      if (found <- model[[i]][[1]]==R_C) {
        if (length(model[[i]]) != len + 1) stop("bug")
        model[[i]] <- c(R_C, rep(0, preceding), model[[i]][-1])
        if (!missing(factor)) model[[length(model) + 1]] <- factor
        break
      }      
    }
  }
  if (!found) {
    bind <- c(list(R_C), rep(0, preceding), rep(1, len))
    splittingC <- function(model, preceding, factor) {
      const <- sapply(model[-1],
                      function(m) (is.numeric(m) && !is.na(m)) ||
                                  (m[[1]] == R_CONST && !is.na(m[[2]]))
                      )
      if (all(const)) {
        model <- c(model[1], if (preceding > 0) rep(0, preceding), model[-1])
        return(list(SYMBOL_MULT, model, if (!missing(factor)) list(factor)))
      }
      for (i in 2:length(model)) {
        vdim <- preceding + (if (i==2) 0 else GetDimension(model[[i-1]]))
        m  <- ReplaceC(model[[i]])
        L <- GetDimension(m)
        model[[i]] <- setvector(m, preceding = vdim, len = L, factor=factor)
      }
      model[[1]] <- SYMBOL_PLUS
      names(model) <- NULL
      L <- GetDimension(model[[length(model)]])
      model <- SetDimension(model, L)
    }


    if (model[[1]] == SYMBOL_MULT)
      model <- c(SYMBOL_MULT, list(bind), model[-1])
    else model <- list(SYMBOL_MULT, bind, model)
    if (!missing(factor)) model[[length(model) + 1]] <- factor
  }
  return(model)
}


GetDimension <- function(model){
  if (model[[1]] == SYMBOL_PLUS) return(max(sapply(model[-1], GetDimension)))
  if (model[[1]] == R_C) return(length(model) - 1)
  if (model[[1]] == R_CONST) return(length(model[[2]]))
  if (model[[1]] == SYMBOL_MULT) 
    return(max(sapply(model, function(m) if (m[[1]]==R_C) length(m)-1 else 1)))
  return(1)
}


SetDimension <- function(model, L){
  if (model[[1]] == SYMBOL_PLUS) {
    return(c(SYMBOL_PLUS, lapply(model[-1], SetDimension, L=L)))
  }
  if (model[[1]] == R_C) {
    if (length(model) <= L) for (i in length(model) : L) model[[i+1]] <- 0
    names(model) <- c("", letters[1:L])
  } else if (model[[1]] == R_CONST) {
    if (length(model[[2]]) < L)
      model[[2]] <- c(model[[2]], rep(0, L - length(model[[2]])))
  } else if (model[[1]] == SYMBOL_MULT) {
    for (i in 2:length(model)) {
      if (model[[i]][[1]]==R_C) {
        if (length(model[[i]]) <= L)
          for (j in length(model[[i]]) : L) model[[i]][[j+1]] <- 0
        names(model[[i]]) <- c("", letters[1:L])
      }
    }
  }
  return(model)
}


splittingC <- function(model, preceding, factor) {
  ## c(model1,...,model_m) -> c(1,0,..0) * model1 + ... +  c(0,..0,1) * model_m
  ## typically, model_i is not a covariance model, but any shape function
  const <- sapply(model[-1],
                  function(m) {
		    ((is.numeric(m) || is.logical(m)) && !is.na(m)) ||
                      (is.list(m) && m[[1]] == R_CONST && !is.na(m[[2]]))
		  })
  if (all(const)) {
    model <- c(model[1], if (preceding > 0) rep(0, preceding), model[-1])
    return(list(SYMBOL_MULT, model, if (!missing(factor)) list(factor)))
  }
  for (i in 2:length(model)) {
    vdim <- preceding + (if (i==2) 0 else GetDimension(model[[i-1]]))
    m <- model[[i]]
    if (is.list(m)) {
      m  <- ReplaceC(m)
      L <- GetDimension(m)
    } else if (is.numeric(m) || is.logical(m)) {
      L <- length(m)
      if (L ==1 && is.na(m)) {
	m <- list(SYMBOL_MULT, list(R_CONST, m)) 
      } else {
	m <- list(R_CONST, m)
      }
    } else stop(m, "not allowed")
    model[[i]] <- setvector(m, preceding = vdim, len = L, factor=factor)
  }
  model[[1]] <- SYMBOL_PLUS
  names(model) <- NULL
  L <- GetDimension(model[[length(model)]])
  model <- SetDimension(model, L)
}


SplitC <- function(model, factor) {
  model <- splittingC(model, 0, factor)
  L <- GetDimension(model)
  return(SetDimension(model, L))
}



AnyIsNA <- function(model) {
  if (is.list(model)) {
    for (i in 1:length(model)) if (AnyIsNA(model[[i]])) return(TRUE)
    return(FALSE)
  } else return((is.numeric(model) || is.logical(model)) && any(is.na(model)))
}


##RReplaceC <-
ReplaceC <- function(model) {
  if (model[[1]] == SYMBOL_PLUS) {
    for (i in 2:length(model)) model[[i]] <- ReplaceC(model[[i]])
  } else if (model[[1]] == SYMBOL_MULT) {
    if (length(model) <= 2) {
      stop("here, products must have at least 2 factors")
    }
    cs <- sapply(model, function(m) m[[1]] == R_C)
    s <- sum(cs)
    if (s > 0) {
      if (s > 1)
        stop("multiplication with '", R_C, "' may happen only once")
      cs <- which(cs)
      C <- model[[cs]]      
      if (!AnyIsNA(model[-cs])) {
        return(SplitC(C, factor=model[-cs]))
      }
    }      
  } else if (model[[1]] == R_C) {
    return(SplitC(model))
  }
  return(model)
}


initRFlikelihood <- function(model, x, y = NULL, z = NULL, T=NULL, grid=NULL,
                             data, params, distances, dim, likelihood,
                             estimate_variance = NA, ignore.trend=FALSE,
                             RFopt, Reg, ...) {
  ## if changed see also initpredict()

  if (!missing(likelihood)) ## composite likelihood etc
    stop("argument 'likelihood' is a future feature, not programmed yet")

  Z <- UnifyData(model=model, x=x, y=y, z=z, T=T, grid=grid,
                 data=data, distances=distances, dim=dim,
                 RFopt=RFopt, params=params, ...)

  model <- ReplaceC(Z$model); ## multivariates c() aufdroeseln
  
  model <- list("loglikelihood", model, data = Z$data,
                estimate_variance=estimate_variance,
                betas_separate = FALSE,
                ignore_trend=ignore.trend)

  rfInit(mode=model, x=Z$C_coords, reg=Reg, RFopt=RFopt)
}




rflikelihood <- function(model, x, y = NULL, z = NULL, T=NULL, grid=NULL,
                         data, params, distances, dim, likelihood,
                         estimate_variance = NA,
                         Reg, RFopt,
                         ...) {
  initRFlikelihood(model=model, x=x, y = y, z = z, T=T, grid=grid,
                   data=data, params=params,
                   distances=distances, dim=dim, likelihood=likelihood,
                   estimate_variance = estimate_variance,
                   RFopt = RFopt,
                   Reg = Reg, ...)

  likeli <- .Call(C_EvaluateModel, double(0), integer(0), Reg)
  info <- .Call(C_get_likeliinfo, Reg)
  globalvariance <- info$estimate_variance
  where <- 1 + globalvariance
  param <- likeli[-1:-where]
  if (length(param) > 0) names(param) <- info$betanames

  return(list(loglikelihood = likeli[1], likelihood = exp(likeli[1]),
              global.variance = if (globalvariance) likeli[where] else NULL,
              parameters = param
              )
         )
}

RFlikelihood <- function(model, x, y = NULL, z = NULL, T=NULL, grid=NULL,
                         data, params, distances, dim, likelihood,
                         estimate_variance = NA,
                         ...) {
  RFopt <- if (missing(likelihood)) internalRFoptions(...)
           else internalRFoptions(likelihood=likelihood, ...)
  Reg <- RFopt$register$likelihood_register
  if (!hasArg("COPY")) on.exit(optionsDelete(RFopt))
  
  rflikelihood(model=model, x=x, y = y, z = z, T=T, grid=grid,
               data=data, params=params, distances=distances, dim=dim,
               likelihood = likelihood,
               estimate_variance = estimate_variance,
               Reg = Reg, RFopt=RFopt,  ...)
}
  

rfInit <- function(model, x, reg, RFopt, NAOK=TRUE) {
  ##  Print("rfinti", model, x, trafo.to.C_UnifyXT(x))
  if (!is.na(RFopt$basic$seed)) set.seed(RFopt$basic$seed)
  .Call(C_Init, as.integer(reg), model, trafo.to.C_UnifyXT(x), NAOK=NAOK)
  ## returns vdim
}


rfdistr <- function(model, x, q, p, n, params, dim=1, ...) {
  ## note: * if x is a matrix, the points are expected to be given row-wise
  ##       * if not anisotropic Covariance expects distances as values of x
  ##
  ## here, in contrast to Covariance, nonstatCovMatrix needs only x

  RFopt <- internalRFoptions(...)
  if (!hasArg("COPY")) on.exit(optionsDelete(RFopt))

  if (!missing(n) && n>10 && RFopt$internal$examples_reduced) {
    message("number of simulations reduced")
    n <- 10
  }

  model<- list("Distr", PrepareModel2(model, params=params, xdim = 1, ## dummy 
                                      ...)$model,
               dim=as.integer(dim));
  if (!missing(x)) model$x <- if (is.matrix(x)) t(x) else x
  if (!missing(q)) model$q <- if (is.matrix(q)) t(q) else q
  if (!missing(p)) model$p <- if (is.matrix(p)) t(p) else p
  if (!missing(n)) {
    if (exists(".Random.seed") && !is.na(RFopt$basic$seed)) {
      .old.seed <- .Random.seed; on.exit(set.seed(.old.seed), add = TRUE) }
    model$n <- n
  }

  rfInit(model=model, x=UnifyXT(matrix(0, ncol=dim, nrow=1), grid=FALSE),
         reg = MODEL_DISTR, RFopt=RFopt)
  res <-  .Call(C_EvaluateModel, double(0), integer(0),
                as.integer(MODEL_DISTR))

  if (RFopt$general$returncall) attr(res, "call") <-
    as.character(deparse(match.call(call=sys.call(sys.parent()))))
  attr(res, "coord_system") <- c(orig=RFopt$coords$coord_system,
                                 model=RFopt$coords$new_coord_system)
  return(res)
}

RFdistr <- function(model, x, q, p, n, params, dim=1, ...) {
   rfdistr(model=model, x=x, q=q, p=p, n=n, params=params, dim=dim, ...)
}
RFddistr <- function(model, x, params, dim=1, ...) {
  if (hasArg("q") || hasArg("p") || hasArg("n")) stop("unknown argument(s)");
  rfdistr(model=model, x=x, params=params, dim=dim, ...)
}
RFpdistr <- function(model, q, params, dim=1, ...) {
  if (hasArg("x") || hasArg("p") || hasArg("n")) stop("unknown argument(s)");
  rfdistr(model=model, q=q, params=params, dim=dim, ...)
}
RFqdistr <- function(model, p, params, dim=1, ...){
  if (hasArg("x") || hasArg("q") || hasArg("n")) stop("unknown argument(s)");
  rfdistr(model=model, p=p, params=params, dim=dim, ...)
}
RFrdistr <- function(model, n, params, dim=1, ...) {
  if (hasArg("x") || hasArg("q") || hasArg("p")) stop("unknown argument(s)");
  rfdistr(model=model, n=n, params=params, dim=dim, ...)
}


rfeval <- function(model, x, y = NULL, z = NULL, T=NULL, grid=NULL,
                   params=NULL, distances, dim, ...,
                   y.ok = FALSE,
                  ##                  dim = ifelse(is.matrix(x), ncol(x), 1),
                  fctncall=c("Covariance", "CovMatrix", "Fctn",
                             FCTN_TYPE_NAMES), reg=MODEL_USER) {

  ## note: * if x is a matrix, the points are expected to be given row-wise
  ##       * if not anisotropic Covariance expects distances as values of x
  ##
  ## here, in contrast to Covariance, nonstatCovMatrix needs only x

  RFopt <-internalRFoptions(...)
  if (!hasArg("COPY")) on.exit(optionsDelete(RFopt))

  if (is.character(fctncall)) fctncall <- match.arg(fctncall)
  if (fctncall != "CovMatrix" && !missing(distances) && !is.null(distances)) {
    if(missing(dim) || length(dim) != 1) {
      warning("'dim' not given or not of length 1, hence set to 1");
      dim <- 1
    }
    if (length(y) != 0 || length(z) != 0 || length(T) != 0 ||
        (!missing(grid) && length(grid) != 0))
      stop("if distances are given 'y', 'z', 'T', 'grid' may not be given")
    x <- (if (is.matrix(distances)) distances else
          cbind(distances, matrix(0, nrow=length(distances), ncol = dim - 1)))
    distances <- NULL
    dim <- NULL
    grid <- FALSE
  }

  x <- UnifyXT(x=x, y=y, z=z, T=T, grid=grid, distances=distances, dim=dim,
               y.ok = y.ok)
  
  if (!is.list(fctncall)) fctncall <- list(fctncall)

  fctncall[[length(fctncall) + 1]] <-
    PrepareModel2(model, params=params, x=x, ...)$model

  rfInit(model=fctncall, x=x, reg = reg, RFopt=RFopt)
  res <- .Call(C_EvaluateModel, double(0), integer(0), as.integer(reg))

#  print(res)
  
  if (RFopt$general$returncall) attr(res, "call") <-
    as.character(deparse(match.call(call=sys.call(sys.parent()))))
  attr(res, "coord_system") <- .Call(C_GetCoordSystem, reg,
              RFopt$coords$coord_system,
              RFopt$coords$new_coord_system)
   return(res)
}



RFcovmatrix <- function(model, x, y = NULL, z = NULL, T=NULL, grid,
                        params, distances, dim,  ...) {  
  rfeval(model=model, x=x, y=y, z=z, T=T, grid=grid, 
         distances=distances, dim=dim, params=params, ..., fctncall="CovMatrix",
         reg=MODEL_COVMATRIX)
}


#extractVal <- function(L, val) {
#  ans <- NULL
#  for (i in 1:length(L)) {
#    p <- L[[i]]
#    if (is.list(p)) ans <- c(ans, extractVal(L, val))
#    else if (is.numeric(p)) {
#      for (j in 1:length(p)) {
#        q <- pmatch(p[j], val)
#        if (!is.na(q)) ans <- c(ans, q)
#      }
#    } ## else string -- does not matter 
#  }
#  return (ans)
#}



covETC <- function(model, x, y = NULL, z = NULL, T=NULL, grid,
                   params, distances, dim, data, vdim=NULL, ..., alpha) {
  ## Fctn darf intern nicht aufgerufen werden!
  ## ausser von den direkten user-Funktionen
  if (hasArg("data")) {
    if (hasArg("model")) {
      RFfit(model=model, x=x, y=y, z=z, T=T,  grid=grid,
        params=params, distances=distances, dim=dim,
        data = data, methods = NULL, emp_alpha = alpha,
        messages.warn_no_fit = FALSE,
        ...)
    } else {
      internalRFoptions(COPY=TRUE)
      rfempirical(x=x, y=y, z=z, T=T, data=data, grid=grid, vdim=vdim,
                  alpha=alpha, ...)
    }
  } else if (hasArg("x") || hasArg("distances")) { ## no data
    if (alpha == COVARIANCE) {
      fctn <- "Covariance"
      reg <- MODEL_COV
    } else if (alpha == VARIOGRAM) {
      fctn <- "Variogram"
      reg <- MODEL_VARIOGRAM
    } else {
      fctn <- if (alpha > 0) list("Pseudomadogram", alpha = alpha)
              else list("Madogram", alpha = -alpha)
      reg <- MODEL_PSEUDO
    }
##    Print(fctn, alpha)
    rfeval(model=model, x=x, y=y, z=z, T=T, grid=grid, params=params,
           distances=distances, dim=dim, ..., fctncall=fctn, reg= reg,
           y.ok=TRUE)
  } else { ## no data, no coordinates
    if (missing(dim) || length(dim) == 0) dim <- 1
    L <- list(...)
    MARGIN <- if (!is.null(L$MARGIN)) L$MARGIN
              else if (dim < 3) NULL else 1:2
    fixed.MARGIN <- if (!is.null(L$fixed.MARGIN)) L$fixed.MARGIN
                    else if (dim < 3) NULL else rep(0, dim - 2)
    plotmethod <- if (is.null(L$plotmethod)) L$plotmethod
                  else plotmethod <- "none"
    L$plotmethod <- L$fixed.MARGIN <- L$MARGIN <- NULL
    ##Print(alpha)
    do.call(calculateRFplot, c(L,
                               list(x=model,
                                    dim=dim,
                                  fctn.type=alpha,
                                  MARGIN = MARGIN,
                                  fixed.MARGIN=fixed.MARGIN,
                                  plotmethod = plotmethod,
                                  params = if (!missing(params)) params)
                               )
            )
  }
}


RFcovariance <- function(model, x, y = NULL, z = NULL, T=NULL, grid,
		  params, distances, dim, data, vdim=NULL, ...)
  covETC(model=model, x=x, y = y, z = z, T=T, grid=grid,
         params=params, distances=distances, dim=dim, data=data, vdim=vdim, ...,
         alpha = COVARIANCE)

RFcov <- RFcovariance

RFvariogram <- function (model, x, y=NULL, z = NULL, T=NULL, grid,
			 params, distances, dim, data, vdim=NULL, ...)
  covETC(model=model, x=x, y = y, z = z, T=T, grid=grid,
         params=params, distances=distances, dim=dim, data=data, vdim=vdim, ...,
         alpha = VARIOGRAM)

RFpseudovariogram <- function(model, x, y=NULL,  z = NULL, T=NULL, grid,
			      params, distances, dim, data, vdim=NULL, ...)
   covETC(model=model, x=x, y = y, z = z, T=T, grid=grid,
         params=params, distances=distances, dim=dim, data=data, vdim=vdim, ...,
         alpha = PSEUDO)

RFmadogram <- function(model, x, y=NULL,  z = NULL, T=NULL, grid, params,
                       distances, dim,data, vdim=NULL, ..., alpha=1){
  if (alpha <= 0 || alpha > 2) stop("alpha must be in (0,2]")
  covETC(model=model, x=x, y = y, z = z, T=T, grid=grid,
         params=params, distances=distances, dim=dim, data=data, vdim=vdim, ...,
         alpha = -alpha)
}

RFpseudomadogram <- function(model, x, y=NULL,  z = NULL, T=NULL, grid,
                             params, distances, dim, data, vdim=NULL, ...,
                             alpha=1){
  if (alpha <= 0 || alpha > 2) stop("alpha must be in (0,2]")
  covETC(model=model, x=x, y = y, z = z, T=T, grid=grid,
         params=params, distances=distances, dim=dim, data=data, vdim=vdim, ...,
         alpha = alpha)
}


RFfctn <- function(model, x, y=NULL,  z = NULL, T=NULL, grid,
                   params, distances, dim, ...) {
  rfeval(model=model, x=x, y=y, z=z, T=T, grid=grid, params=params,
                distances=distances, dim=dim, ..., fctncall="Fctn",
                reg=MODEL_FCTN )
}

RFcalc <- function(model, params, ...) {
  if (is.numeric(model)) return(model)
  rfeval(model=model, x=0, params=params, ...,
         coord_system="cartesian", new_coord_system="keep", spConform = FALSE,
         fctncall="Fctn",reg=MODEL_CALC)
}


######################################################################
######################################################################


rfDoSimulate <- function(n = 1, reg, spConform) {
  stopifnot(length(n) == 1, n>0, is.finite(n))
  RFopt <- getRFoptions()
  if (missing(spConform)) spConform <- RFopt$general$spConform

  if (RFopt$gauss$paired && (n %% 2 != 0))
    stop("if paired, then n must be an even number")

  info <- RFgetModelInfo_register(reg, level=3)

  vdim <- info$vdim
  totpts <- info$loc$totpts
  if (is.null(totpts) || totpts <= 0)
    stop("register ", reg, " does not look initialized")

  result <- .Call(C_EvaluateModel, as.double(0), as.integer(n),
                  as.integer(reg)) #userdefined,
  if (!spConform) return(result)
  
  prep <- prepare4RFspDataFrame(x=info, RFopt=RFopt) 
  res2 <- conventional2RFspDataFrame(result,
                                     coords=prep$coords,
                                     gridTopology=prep$gridTopology,
                                     n=n,
                                     vdim=info$vdim,
                                     T=info$loc$T,
                                     vdim_close_together
                                     =RFopt$general$vdim_close_together)
  return(res2)
}


    
RFsimulate <- function (model, x, y = NULL, z = NULL, T = NULL, grid=NULL,
                        distances, dim, data, given = NULL, err.model,
                        params, err.params, n = 1, ...) {


  ##  str(model, max=2)
  ## Print("--|+++++--------GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG")
  str(list(...));#ddddd
  if (length(list(...)) > 15) stop("")
   
  if (!missing(model) && is(model, CLASS_FITLIST)) {
    nm <- names(model[METHOD_PREFLIST])
    model <- model[!is.na(nm)]
    if (!is(model, CLASS_CLIST))
      stop("To continue with the output of 'RFfit' ", #use 'predict' or
           "give the component explicitely, e.g., 'RFsimulate(model=",
           deparse(substitute(model)), "[\"",
           "\"], ...)")
  }
  mc <- as.character(deparse(match.call()))

 
### preparations #########################################################  
  stopifnot(is.numeric(n)) ## check whether n is correctly given -- otherwise
  ##            user gets a poor message about expected_number_simu
  if (!missing(distances) && length(distances)  > 0)
    RFoptOld <- internalRFoptions(xyz_notation=length(y)!=0,
                                   expected_number_simu=n, ..., 
                                   general.spConform = FALSE, FORMER=TRUE)
  else {
    RFoptOld <- internalRFoptions(xyz_notation=length(y)!=0,
                                   expected_number_simu=n, ..., FORMER=TRUE)
  }
  
  RFopt <- RFoptOld[[2]]

  if (exists(".Random.seed") && !is.na(RFopt$basic$seed)) {
    .old.seed <- .Random.seed; on.exit(set.seed(.old.seed), add = TRUE) }

  if (n>2 && RFopt$internal$examples_reduced) {
    message("number of simulations reduced")
    n <- 2
  }

  cond.simu <- !missing(data) && length(data) > 0
  reg <- if (cond.simu) RFopt$registers$predict_register
         else RFopt$registers$register

  ### simulate from stored register ########################################
  mcall <- as.list(match.call(expand.dots=FALSE))
  if (length(mcall)==1 ||
      length(mcall)==2 && !is.null(mcall$n) ||
      length(mcall)==3 && !is.null(mcall$n) && "..." %in% names(mcall)) {
    if (cond.simu) {
      ## TO DO
      stop("repeated performance of conditional simulation not programmed yet")
    } else {
      res <- rfDoSimulate(n=n, reg=reg, spConform=RFopt$general$spConform
                          #userdefined=userdefined
                          )
      if (RFopt$general$returncall) attr(res, "call") <- mc
      attr(res, "coord_system") <- .Call(C_GetCoordSystem, reg,
                                     RFopt$coords$coord_system,
                                     RFopt$coords$new_coord_system)
      return(res)
    }
  }

### preparations #########################################################
  stopifnot(!missing(model) && !is.null(model))

  if (!missing(x)) x <- UnifyXT(x=x, y=y, z=z, T=T, grid=grid)
  
  model.orig <- model
  pM2 <- PrepareModel2(model, params=params, x = if (missing(x)) NULL else x,
                       ..., dont.add.data=!cond.simu)
   model <- pM2$model
 
   ### conditional simulation ###############################################
  if (cond.simu) {
    if (n < 1) return(NULL)
    if (missing(err.model) || length(err.model) ==0) err.model <-NULL
    if (isSpObj(data)) data <- sp2RF(data)
    stopifnot(missing(distances) || is.null(distances))
    res <- switch(GetProcessType(model),
                  RPgauss =  rfCondGauss(model=model.orig, x=x, n=n, data=data,
                                         given=given, err.model=err.model,
                                         params = params,
                                         ## next line to ensure that this part
                                         ## matches with predictGauss
                                         ...),
                  stop(GetProcessType(model),
                       ": conditional simulation of the process not programmed yet")
                  )
    
  } else { ## unconditional simulation ####    
    if(!missing(err.model) && !is.null(err.model))
      warning("error model is unused in unconditional simulation")
    warn.seed.not.na(RFoptOld)

    ## ACHTUNG! KRIGING ruft auf rfInit("Simulate") auf !! Abgleich
    ##          gegebenenfalls notwendig
    rfInit(model=list("Simulate",
                      setseed=eval(parse(text="quote(set.seed(seed=seed))")),
                      env=.GlobalEnv, model), x=x, reg=reg, RFopt=RFopt)
    if (n < 1) return(NULL)
    res <- rfDoSimulate(n=n, reg=reg, spConform=FALSE)
  } # end of uncond simu

  
  ## output: RFsp   #################################
  if (((!missing(x) && length(x) != 0)) ## not imputing
      && RFopt$general$spConform) {

    info <- RFgetModelInfo_register(reg, level=3)
    
    if (length(res) > 1e7) {
      message("Too big data set (more than 1e7 entries) to allow for 'spConform=TRUE'. So the data are returned as if 'spConform=FALSE'")
      return(res)
    }
    
    prep <- prepare4RFspDataFrame(x=x, RFopt=RFopt,
				  coordnames = pM2$data.names[pM2$is.x])

     
    res <- conventional2RFspDataFrame(data=res, coords=prep$coords,
                                      gridTopology=prep$gridTopology,
                                      n=n,
                                      vdim=info$vdim[1],
                                      T=x$T,
                                      vdim_close_together=
                                      RFopt$general$vdim_close_together)
    if (is.raster(x)) {
      res <- raster::raster(res)
      raster::projection(res) <- raster::projection(x)
    }

    attributes(res)$varnames <- pM2$data.names[pM2$is.var]
  }

   
  if (RFopt$general$returncall) attr(res, "call") <- mc
  attr(res, "coord_system") <- .Call(C_GetCoordSystem,
                                     as.integer(reg),
                                     RFopt$coords$coord_system,
                                     RFopt$coords$new_coord_system)

  return(res)
}

					# RFreplace <- function(model, by) { }


RFplot <- function(x, ...) {
  if (missing(x)) {
    cat("Empirical Variogram:\n"); print(args(RFplotModel))## ok
    cat("\nSimulations:\n"); print(args(RFplotSimulation))## ok
    cat("\nModels:\n(Note that x can be a list of models; or further models can be given\nas additional arguments whose names starting with 'model')\n");
    print(args(RFplotModel)) ## ok
  } else if (is(x, CLASS_FITLIST) || is(x, CLASS_EMPIR)) {
    RFplotEmpVariogram(x=x, ...)
  } else if (is(x, "RFdataFrame") || is(x, "RFspatialGridDataFrame") ||
           is(x, "RFspatialPointsDataFrame")) RFplotSimulation(x=x, ...)
#  else if (is(x, CLASS_PLOT)) RFplotModel(x, ...)
  else RFplotModel(x=x, ...)
}
