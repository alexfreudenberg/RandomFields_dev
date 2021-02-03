# Authors 
## Martin Schlather, schlather@math.uni-mannheim.de
##
##
## Copyright (C) 2015 -- 2017 Martin Schlather
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of th GNU General Public License
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


## source("modelling.R")


initpredict <- function(conditioning, C_coords, grid=NULL, given, model, data, 
                        ignore.trend,
                        RFopt, Reg, ...) {
  conditioning <- ReplaceC(conditioning) ## multivariates c() aufdroeseln

  predict <- list("predict",
                  conditioning=conditioning,
                  predict=if (missing(model)) conditioning else ReplaceC(model),
                  data = data,
                  estimate_variance=NA,
                  betas_separate = FALSE,
                  ignore_trend=ignore.trend,
                  given=given)
  rfInit(model=predict, x=C_coords, reg=Reg, RFopt=RFopt)
  ## naechste Zeile notwendig!?!
  ##  if (!ignore.trend) .Call(C_EvaluateModel, double(0), as.integer(1), Reg)
  NULL
}



ModelAbbreviations <- function(model){
  L <- length(model)
  N <- character(L)
  for (mm in 1:L) {
    m <- unlist(model[[mm]])
    dollar <- m==RM_S[1] | m==RM_S[2]
    idx <- which( names(m) != "" & !dollar) ## names(m):names of parameters
    for (s in which(dollar)) {
      z <- idx[min(which(idx > s))] ## to which model belongs the dollar param?
      m[s] <- m[z] ## put the name of subsequent model in position of dollar pos
      m <- m[-z] ## delete subseqeuent model
    }
    idx <- names(m) != ""
    m[idx] <- paste(names(m)[idx], m[idx], sep="=")
    N[mm] <- paste(abbreviate(m), collapse =".")
  }
  return(N)
}


ModelParts <- function(model, effects, complete) { ## model immer schon aufbrtt
  ## wird in kriging und in RFfit.R aufgerufen!!
  
  ## model <- P repareModel2(model) 13.7.19
  model <- if ((model[[1]] %in% SYMBOL_PLUS)) model[-1] else list(model)
  if (complete) {
    return(model[effects >= RandomEffect])
  } else {
    err <- model[effects == ErrorEffect]
    if (length(err) == 1) err <- err[[1]]
    else if (length(err) == 0) err <- NULL
    else err <- c(SYMBOL_PLUS, err)
    
    m <- model[ effects == RandomEffect ]
    if (length(m) > 0) {
      if (length(m) == 1) m <- m[[1]] else m <- c(SYMBOL_PLUS, m)
      return(list(model=m, err.model=err))
    } else {
      return(list(model=err, err.model=NULL))
    }
  }
}

"
xRFranef <- function(fit, method='ml', OP='@') {
  Z <- do.call(OP, list(object, 'Z')
  model <- fit[method]
  parts <- ModelParts(model, effect=GetModelInfos(Z)$effect, complete=TRUE)
  
  ans <- lapply(parts, function(parts)
    if (TRUE)
      RFinterpolate(COPY=FALSE,
                    model, x = Z$coord, data = Z$data, given = Z$coord,
                    err.model=parts)
    else ## oder?! 29.12.20
      RFinterpolate(COPY=FALSE,
                    parts, x = Z$coord, data = Z$data, given = Z$coord,
                    err.model=model) ## ?? 29.12.20
    )
  return(ans)
}
"

FinImputIntern <- function(data, simu, coords, coordnames, data.column, vdim,
                           spConform, fillall=FALSE) {
  data <- data.matrix(data)
  n <- length(data) / (vdim * coords$totpts)
  if (is(data, "RFsp")) {
    if (spConform) {
      data@data[ , ] <- as.vector(simu)
      return(data)
    } else {
      values <- as.matrix(data@data)
      values[is.na(values) | fillall] <- simu
      return(cbind(coordinates(data), values))
    }
  } else { ## not RFsp
    if (coords$grid) {
      ## to do
      stop("not programmed yet")
    } else {
      ##  coords <- all$x
      colnames(coords$x) <- coordnames
      
      values <- data[, data.column]
      values[is.na(values) | fillall] <- simu
      
      if (!spConform)  return(cbind(coords$x, values))
      
      tmp.all <- conventional2RFspDataFrame(data=values, coords=coords$x,
                                            gridTopology=NULL,
                                            n=n, vdim=vdim,
                                            vdim_close_together=FALSE)
      if (is(tmp.all, "RFspatialPointsDataFrame"))
        Try(tmp.all <- as(tmp.all, "RFspatialGridDataFrame"))
      if (is(tmp.all, "RFpointsDataFrame"))
        Try(tmp.all <- as(tmp.all, "RFgridDataFrame"))
    }
    return(tmp.all)
  }
}


FinishImputing <- function(data, simu, Z, spConform, fillall) {
  ## to do: grid
  is.var <- Z$is.var.orig
  vdim <- Z$vdim
  if (is.list(data)) {
    for (i in 1:length(data))
      data[[i]] <- FinImputIntern(data=data[[i]], simu=simu[[i]],
                                  coords=Z$coords[[i]], coordnames=Z$coordnames,
                                  data.column=is.var, vdim=vdim,
                                  spConform = spConform, fillall=fillall)
    return(data)
  }

  return(FinImputIntern(data=data, simu=simu, coords=Z$coords[[1]],
                        coordnames=Z$coordnames, data.column=is.var, 
                        vdim=vdim, spConform=spConform, fillall=fillall))

}



ExpandGrid <- function(x) {
  #### ACHTUNG! ZWINGENDE REIHENFOLGE
  if (x$grid) { # 0
    x$x <-
      as.matrix(do.call(expand.grid,
                        lapply(apply(cbind(x$x, x$T), 2,
                                     function(x) list(seq(x[1],by=x[2],length.out=x[3]))), function(x) x[[1]])))
  } else if (x$has.time.comp) {
    dim.x <- if (is.vector(x$x)) c(length(x$x), 1) else dim(x$x)
    x$x <- cbind(matrix(rep(t(x$x), times=x$T[3]),
                            ncol=dim.x[2], byrow=FALSE),
                     rep(seq(x$T[1], by=x$T[2],
                             length.out=x$T[3]), each=dim.x[1]))
  }
  if (length(x$y) > 0) stop("no expansion within a kernel definition")
#  x$y <- double(0) #1 
  x$T <- double(0) #2
  x$grid <- FALSE  #3
#  x$spatialdim <- ncol(x$x) #4
  x$has.time.comp <- FALSE           #5
#  x$di st.given <- FALSE      #6
  x$totpts <- nrow(x$x)   #7
  return(x)
}



rfPrepareData <- function(model, x, y=NULL, z=NULL, T=NULL,
                          distances=NULL, dim, grid,
                          data, given=NULL, params,
                          RFopt, reg, err = NULL,
                          err.params, 
                          ...) {
  ## NOTE: if err is a list of list, the err is
  ## NOT added to the conditioning model, but both are taken as they are.
  ## AND these two models exchange their meaning to
  ## 
  ## internal behaviour to get the random effects with the
  ## existing algorithm!!

                                        #  apart_linear_fit <- RFopt$krige$apart_linear_fit ## might be an option to the user

  if (!missing(distances) && length(distances)>0)
    stop("option distances not programmed yet.")

  missing.x <- missing(x) || length(x) == 0
  imputing <- missing.x && (missing(distances) || length(distances) == 0)
  new <- NULL
  
  if (!imputing) {
    ## here, other values of covariates given in ... might be used
    new <- UnifyXT(x, y, z, T, grid=grid, distances=distances, dim=dim)
  }

  ## NOTE: new-x is for prediction; model-x is "given"
  ## only if not "given", then x=new is used, assuming
  ## that given and x are the same; or imputing
  further.models <- NULL
  genuine.err <- length(err) >= 1 && !is.numeric(err) && !is.logical(err)
  ##             from RFLight: err=NA 
  if (genuine.err) {
    add <- is(err, CLASS_CLIST)
    if (missing(err.params) || length(err.params) == 0) further.models <- err
  }

  if (length(given) == 0) {
    ## so either within the data or the same the x-values (some kind of imputing)
    Z <- UnifyData(model=model, data=data, ## ergo, x muss in data sein!
                   RFopt=RFopt, params=params, dont.add.data = FALSE,
                   further.models=further.models, ...)
     
    if (Z$matrix.indep.of.x.assumed) {
      if (imputing) stop("coordinates cannot be detected") ## if (missing.x) ?
      Z <- UnifyData(model=model, x=new, data=data, RFopt=RFopt, params=params,
                     dont.add.data = FALSE, further.models=further.models,
                     ...)
    }
  } else {
    Z <- UnifyData(model=model, data=data, x=given, RFopt=RFopt,
                   params=params, dont.add.data = FALSE,
                   further.models=further.models, ...)
  }

  
  model <- Z$model
  conditioning <- if (length(Z$further.models) > 0) Z$further.models[[1]]
                  else model

  if (length(Z$data) != 1) stop("exactly one data set expected.", CONTACT)
  dim_data <- base::dim(Z$data[[1]]) ## loc x (vdim * repet)
  Z$data[[1]] <- as.double(Z$data[[1]])
  repet <- Z$repetitions
  new.dim_data <- c(prod(dim_data) / repet, repet)# (loc * vdim) x repet
  base::dim(Z$data[[1]]) <- new.dim_data # (loc * vdim) x repet
  fit.these.na <- is.na(Z$data[[1]])# (loc * vdim) x repet; boolean, no NAs
  any.na.in.repeated <- rowSums(fit.these.na) > 0# (loc * vdim); integer, no NAs
  base::dim(Z$data[[1]]) <- dim_data # orig: loc x (vdim * repet)
  base::dim(any.na.in.repeated) <-
    c(length(any.na.in.repeated) / Z$vdim , Z$vdim)# loc x vdim
  any.na.at.loc <- rowSums(any.na.in.repeated > 0) > 0# loc x 1; integer, no NAs
  any.na <- any(any.na.at.loc) 

  if (any.na && Z$coords[[1]]$dist.given)
    stop("missing values not programmed yet for given distancs")

  if (imputing) {
      if (RFopt$krige$fillall || !any.na) {
        fit.these.na <- rep(TRUE, length(fit.these.na)) ## nur
        any.na <- TRUE
    ## um Daten im Endergebnis einzutragen
        new <- Z$coords[[1]]
      } else {
        new <- ExpandGrid(Z$coords[[1]])
        new$x <- new$x[any.na.at.loc, , drop=FALSE]
      }
  } else {
    if (Z$tsdim != new$spatialdim + new$has.time.comp)
      stop("coodinate dimension of locations with and without data, ",
           "respectively, do not match.")
    if (RFopt$general$na_rm_lines) { ## i.e. rm complete locations
      Z$data[[1]] <- Z$data[[1]][!any.na.at.loc, , drop=FALSE]
      Z$coords[[1]] <- ExpandGrid(Z$coords[[1]])
      Z$coords[[1]]$x <- Z$coords[[1]]$x[!any.na.at.loc,  , drop=FALSE]
    } else if (!Z$coords[[1]]$grid && !Z$coords[[1]]$has.time.comp)
      Z$coords[[1]] <- ExpandGrid(Z$coords[[1]]) ## to do in variogram.cc
  }

  if (genuine.err) {
    ## 29.12.20 the definitions have completely changed to version 4.0
    ## formerly list of err have been interpreted completely differently
    if (!missing(err.params) && length(err.params) > 0) {
      ## UnifyData muss aufgerufen werden, da covariaten auch im
      ## error model sein koennten (z.B. die Varianz)
      if (!add) err <- err[[1]] ## err has covering list()
      setRFoptions(coords.coordnames = Z$coordnames, coords.varnames=Z$varnames)
      conditioning <- UnifyData(model=err, data=data, x=given, RFopt=RFopt,
                                params=params, dont.add.data = FALSE,...)$model
    }
    if (add && !RFopt$basic$skipchecks) {
      linpart <- RFlinearpart(COPY=FALSE, model=err,
                              params=err.params, x=Z$coords, set=1, ...)
      if (length(linpart$X) > 0 || any(linpart$Y != 0))
        stop("a trend is not allowed for the error model.")
    }
    if (add) conditioning <- list(SYMBOL_PLUS, model, conditioning)
  }

  return(list(Z=Z,  ## given coordinates
              X=new, ## new coordinates
              conditioning = conditioning, # model for the matrix to be
              ##                            inverted in SK
              model=model, # model for the vector in SK
              imputing=imputing,
              any.na.at.loc = any.na.at.loc,
              fit.these.na =fit.these.na))
}


RFinterpolate <- function(model, x, y=NULL, z=NULL, T=NULL, grid=NULL,
                          distances, dim, data, given=NULL, params,
                          err.model=NULL, err.params,                          
                          ignore.trend=FALSE, ...) {
  
  RFopt <- internalRFoptions(...,
                             return_variance = FALSE) ## currently!! TO DO !!
  if (!hasArg("COPY")) on.exit(optionsDelete(RFopt))
  
  if (length(err.model) > 1 && is.list(err.model) &&
      is(err.model[[1]], CLASS_CLIST)) {
    stop("this variation of RFinterpolate has to be progreammed yet")
    names <- ModelAbbreviations(err.model)
    ans <- vector("list", length(err.model))    
    for (i in 1:length(ans)) {
      ## check RandomFieldsLight, warum err.model zu model wird
      ## und umgekehrt -- klaeren warum!      
      ans[[i]] <- RFinterpolate(COPY = FALSE,
                                model=err.model[[i]], x=x, y=y, z=z,
                                T=T, grid=grid,
                                distances=distances, dim=dim,
                                data=data, given=given, params=params,
                                err.model=model,
                                err.params=err.params[[i]], 
                                ignore.trend=ignore.trend,
                                ...)
      attr(ans[[i]], "model") <- names[i]
    }
    return(ans)
  }

  err.NA <- (length(err.model) > 0 &&
             (is.numeric(err.model) || is.logical(err.model) ))
 
  #Print(err.NA, is.list(data) && length(data) > 1 && !err.NA)


  if (!is.data.frame(data) && is.list(data) && length(data) > 1 && !err.NA) {
    ## err.NA treated below, should be the outer loop
    if (is.list(x)) { # y,z,T, must be NULL, but this not checked here
      ans <- vector("list", length(data))
      len <- length(x)
      for (i in 1:length(data))
        ans[[i]] <- RFinterpolate(COPY = FALSE,
                                  model=model, x=x[[((i-1) %% len) + 1]],
                                  y=y, z=z, T=T, grid=grid,
                                  distances=distances, dim=dim,
                                  data=data[[i]], given=given, params=params,
                                  err.model=err.model,
                                  err.params=err.params, 
                                  ignore.trend=ignore.trend, ...)
    } else {      
      for (i in 1:length(data))
        ans[[i]] <- RFinterpolate(COPY = FALSE,
                                  model=model, x=x,
                                  y=y, z=z, T=T, grid=grid,
                                  distances=distances, dim=dim,
                                  data=data[[i]], given=given, params=params,
                                  err.model=err.model,
                                  err.params=err.params, 
                                  ignore.trend=ignore.trend, ...)
    }
    return(ans)
  }

  
  if (!missing(distances) && length(distances) > 0)
    stop("'distances' not programmed yet.")
  
  
  boxcox <- .Call(C_get_boxcox, NULL)


  ## eingabe wird anstonsten auch als vdim_close erwartet --
  ## dies ist nocht nicht programmiert! Ausgabe ist schon programmiert
  ## CondSimu ist auch noch nicht programmiert
  if (RFopt$general$vdim_close_together)
    stop("'vdim_close_together' must be FALSE")

  
  return_variance <- RFopt$krige$return_variance

  reg <- RFopt$registers$predict_register
  all <- rfPrepareData(model=model, x=x, y=y, z=z, T=T,
                       distances=distances, dim=dim, grid=grid,
                       data=data, given=given, params=params, RFopt=RFopt,
                       reg=reg, err = err.model, err.params=err.params,
                       ...)

  if (err.NA) {# GetModelEffects
    if (length(err.NA) != 1 || !is.na(err.NA)) stop("'err.model' invalid.")
    model <- list(ModelParts(model, effects=GetModelEffects(all$Z),
                             complete = FALSE)$model)    
    ans <- vector("list", length(err.model))    
    for (i in 1:length(ans)) {
      ## check RandomFieldsLight, warum err.model zu model wird
      ## und umgekehrt -- klaeren warum!
      ans[[i]] <- RFinterpolate(COPY=FALSE,
                                model=err.model[[i]], x=x, y=y, z=z,
                                T=T, grid=grid,
                                distances=distances, dim=dim,
                                data=data, given=given, params=params,
                                err.model=model,
                                err.params=err.params[[i]], 
                                ignore.trend=ignore.trend, ...)
    }
    return(ans)
  }

  imputing <- all$imputing
  tsdim <- as.integer(all$Z$tsdim)
  repet <- as.integer(all$Z$repetitions)
  vdim <- all$Z$vdim
  if (!imputing) {
    coordnames.incl.T <-
      c(if (!is.null(all$Z$coordnames)) all$Z$coordnames else
        paste(COORD_NAMES_GENERAL[1], 1:all$Z$spatialdim, sep=""),
        if (all$Z$has.time.comp) COORD_NAMES_GENERAL[2] else NULL)
    if (all$X$grid) {
      coords <- list(x=NULL, T=NULL)
      xgr <- cbind(all$X$x, all$X$T)

      colnames(xgr) <- coordnames.incl.T
      gridTopology <- sp::GridTopology(xgr[1, ], xgr[2, ], xgr[3, ])
      ## bis 3.0.70 hier eine alternative
    } else {
      coords <- list(x=all$X$x, T=all$X$T)
      ## wenn bei gegeben unklar was zu tun ist. Ansonsten
      if (length(coords$T) == 0)  colnames(coords$x) <- coordnames.incl.T
      gridTopology <- NULL
    }
  }
  nx <- all$X$totpts
  dimension <-
    if (all$X$grid) c(if (length(all$X$x) > 0) all$X$x[3, ],
                      if (length(all$X$T) > 0) all$X$T[3]) else nx # to do:grid
  newdim <- c(dimension, if (vdim>1) vdim, if (repet>1) repet)

  if (imputing && return_variance) {
    return_variance <- FALSE
    warning("with imputing, currently the variance cannot be returned")
  }
    
  exact <- RFopt$general$exact
  maxn <- RFopt$krige$locmaxn
  ngiven <- as.integer(all$Z$coords[[1]]$totpts) ## number of given points
  split <- RFopt$krige$locsplitn[NEIGHB_SPLIT+1]
  
  split <- ngiven > maxn || (!is.na(exact) && !exact && ngiven > split)
  data <- boxcoxIntern(all$Z$data[[1]], boxcox=boxcox, vdim=vdim)
    
  Res <- if (imputing) data else matrix(nrow=nx, ncol=repet * vdim)
  sigma2 <- NULL  ## currently just a dummy
  
  all.at.once <- TRUE ## vielleicht in Zukunft auch FALSE moeglich:
  ## Vorteil: linearer part wird en bloc geschaetzt und nicht als
  ##          unabhaengige Teile
  ## Nachteil: grosse Matrix muss fuer den Schaetzer des linearen Term
  ##           invertiert werden
  ## Somit: nur interessant wenn all$X sehr viele Punkte umfasst
  ##        und Z$coord hoch aber nicht zu hoch (5000er Bereich)
  splitted.data <- data
  splitted.X <- all$X
  splitted.Z <- all$Z$coords
  if (split) {
    ## to do:
    all$X <- ExpandGrid(all$X) ## to  do
    all$Z$coords[[1]] <- ExpandGrid(all$Z$coords[[1]]) ## to  do
    
    ## neighbourhood kriging !
    if (!is.na(exact) && exact)
      stop("number of conditioning locations too large for an exact result.")
    if (ngiven > maxn && is.na(exact) &&
        RFopt$basic$printlevel>=PL_IMPORTANT)
      message("performing neighbourhood kriging")
    
    ## calculate the boxes for the locations where we will interpolate
    idx <- GetNeighbourhoods(
        model = all$model, # [[1]],
        Z=all$Z,
        X=all$X, ## newlocations; to do: grid
        splitfactor=RFopt$krige$locsplitfactor,
        maxn=RFopt$krige$locmaxn,
        split_vec = RFopt$krige$locsplitn,
        shared = TRUE
    )
    totalparts <- length(idx[[2]])      
    givenidx <- lapply(idx[[2]], function(p) unlist(idx[[1]][p]))
    if (totalparts > 1) {
      if (all.at.once) {
        setRFoptions(general.pch="")
        splitted.data <- lapply(givenidx, function(p) data[p, , drop=FALSE])
        splitted.X <- splitUnifyXT(x=splitted.X, split=idx[[3]])
        splitted.Z <- splitUnifyXT(x=splitted.Z, split=givenidx)
      } # else: gar nichts machen -- wird bei for-Schleifen bzw C gemacht
    } else all.at.once <- TRUE
    pr <- totalparts > 1 && RFopt$general$pch != "" &&RFopt$general$pch != " "
  } else {
    totalparts <- 1
    idx <- list(list(TRUE), list(TRUE), list(TRUE))
    pr <- FALSE
  }

  initpredict(conditioning=all$conditioning,# including error structure
              model=all$model,
              grid=FALSE,
              C_coords = C_UnifyXT(splitted.X),
              given = C_UnifyXT(splitted.Z),
              data=splitted.data,
              ignore.trend = ignore.trend,
              Reg=reg, RFopt=RFopt)

  .Call(C_set_boxcox, c(Inf, 0), reg)

##  Print("x")
  
  for (p in 1:totalparts) {
    ##    Print(p, totalparts)
    if (pr && p %% 5==0) cat(RFopt$general$pch)
    res <- if (!split || all.at.once)
             ## EvaluateModel/gauss_predict schafft nur 1 Datenssatz zu
             ## verarbeiten, insbesondere da EvaluateModel nur
             ## vector of doubles zurueckliefern kann.
             ## Rechenzeittechnisch sollte aber der wiederholte Aufruf
             ## von R nicht so viel ausmachen
             .Call(C_EvaluateModel, as.double(0), p, reg)
           else ## nur in ansatzen programmiert :
             .Call(C_EvaluateModel,
                     as.double(c(length(givenidx[[p]]),
                                 length(idx[[3]][[p]]))),
                   as.integer(c(givenidx[[p]], idx[[3]][[p]])), reg)

    whereto <- idx[[3]][[p]]
    isNA <- is.na(Res[whereto, ])

  ##  Print("yx")
  

# Print(whereto, isNA); Print(Res, Res[whereto, ][isNA]); Print(res, res[isNA])
    
    Res[whereto, ][isNA] <- res[isNA]        
  } 
  if (pr) cat("\n")

  Z <- all$Z ## achtung! oben kann sich noch all$Z aendern!
  X <- all$X

  if (FALSE) {## jonas
    if (return_variance && length(newdim <-
                                    c(dimension, if (vdim>1) vdim)) > 1)
      base::dim(sigma2) <- newdim
  }

  if (length(newdim)>1) base::dim(Res) <- newdim
  else Res <- as.vector(Res)
  if (!is.null(Z$varnames)) attributes(Res)$varnames <- Z$varnames
  Res <- boxcoxIntern(data=Res, boxcox = boxcox, inverse=TRUE, vdim=vdim)
  
 
  spConform <- RFopt$general$spConform
  if (imputing) {
    Res <- FinishImputing(data=data, simu=Res, Z=Z,
                          spConform=spConform,
                          fillall = RFopt$krige$fillall) ## to do : grid
    if (return_variance){
      var <- FinishImputing(data=data, simu=sigma2, Z=Z,
                            spConform=spConform,
                            fillall = RFopt$krige$fillall)# to do : grid
      if (spConform) {
        names(var@data) <- paste("var.", names(var@data), sep="")
          Res@.RFparams$has.variance <- TRUE
        Res <-  cbind(Res, var)
      } else Res <- list(Res, var)
    }
    return(Res)
  } else {
    if (!spConform) {
      if (vdim > 1 && RFopt$general$vdim_close_together) {
        Resperm <- c(length(dimension)+1, 1:length(dimension),
                     if(repet>1) length(dimension)+2)
        Res <- aperm(Res, perm=Resperm)
        if (return_variance)
          sigma2 <- aperm(sigma2, perm=Resperm[1:(length(dimension)+1)])
      }
      return(if (return_variance) list(estim = Res, var = sigma2) else Res)
    }
    
    Res <- conventional2RFspDataFrame(Res, coords=coords$x,
                                      gridTopology=gridTopology,
                                      n=repet, vdim=vdim, T = coords$T,
                                      vdim_close_together =
                                        RFopt$general$vdim_close_together)
    
    if (return_variance){
      var <- conventional2RFspDataFrame(sigma2, coords=coords$x,
                                        gridTopology=gridTopology,
                                        n=1, vdim=vdim, T = coords$T,
                                        vdim_close_together =
                                          RFopt$general$vdim_close_together)
      names(var@data) <- paste("var.", names(var@data), sep="")
      Res@.RFparams$has.variance <- TRUE
      Res <-  cbind(Res, var)
    }
  }
  
#  Res@.RFparams$krige.method <-
#    c("Simple Kriging", "Ordinary Kriging", "Kriging the Mean",
#      "Universal Kriging", "Intrinsic Kriging")[krige.meth.nr]
  

  ## Res@.RFparams$var <- sigma2 ## sehr unelegant.
  ## * plot(Res) sollte zwei Bilder zeigen
  ## * var(Res) sollte sigma2 zurueckliefern
  ## * summary(Res) auch summary der varianz, falls vorhanden
  ## * summary(Res) auch die Kriging methode

  if (is.raster(x)) {
    Res <- raster::raster(Res)
    raster::projection(Res) <- raster::projection(x)
  }
  
  return(Res)
}

rfSimulate <- function(model, x, reg, RFopt, n) {
  rfInit(model=list("Simulate", model, env=.GlobalEnv,
                    setseed=eval(parse(text="quote(set.seed(seed=seed))"))),
           x=x, reg=reg, RFopt=RFopt)
  res <- rfDoSimulate(n=n, reg=reg, spConform=FALSE)
  return(res)
}

rfCondGauss <- function(model, x, n=1,
                        data,   # first coordinates, then data
                        given=NULL, ## alternative for coordinates of data
                        params=NULL,
                        err.model=NULL, err.params, ...) { # ... wegen der Variablen
  .Call(C_setlocalRFutils, NA, NULL)
  RFopt <- getRFoptions()
   
  reg <- RFopt$registers$predict_register
  boxcox <- .Call(C_get_boxcox, NULL)

  if (length(err.model) > 0 &&
      (is.list(err.model) && !is.character(err.model[[1]]))) {
    stop("TO DO.", CONTACT)
  }

  err.NA <- (length(err.model) > 0 &&
             (is.numeric(err.model) || is.logical(err.model) ))

  if ((is.list(data) && !is.data.frame(data)) && length(data) > 1 && !err.NA) {
    stop("TO DO;", CONTACT)
  }

  all <- rfPrepareData(model=model, x=x, data=data, given=given, params=params,
                       RFopt=RFopt, reg=reg,
                       err = err.model, err.params=err.params,...)

  if (err.NA) {     
    if (length(err.NA) != 1 || !is.na(err.NA)) stop("'err.model' invalid.")
    stop("TO DO:", CONTACT)
 }
  
  Z <- all$Z
  X <- all$X
  simu.grid <- X$grid
  tsdim <- Z$tsdim
  vdim <- Z$vdim

  data <- boxcoxIntern(Z$data[[1]], boxcox=boxcox, vdim=vdim)
  if (all$Z$repetitions != 1)
     stop("conditional simulation allows only for a single data set")
  
  txt <- "kriging in space time dimensions>3 where not all the point ly on a grid is not possible yet"
  ## if 4 dimensional then the last coordinates should ly on a grid

  ## now check whether and if so, which of the given points belong to the
  ## points where conditional simulation takes place

  duplicated_loc.old <- duplicated_loc <- 
    which(RFopt$general$duplicated_loc == DUPLICATEDLOC_NAMES) - 1
  if (duplicated_loc.old == DUPLICATEDLOC_REPEATED && !hasArg(err.model))
    stop("in case of 'duplicated_locations == \"error\"' an error model must be given.")
  if (duplicated_loc == DUPLICATEDLOC_ERROR) {
    duplicated_loc <- if (length(err.model) > 0) DUPLICATEDLOC_REPEATED
                      else DUPLICATEDLOC_SCATTER
    setRFoptions(general.duplicated_locations =
                   DUPLICATEDLOC_NAMES[duplicated_loc + 1])
  }

  coords <- ExpandGrid(Z$coords[[1]]) ## conditioning locations
  simu <- NULL
   if (simu.grid) { ## simulation on a grid
    xgr <- cbind(X$x, X$T)
    l <- ncol(xgr)
    ind <- 1 + (t(coords$x) - xgr[1, ]) / xgr[2, ] 
    index <- round(ind)
    outside.grid <- (duplicated_loc == DUPLICATEDLOC_REPEATED ||
                     duplicated_loc == DUPLICATEDLOC_RISKERROR) | 
      apply((abs(ind-index)>RFopt$general$gridtolerance) | (index<1) |
            (index > 1 + xgr[3, ]), 2, any)

    if (any(outside.grid)) {     
      ## at least some data points are not on the grid:
      ## simulate as if there is no grid

      if (!all(outside.grid))
        stop("locations of given data and those to be predicted partially match. This cannot be treated currently with the current value of the option 'duplicated_locations'. See man pages for alternative values.") ## TO DO
      
      simu.grid <- FALSE
 
      if (l>3) stop(txt)
      xx <- if (l==1) ## dim x locations
             matrix(seq(from=xgr[1], by=xgr[2], len=xgr[3]),
                        nrow=1)
            else eval(parse(text=paste("t(expand.grid(",
                            paste("seq(from=xgr[1,", 1:l, 
                                  "], by=xgr[2,", 1:l,
                                  "], len=xgr[3,", 1:l, "])", collapse=","),
                         "))")))  
      ll <- eval(parse(text=paste("c(",
                   paste("length(seq(from=xgr[1,", 1:l, 
	                 "], by=xgr[2,", 1:l, 
		         "], len=xgr[3,", 1:l, "]))",
                         collapse=","),
                   ")")))

      new.index <- rep(0,ncol(index))
      ## data points that are on the grid, must be registered,
      ## so that they can be used as conditioning points of the grid
      if (!all(outside.grid)) {
        new.index[!outside.grid] <- 1 +
          colSums((index[, !outside.grid, drop=FALSE]-1) *
                  cumprod(c(1, ll[-length(ll)])))
      }
      index <- new.index
      new.index <- NULL
    } else {
      ##      Print("GRID")
      ## data points are all lying on the grid
      simu <- rfSimulate(model=all$conditioning,
                         x=UnifyXT(X$x, T=X$T, grid=X$grid),
                         n=n, reg=reg, RFopt=RFopt)
 
      ## for all the other cases of simulation see, below
      if (is.vector(simu)) dim(simu) <- c(length(simu), 1)
      else if (!is.matrix(simu)) { ## 3.1.26 is programmed differently
        nvdim <- (vdim > 1) + (n > 1)
        if (nvdim > 0) {
          d <- dim(simu)
          last <- length(d) + 1 - (1 : nvdim)
          dim(simu) <- c(prod(d[-last]), prod(d[last]))
        } else dim(simu) <- c(length(simu), 1)
      }
      
      cum <- cumprod(c(1, xgr[3, -ncol(xgr)]))
      index <- as.vector(1 + cum %*% (index - 1))

      if (is.vector(simu)) dim(simu) <- c(length(simu), 1)
      total <- dim(simu)[1]
      simu.given <-  do.call("[", c(list(simu, index, drop=FALSE),
                                    as.list(rep(TRUE,length(dim(simu))-1))))
     }
  } else { ## not simu.grid
    xx <- t(X$x)  ## dim x locations
   
    ## the next step can be pretty time consuming!!!
    ## to.do: programme it in C
    ##
    ## identification of the points that are given twice, as points to
    ## be simulated and as data points (up to a tolerance distance !)
    ## this is important in case of nugget effect, since otherwise
    ## the user will be surprised not to get the value of the data at
    ## that point
    one2ncol.xx <- 1:ncol(xx)
    if (duplicated_loc == DUPLICATEDLOC_REPEATED)
      index <-rep(0, nrow(coords$x))
    else if (duplicated_loc == DUPLICATEDLOC_SCATTER) {
      ## nothing to be done since already scattered at unifyData
    } else {
      index <- apply(coords$x, 1, function(u){
               i <- one2ncol.xx[colSums(abs(xx-u)) <RFopt$general$gridtolerance]
               if (length(i)==0) return(0)
               if (length(i)==1) return(i)
               return(NA)
      })
    }
  }

    if (any(index == 0) && any(index != 0))
      stop("locations of given data and those to be predicted partially match. This cannot be treated currently with the current value of the option 'duplicated_locations'. See man pages for alternative values.") ## TO DO

  
  if (!simu.grid) {
    ## otherwise the simulation has already been performed (see above)
    tol <- RFopt$general$gridtolerance * nrow(xx)
    if (any(is.na(index)))
      stop("identification of the given data points is not unique - `tol' too large?")
    if (any(notfound <- (index==0))) {
      index[notfound] <- (ncol(xx) + 1) : (ncol(xx) + sum(notfound))
    }

    xx <- rbind(t(xx), coords$x[notfound, , drop=FALSE])
    simu <- rfSimulate(model=all$conditioning, x=UnifyXT(xx, grid=FALSE), n=n,
                       RFopt = RFopt, reg=reg)


    if (is.vector(simu)) dim(simu) <- c(length(simu), 1)
    simu.given <-  do.call("[", c(list(simu, index, drop=FALSE),
                                      as.list(rep(TRUE, length(dim(simu))-1))))

    simu <- do.call("[", c(list(simu, 1:X$totpts, drop=FALSE),
                           as.list(rep(TRUE, length(dim(simu))-1))))
  }

  ## to do: als Naeherung bei UK, OK:
  ## kriging(data, method="A") + simu - kriging(simu, method="O") !

  d <- dim(data)
  data <- as.vector(data) - simu.given
  dim(data) <- c(d, length(simu.given) / prod(d))
  
  stopifnot(length(X$y)==0, length(X$z)==0)

  setRFoptions(general.duplicated_locations =
                DUPLICATEDLOC_NAMES[duplicated_loc.old + 1])

##  Print(all)
  ## TO DO RFinterpolate und koenn(t)en dasselbe register verwenden,
  ## d.h. rfSimulate die initialisierung von RFinterpolate.
  ## dazu muss bei RFinterpolate die  initialisierung intern ausgelagert werden.
  ## und von condsimu aufgerufen werden.

  ## ACHTUNG! alles wird hier ueberschrieben, da bei predict_register verwenden!
  interpol <- RFinterpolate(COPY=FALSE, x=X, # OK
                            model=all$model,
                            err.model = list(all$conditioning),
                            given = coords,
                            data = data,                            
                            spConform=FALSE, ignore.trend = TRUE,
                            boxcox = c(Inf, 0)
                            )
  
  simu <- as.vector(simu) + as.vector(interpol)
  dim(simu) <- c(if (X$grid) X$x[3,] else X$totpts,
                 if (vdim>1) vdim, if (n > 1) n)
  
  simu <- boxcoxIntern(data=simu, boxcox = boxcox, inverse=TRUE, vdim=vdim)
  
  if (all$imputing) {
    return(FinishImputing(data=Z$data[[1]], simu=simu, Z=Z,
                          spConform=RFopt$general$spConform,
                          fillall = RFopt$krige$fillall))
  }

  attributes(simu)$varnames <- Z$varnames
  attributes(simu)$coordnames <- Z$coordnames


  return(simu)  
}
