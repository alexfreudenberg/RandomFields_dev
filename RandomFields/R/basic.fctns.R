
## Authors 
## Martin Schlather, schlather@math.uni-mannheim.de
##
##
## Copyright (C) 2015 -- 2017  Martin Schlather
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




## FUNCTION-STARP***********************************************************************************
## NAME		extract VarNames
## REQUIRE	$model is a formula or a RMmodel
## ENSURE		the return value is a vector with the names of the response variables of the formula 
#			or NULL if no left side given
## SEE		
## AUTHOR		Sebastian Gross 
## DATE		26.08.2011; 2014 Martin Schlather modified
## FUNCTION-END*************************************************************************************



##OneTo <- function(n) return(if (n < 1) NULL else 1:n)
OneTo <- function(n) return(if (length(n) > 1) stop("invalid end of for loop") else if (n < 1) NULL else 1:n)

TwoTo <- function(n) return(if (length(n) > 1) stop("invalid end of for loop") else if (n < 2) NULL else 2:n)

Try <- function(expr) {
  z <- tryCatch(expr, error = function(e) e)
  if (is.list(z) && !is.null(z$message) && !is.null(z$call))
    class(z) <- CLASS_TRYERROR
  z
}

rawTry <- function(expr) tryCatch(expr, error = function(e) e)
rawError <- is.list

XXXextractVarNames <- function(model) {
  if (missing(model) || length(model) <=2 || !is(model, "formula"))
    return(character(0))
  tmp <- as.character(model)[2]
  tmp <- strsplit(tmp, "c\\(")[[1]]
  tmp <- paste(tmp, sep="", collapse="")
  tmp <- strsplit(tmp, "\\)")[[1]]
  tmp <- paste(tmp, sep="", collapse="")
  varnames <- strsplit(tmp, ", ")[[1]]

  ## ignore numeric formated varnames
  i<- 1
  while (i <= length(varnames)) {
    if (regexpr("^[[:digit:]]", varnames[i]) > 0) varnames <- varnames[-i]
    else i <- i+ 1
  }
  if (length(varnames) == 0) return (NULL)
  return (varnames)
}


earth_coord_names<- function(names, opt) {
  ## Earth coordinates + possibly radius
  n <- substr(tolower(names), 1, 6)
  nc <- nchar(n)
  lon <- lat <- logical(length(n))
  for (i in 1:length(n)) {
    lon[i] <- substr(opt$earth_coord_names[1], 1, nc[i]) == n[i]
    lat[i] <- substr(opt$earth_coord_names[2], 1, nc[i]) == n[i]
  }
  lonORlat <- lon | lat  
  earth <- all(nc[lonORlat] >= 2) && sum(lon==1) && sum(lat == 1)
  if (!earth || (lo <- which(lon)) > (la <- which(lat))) return(NULL)
  c(lo, la)
}

cartesian_coord_names <- function(names, opt) {
  n <- tolower(names)
  coords <- tolower(opt$cartesian_names[c(4, 1:3)]) # c("T", "x", "y", "z")
  Txyz <- outer(n, coords, "==")
  cs <- colSums(Txyz)
  if (any(cs > 1) || sum(cs[1:2]) == 0 || any(diff(cs[-1]) > 0))
    return (integer(0))
  Txyz <- Txyz[, c(2:4, 1), drop=FALSE]
  ord <- apply(Txyz, 2, function(x) which(x > 0))
  ord <- order(unlist(ord))
  rs <- which(rowSums(Txyz) > 0)
  return(rs[ord])
}


general_coord_names <- function(names) {
  n <- substr(tolower(names), 1, 5)
  return(which(n == "coord"))
}


extractFromNames <- function(what="var", RFopt, cn, ncol=length(cn)) {
  ## this function partially interprets RFoptions$coords$varnames/varidx,
  ## namely when varnames is not NA
  ##  Print(cn)

  varidx <- RFopt$coord[[paste0(what, "idx")]]
  if (!is.na(varidx[1])) {
    is.var <- varidx
    if (is.na(is.var[2])) is.var[2] <- ncol
    ans <- is.var[1] : is.var[2]
    if (length(cn) != 0) names(ans) <- cn[ans]
    return(ans)
  }
  if (length(cn) == 0) return(NULL)
  
  varnames <- RFopt$coords[[paste0(what, "names")]]
  ans <- NULL
  if (what == "var") { 
    if ((vdim <- length(varnames)) >  0) {
      if (RFopt$general$vdim_close_together)
        stop("'vdim_close_together' must be FALSE")
      l <- vector("list", vdim)
      for (v in 1:vdim)
        l[[v]] <- substring(cn, 1, nchar(varnames[v])) == varnames[v]
      repet <- sapply(l, sum)
      if (repet[1] > 0) {
        if (any(repet != repet[1]))
          stop("detected repetitions are not equal for all components")
        m <- matrix(unlist(l), ncol=vdim)
        if (any(rowSums(m) > 1))
          stop("names of multivariate components not unique")        
        ans <- apply(m, 2, which)
        ans <- if (is.vector(ans)) as.matrix(ans) else t(ans)
        rownames(ans) <- varnames
      }
    }
    if (length(ans) == 0) {
      cnlow <- tolower(cn)
      ans <- which(substr(cnlow, 1, 4) == "data" |
                       substr(cnlow, 1, 4) == "value" |
                       substr(cnlow, 1, 8) == "variable")
      ##  Print(cn, cnlow , cn[ans], ans, strsplit(cn[ans], ".", fixed=TRUE))

      if (length(ans) == 0) return(NULL)
      s <- strsplit(cn[ans], ".", fixed=TRUE)
      if (length(s[[1]]) > 1) { ## check for repeated measurments
        cn <- sapply(s, function(x) paste0(x[-length(x)], collapse=".")) # cut last part
        if ((rep <- sum(cn[1] == cn)) > 1) {
          size <- ncol / rep
          if (size == as.integer(size) && all(cn == cn[1 : size])) {
            dim(ans) <- c(size, rep)
            rownames(ans) <- cn[1 : size]
          }
        }
      }
    }
  } else { ## "coord"
    ans <- which(cn %in% varnames)
    if (length(ans) == 0 &&
        length(ans <- earth_coord_names(cn, opt=RFopt$coords)) != 2 &&
        length(ans <- cartesian_coord_names(cn, opt=RFopt$coords)) == 0)
      ans <- general_coord_names(cn)

    if (length(ans) > 1 && length(unique(cn[ans])) < length(ans))
      stop("column names in the data frame that indicate coordinates may not be doubled")
  }
  
 ## Print(what, ans, cn, length(ans) > 0  && !is.matrix(ans))
  stopifnot(all(ans <= length(cn)))
  if(length(ans) > 0  && !is.matrix(ans)) names(ans) <- cn[ans]  
##  print(ans)
  return(ans)
}



S <- function(x) if (length(x) > 1) "s" else ""
ARE <- function(x) if (length(x) > 1) "are" else "is"


data.columns <- function(data, model=NULL, xdim=0, x=NULL,
                         RFopt=getRFoptions(),
                         force=FALSE, halt=TRUE, 
			 vdim=0, ...,
                         params=NULL, 
                         only.data=FALSE,
                         return_transform = FALSE) {
##  Print(data, force, model, xdim)
  ##  Print(halt, vdim,  xdim)
  ##  Print(list(...))

##  Print(data, x)
  
  ## out of the model only the trend part is needed here
  ## so model can be the trend only
  ##
  ## this function interprets RFoptions$coords$varnames/varidx or model
  ## also if varnames is NA
  ## gives the columns in the original data

  PL <-  as.integer(RFopt$basic$printlevel)
  coord.opt <- RFopt$coords
  orig.data <- data

  if (length(xdim) == 0) xdim <- 0
  data <- if (is.list(data) && !is.data.frame(data)) data[[1]] else data
  data <- if (is(data, "RFsp")) data@data else data
  cn <- colnames(data)
  data.info <- "safe"
  if (x.given <- hasArg("x") || hasArg("distances")) {
    ## besser mit ...names, so dass missing(param) nicht zu einem Fehler fuehrt
    L <- list(...)
    x.given <- length(L$x) + length(L$distances) > 0
  }
  searching.for.x <- !only.data && !x.given
  if (length(vdim) == 0) vdim <- 0
  repet <- 0
  if (!missing(model) && !is.null(model)) {
    M <- PrepareModel2(model=model, params=params,
                       x=x, ...,  data=orig.data, xdim=xdim,
                       return_transform = return_transform)
    components <- c("is.x", "is.factor", "is.var", "is.covariate",
                    "is.unclear")
    for (i in components) assign(i, c(M[[i]]))# resolve in case of matrix
    if (length(is.x) > 0 && x.given) stop("bug in 'data.matrix'")
    if (M$repet > 0) repet <- M$repet
    if (M$vdim > 0) {
      if (vdim > 0 && vdim != M$vdim)
        stop("the given value for 'vdim=", vdim,
             "' does not match the detected value: vdim=", M$vdim)
      vdim <- M$vdim
    }

    for (i in components) {
      if (length(M[[i]]) > 0) {
        if (i == "is.unclear")
          Note("detection", "the model does not use '",
               paste(M[[i]], collapse="', '"), "'.")
        else  Note("detection", "the model uses the following ",
                   switch(i, "is.x"="coordinate",
                          "is.factor"="factor",
                          "is.var"="variable",
                          "is.covariate"="covariate",
                          stop(CONTACT)),
                   S(M[[i]]),  ": '",
                   paste(M$data.names[M[[i]]], collapse="', '"), "'.")
      }
    }
  } else { ## model not given
    is.var <- extractFromNames(RFopt=RFopt, cn=cn, ncol=ncol(data))
    if (is.matrix(is.var)) {
      if (vdim != 0 && vdim != nrow(is.var))
        stop("multivariate dimension of data unclear")
      repet <- ncol(is.var)
      vdim <- nrow(is.var)
      is.var <- as.vector(is.var)
    }
    is.x <- if (!searching.for.x) NULL
            else extractFromNames("coord", RFopt=RFopt, cn=cn, ncol=ncol(data))

    is.unclear <- if (length(cn) > 0) cn[-c(is.var, is.x)]
    is.factor <- is.covariate <- NULL

    M <- list()
  }
  
  if (length(is.var) == 0) {
    is.var <- if (ncol(data) == 1) 1L
              else if (vdim > 0 && ncol(data) == vdim) 1:vdim
  }

  if ((n <- length(is.var)) > 0) {
    d <- min(c(is.var, which(cn == "")))
    m <- max(0, is.x, is.covariate, is.factor)
    if (m < d) {
      if (is.var[1] == 1 && vdim  == 0) is.var <- d : ncol(data) 
      else if ((ncol(data) - d + 1) %% vdim == 0 &&
               (max(is.var) - min(is.var) + 1) <= vdim) {
        if (repet == 0) repet <- (ncol(data) - d + 1) /  vdim
        is.var <- is.var + ((0:(repet-1)) * vdim)
      }
      data.info <-  "1stsafe"      
      if (n != length(is.var)) {
        if (length(is.var) < 20) {
          Note("detection", "The column",S(is.var), " ",
               if (length(cn) == 0) paste(is.var, collapse=",")
               else paste("'", cn[is.var],"'", sep="",collapse=", "),
               " ", ARE(is.var), " considered as ",
               if (data.info=="1stsafe") "potential ",
               "data column",S(is.var), ".")
        } else {
          nnn <- 2
          show.data <- is.var[c(1:nnn, length(is.var) + 1 - (nnn:1))] 
          if (length(cn) > 0) show.data <- paste0("'", cn[show.data],"'")
          Note("detection", "The columns ",
               paste(show.data[1:nnn], collapse=","), ", ..., ",
             paste(show.data[-1:-nnn], collapse=","),
             " are considered as ", if (data.info=="1stsafe") "potential ",
             "data columns.")
        }
      }
    }
  }

  if (isRFsp <- is(data, "RFsp")) {
    xdim <- if (is(data, "SpatialPointsDataFrame") ||
		is(data, "RFpointsDataFrame")) ncol(data@coords)
            else if (is(data, "SpatialGridDataFrame") ||
		     is(data, "RFgridDataFrame")) length(data@grid@cells.dim)
	    else stop("unknown class for 'data'")  
  } else if (searching.for.x) {
    if (xdim > 0 && xdim >= ncol(data)) stop("not enough columns in 'data'.")
    if ((len <- length(is.x)) > 0) {
      if (length(is.var) > 0 && data.info=="safe")
        is.x <- setdiff(is.x, is.var)
      if (length(is.x) < len)
        Note("detection", "column", S(is.x), " ", paste(is.x, collapse=",")," ",
             ARE(is.x), " assumed to define the coordinate",S(is.x),
             ", since the following keyword",S(is.x),
             "for coordinates ", ARE(is.x),
             " recognized: ", paste(names(is.x), collapse=","), "\n")
    } else if (PL >= PL_SUBIMPORTANT)
      Note("detection", "None of the column names belong to the default ",
           "coordinate names of RandomFields.")
  }


  ######### NOTPROGRAM ###############
  ##
  if (length(is.var)==0 && length(is.x)==0) {
    data.info <- "guess"
    if (is.null(cn) && !x.given) {
      if (!force) {
          if (halt)
            stop(if (is.null(cn)) "colnames of data argument must contain"
                 else 'no colname starts with', ' "data" or "variable"')
          else return(list(is.var=TRUE, is.x=NULL));
      }
      is.var <- (xdim+1):ncol(data)
    } else {
      if (!searching.for.x) is.var <- 1:ncol(data)
      else stop("data name(s) could not be found in the column names")
    }
    if (searching.for.x)
      Note("detection", "taking column", S(is.var),
           paste(is.var, collapse=","),
           " as data column", S(is.var))
  }

  if (!x.given) {
    if (length(is.x) == 0) {
      if (searching.for.x && !isRFsp) {
        is.x <- (1:ncol(data))[-c(is.var, is.factor)]
        if (xdim > 0) {
          if (length(is.x) < xdim)
            stop("not enough columns for coordinates found ")
          if (xdim < length(is.x) && PL >= PL_SUBIMPORTANT)
            Note("detection", "column(s) ",
                    paste("'", is.x[-1:-xdim], "'", , sep="", collapse=", "),
                    "in data frame unused.\n")
          is.x <- is.x[1:xdim]
        }
        if (length(is.x) > 0)
          Note("detection", "Using ",
               if (length(cn)>0) paste("'", cn[is.x], "'", sep="",collapse=", ")
               else paste("column(s) ", paste(is.x, collapse=",")),
               " as coordinate(s).")
      }
    } else if (xdim > 0 && xdim != length(is.x))
      stop("expected dimension of coordinates does not match the found coordinates")
  }

  
  if (length(is.var) == 0) {  #  (all(is.na(info$varnames))) 
    is.var <- (1:ncol(data))[-is.x]
    if (length(is.var) == 0) stop("no columns for data found")
  } else {
    if (any(is.x %in% is.var)) stop("column names and data names overlap.")
    if ( (isRFsp && length(is.var) < ncol(data)) ||
	 (!isRFsp && length(is.x) + length(is.var) < ncol(data)) ) {
      if (PL >= PL_SUBIMPORTANT)
	Note("detection", "column(s) ",
             paste((1:ncol(data))[-c(is.x, is.var, is.factor)],
                   collapse=", "), " unused.\n")
    }
  }

  if (length(vdim) > 0 && vdim > 1 && length(is.var) %% vdim != 0)
    stop("recognized data columns not a multiple a multiple does not match the multivariate model")

  if (length(cn) > 0) {
    if (length(is.var) > 0) names(is.var) <- cn[is.var] 
    if (length(is.x) > 0) names(is.x) <- cn[is.x]
  }

  #  Print(is.var=is.var, is.x=is.x, is.factor=is.factor,vdim=vdim, repet = repet, data.info=data.info, model=M)

  return(list(is.var=is.var, is.x=is.x, is.factor=is.factor,
              vdim=vdim, repet = repet,
              data.info=data.info, model=M))
}


    
SystemCoordNames <- function(locinfo, opt) {
  ## this function tries to combine all information available on
  ## coordinate names and variable names und returns the names if
  ## ensured that the names are the correct one.
  has.time.comp <- locinfo$has.time.comp
  tsdim <- locinfo$spatialdim + locinfo$has.time.comp
 
  system <- opt$coord_system
  if (system == "earth") {
    coordnames <- if (tsdim == 4) opt$earth_coord_names
		  else if (tsdim == 2) opt$earth_coord_names[1:2]
		  else if (tsdim == 3) c(opt$earth_coord_names[1:2], "HeightOrTime")
  } else if (system == "cartesian" && tsdim <= 4) {
    coordnames <- opt$cartesian_names[1:tsdim]
    if (has.time.comp) coordnames[tsdim] <- opt$cartesian_names[3 + 1]
  } else {
    coordnames <- paste(COORD_NAMES_GENERAL[1], 1:tsdim, sep="")
    if (has.time.comp) coordnames[tsdim] <- COORD_NAMES_GENERAL[2]
  }
  return(coordnames)
}



search.model.name <- function(cov, name, level) {
  if (length(name) == 0 || length(cov) ==0) return(cov);
  if (!is.na(pmatch(name[1], cov))) return(search.model.name(cov, name[-1], 1))

  for (i in 1:length(cov$submodels)) {
    found <- search.model.name(cov$submodels[[i]], name, 1)
    if (!is.null(found)) return(found)      
  }
  found <- search.model.name(cov$internal, name, 1)
  if (!is.null(found)) return(found)
  if (level == 0) stop("model name not found")
  return(NULL)
}






vectordist <- function(x, diag=FALSE) {
  storage.mode(x) <- "double"
  res <- .Call(C_vectordist, t(x), diag)
  dimnames(res) <- list(dimnames(x)[[2]], NULL)
  return(t(res));
}



my.legend <- function(lu.x, lu.y, zlim, col, cex=1, ...) {
  ## uses already the legend code of R-1.3.0
  cn <- length(col)
  if (cn < 43) {
    col <- rep(col, each=ceiling(43 / cn))
    cn <- length(col)
  }
  filler <- vector("character", length=(cn-3)/2)
  legend(lu.x, lu.y, y.intersp=0.03, x.intersp=0.1, 
         legend=c(format(zlim[2], dig=2), filler,
             format(mean(zlim), dig=2), filler,
             format(zlim[1], dig=2)),
         lty=1, col=rev(col),cex=cex, ...)
}

