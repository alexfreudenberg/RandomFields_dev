
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

R.p <- function(proj, new, factor) iR.p(proj=proj, new=new, factor=factor)
R.p <- copyProp(R.p, iR.p)


RRdistr <- function(name, nrow, ncol,  ## ddistr, pdistr, qdistr, rdistr,
		    envir,  ...) {

#  Print(missing(name), if (!missing(name)) name, missing(nrow), if (!missing(nrow)) nrow, missing(ncol), if (!missing(ncol)) ncol,  missing(envir), ...)
  
  if (!missing(name)) {
    u <- rawTry(isS4(name) && is(name, CLASS_CLIST))
    if (is.logical(u) && u) return(name)
  }

  par.general <- submodels <- par.model <- list()

  if (FALSE)
  if (!missing(ddistr) && length(ddistr) != 0) {
    stopifnot(!missing(ddistr), !missing(pdistr), !missing(qdistr),
	      !missing(rdistr), !missing(envir))
    par.model <- c(list(name=name, nrow=nrow, ncol=ncol,
			ddistr=ddistr, pdistr=pdistr, qdistr=qdistr,
			rdistr=rdistr, envir=envir),
		   list(...))
    model <- new(CLASS_CLIST, name = RM_DISTR[1],
		 submodels = submodels, 
		 par.model = par.model, par.general = par.general)
    return(model) 
  }

  ll <- as.list(substitute(name))
  name <- ll[[1]]
  if (is.character(name)) ll <- list(...)
  else {
    name <- as.character(name)
    ll <- ll[-1]
  }

  if (length(ll) > 0) {
    par.names <- names(ll)
    if (length(par.names) == 0 || any(par.names == ""))
      stop("In distribution families\n  (i) all parameters must be named;\n (ii) in more complicated models all distributions taken from R must be\n      encapsulated by 'RRdistr';\n(iii) use 'RRloc' to modifiy scale and location.\nIf the call of a distribution family was not the intension, the error could be\n  (i) missing ~")
    num <- sapply(ll, function(x) 
      is.numeric(x) || is.symbol(x) || !(is(Try(eval(x)), "try-error"))
      )
    if (!all(num)) {
      subs <- ll[!num]
      n <- names(subs)
      for (i in 1:length(subs)) {
        if (!is.language(subs[[i]]))
          stop("type of parameter (function, constant) cannot be determined")
        par.model[[n[i]]] <-
          if (substr(deparse(subs[[i]]), 1, 1)=='R') eval(subs[[i]]) else
              do.call(RM_DISTR[1], list(subs[[i]]))
      }
    }
    if (any(num)) {
      ll <- ll[num]
      n <- names(ll)
      for (i in 1:length(ll)) par.model[n[i]] <- eval(ll[[i]])
    }
  }
  
  var <- c('x', 'q', 'p', 'n')
  fctnames <-  c('d', 'p', 'q', 'r')
  pm <- list()
  pm[['name']] <- name
  if (hasArg(ncol)) pm[['ncol']] <- ncol
  if (hasArg(nrow)) pm[['nrow']] <- nrow
  for (ii in 1:length(fctnames)) {
    i <- fctnames[ii]
    iname <- paste(i, name, sep="")
    if (!exists(iname)) {
      stop("'", name,
     "' is considered as a name for a distribution family,\n  but the function '",
           iname, "' is not visible.")
    }
    pm[[paste(i, "distr", sep="")]] <-
      eval(parse(text=paste("quote(", iname, "(",
                   var[ii], "=", var[ii],
                   if (length(ll) > 0 && length(par.names)>0)
                      paste(",",
                            paste(par.names, "=", par.names, collapse=", ")),
                   "))", sep="")))
  }
  pm[['envir']] <- if (hasArg(envir)) envir else new.env()
  par.model <- c(pm, par.model)
    
  model <- new(CLASS_CLIST, name = RM_DISTR[1], submodels = submodels, 
               par.model = par.model, par.general = par.general)
  return(model) 
}

RRdistr <- new(CLASS_RM,
               .Data = RRdistr,
               type = TYPE_NAMES[RandomType + 1],
               domain = DOMAIN_NAMES[DOMAIN_MISMATCH + 1],
               isotropy = ISO_NAMES[ISO_MISMATCH + 1],
               operator = FALSE,
               monotone =  MONOTONE_NAMES[NOT_MONOTONE + 1 - MON_UNSET],
               finiterange = FALSE,
               simpleArguments = FALSE,
               maxdim = Inf,
               vdim = 1
               )


GetSymbols <- function(ll) {
  idx <- sapply(ll, function(x) is.symbol(x) || !is.language(x))
  symbols <- as.character(ll[idx])
  if (!all(idx)) 
    for (i in which(!idx)) {
      symbols <- c(symbols, GetSymbols(as.list(ll[[i]])))
    }
  return(symbols)
}
       

RMuser <- function(type, domain, isotropy, vdim, beta,
                   coordnames = c("x", "y", "z", "T"),
                   fctn, fst, snd, envir,
                   var, scale, Aniso, proj) {
	submodels <- par.general <- par.model <- list() 
	
	if (!hasArg(type)) type <- TYPE_NAMES[ShapeType + 1]
        if (is.numeric(type)) par.model[['type']] <- type
        else if (is.character(type))
          par.model[['type']] <- pmatch(type, TYPE_NAMES) - 1

        
        if (any(par.model[['type']] < ProcessType))
          message("It is likely that the function you are defining is already available in 'RandomFields', or hasn't got the claimed property ('",
                  TYPE_NAMES[1+par.model[['type']]],
                  "'). (If you are not sure whether your function is positive/negative definite, please contact ", AUTHOR, ")\nUsing predefined functions leads to (much!) shorter computing times (up to a factor 100).\nSee ?RM for an overview over the implemented models. Further,\nsome simulation methods do not work at all for user defined functions.")
        else if (any(par.model[['type']] == TrendType))
          message("Please make sure that the defined function is not available in 'RandomFields'.\nUsing predefined functions leads to (much!) shorter computing times (up to a factor 100).\nSee ?RMmodelsTrend  for an overview over the implemented models. Further,\nsome simulation methods do not work at all for user defined functions.");

	if (hasArg(domain)) {
	  if (is.numeric(domain)) par.model[['domain']] <- domain
	  else if (is.character(domain))
                   par.model[['domain']] <- pmatch(domain, DOMAIN_NAMES) - 1
	  else stop("wrong value for 'domain'")
	}
	if (hasArg(isotropy)) {
	  if (is.numeric(isotropy)) par.model[['isotropy']] <- isotropy
	  else if (is.character(isotropy))
                   par.model[['isotropy']] <- pmatch(isotropy, ISO_NAMES) - 1
	  else stop("wrong value for 'isotropy'")
	}

	if (hasArg(vdim)) {
	  if (is.numeric(vdim)) par.model[['vdim']] <- vdim
	  else stop("wrong value for 'vdim'")
	}
	if (hasArg(beta)) {
	  if (is.numeric(beta) || is.language(beta) || 
	      is.list(beta))
	     par.model[['beta']] <- beta
	  else if (substr(deparse(substitute(beta)), 1, 1)=='R')
	    submodels[['beta']] <- beta
	  else submodels[['beta']] <- RRdistr(beta)
	}
 	if (hasArg(fctn)) {
          f <- substitute(fctn)
          par <- coordnames %in% GetSymbols(as.list(as.list(f)[-1]))
          par.model[['fctn']] <- f
 	} else stop("'fctn' not given")
       
	if (hasArg(fst)) {
          f <- substitute(fst)
          if (any(xor(par,
                      coordnames %in% GetSymbols(as.list(as.list(f)[-1])))))
            stop("the variables in 'fst' do not match the ones in 'fctn'")
          par.model[['fst']] <- f
	} else if (hasArg(snd)) stop("'fst' not given")
        
	if (hasArg(snd)) {         
	  f <- substitute(snd)          
          if (any(xor(par,
                      coordnames %in% GetSymbols(as.list(as.list(f)[-1])))))
            stop("the variables in 'snd' do not match the ones in 'fctn'")
	  par.model[['snd']] <- f
	}

        par.model[['coordnames']] <- which(par)    
	par.model[['envir']] <- if (hasArg(envir)) envir else new.env()
	par.general[['var']] <-if (hasArg(var)) var else NO_DOLLAR_VALUE
	par.general[['scale']] <-if (hasArg(scale)) scale else NO_DOLLAR_VALUE
	par.general[['Aniso']] <-if (hasArg(Aniso)) Aniso else NO_DOLLAR_VALUE
	par.general[['proj']] <-if (hasArg(proj)) proj else NO_DOLLAR_VALUE

	model <- new(CLASS_CLIST, name = RM_USER[1], 
			submodels = submodels, 
			par.model = par.model, par.general = par.general)
	return(model) 
 } 

RMuser <- new(CLASS_RM,
              .Data = RMuser,
              type = TYPE_NAMES[PosDefType + 1],
              domain = DOMAIN_NAMES[PREVMODEL_D + 1],
              isotropy = ISO_NAMES[PREVMODEL_I + 1],
              operator = FALSE,
              monotone =  MONOTONE_NAMES[NOT_MONOTONE + 1 - MON_UNSET], # [MON_PARAMETER]
              finiterange = TRUE,
              simpleArguments = FALSE,
              maxdim = Inf,
              vdim = -1
              )



RMdeclare <- function(...) {
  CL <- match.call()
  cl <- as.list(CL)
  model.name <- as.character(cl[[1]])
  cl <- cl[-1]
  right <- sapply(cl, as.character)
  left <- names(cl)
  if (length(left) == 0) left <- rep("", length(right))
  ok <- (left == right) | (left == "")
  ## Print(model.name, cl, right, left, ok)
  if (!all(ok)) {
    left <- left[!ok]
    stop("argument", S(left),  " ",  paste0("'", left, "'", collapse =","),
         " do", if (length(left) == 1) "es",
         " not follow the requirements of 'RMdeclare'.")
  }
  par.model <- list()
  par.model[right] <- NA
  ##  Print(left, right, par.model)
  left <- left[left == right]
  for (i in OneTo(length(left))) {
##    Print(left[i], base::...elt(i))
    par.model[[left[i]]] <- base::...elt(i)
  }
##   Print(left, par.model)
  model <- new(CLASS_CLIST,
               name = model.name, ## == RMdeclare
               submodels = list(), 
               par.model = par.model, par.general = list())
  
  return(model) 
}

RMdeclare <- new(CLASS_RM,
               .Data = RMdeclare,
               type = TYPE_NAMES[TrendType + 1],
               domain = DOMAIN_NAMES[XONLY+ 1],
               isotropy = ISO_NAMES[PREVMODEL_I + 1],
               operator = FALSE,
               monotone =  MONOTONE_NAMES[NOT_MONOTONE + 1],
               finiterange = FALSE,
               simpleArguments = TRUE,
               maxdim = Inf,
               vdim = -2
               )



RMcovariate <- function(formula=NULL, data, x, y=NULL, z=NULL, T=NULL, grid,
                        raw, addNA, factor) {
  if (!missing(factor)) {
    if (!missing(addNA) && addNA)
      stop("'addNA' and 'factor' may not be given at the same time.")
    isna <- is.na(factor)
    if (any(xor(isna[1], isna))) stop("If 'factor' has NAs then all of the values must be NAs")
  }
  if (missing(data)) {
    if ("formula" %in% class(formula))
      stop("data argument 'data' has not been given")
    data <- formula
    formula <- NULL
  }

  N <- "data.frame"
  if (is.matrix(data)) {
    if (!is.null(colnames(data)) && !is.null(formula)) data <- as.data.frame(data)
  } else if (is.factor(data)) {
    if (!is.null(formula)) stop("if 'data' is a factor, a formula may not be given")
    data <- data.frame(own.factor = data)
    N <- "factor"
  }

  if (is.data.frame(data)) {
    if (is.null(formula)) 
      formula <- eval(parse(text=paste("~", paste(names(data), collapse="+"))))
    data <- model.matrix(object=formula, data=data)
    if (getRFoptions(getoptions_="basic")$printlevel > 0)
      message("Intercept added to a model component that is a ", N, ".")
  }
  
  Call <- iRMcovariate
  if (missing(x) && length(T)==0) {
    if (length(y)!=0 || length(T)!=0 || !missing(grid))
      stop("y, z, T, grid may only be given if 'x' is given")
    ans <- Call(data=data, raw=raw, addNA=addNA)
  } else {
    PL <- getRFoptions(getoptions_="basic")$printlevel
    .Call(C_setlocalRFutils, NULL, 0)
    new <- C_UnifyXT(x=x, y=y, z=z, T=T, grid=grid)
    .Call(C_setlocalRFutils, NULL, PL)
    ans <- Call(data=data, x=new, raw=raw, addNA=addNA)
  }
  ans
}
RMcovariate <- copyProp(RMcovariate, iRMcovariate)
 
RMfixcov <- function(M, x, y=NULL, z=NULL, T=NULL, grid, var, proj, raw#, norm
                     ) {
  Call <- iRMfixcov
  if (missing(x) && length(T)==0) {
    if (length(y)!=0 || length(T)!=0 || !missing(grid))
      stop("y, z, T, grid may only be given if 'x' is given")
    Call(#norm=norm,
         M=M, raw=raw, var=var, proj=proj)
  } else {
    PL <- getRFoptions(getoptions_="basic")$printlevel
    .Call(C_setlocalRFutils, NULL, 0)
    new <- C_UnifyXT(x, y, z, T, grid)
    .Call(C_setlocalRFutils, NULL, PL)
    Call(#norm=norm,
         M=M, x=new, raw=raw, var=var, proj=proj)
  }
}
RMfixcov <- copyProp(RMfixcov, iRMfixcov)


RMcov <- function(gamma, x, y=NULL, z=NULL, T=NULL, grid, a,
                  var, scale, Aniso, proj, raw) {
  Call <- iRMcov
  if (missing(x) && length(T)==0) {
    if (length(y)!=0 || length(T)!=0 || !missing(grid))
      stop("y, z, T, grid may only be given if 'x' is given")
    x <- list(1) ## origin
  } else if (is.character(x)) {
    if (!is.null(y) || !is.null(z) || !is.null(T) || !missing(grid))
      stop("If x is a character, the arguments y,z,T,grid may not be given.")
    x <- pmatch(x, RMCOV_X)
    if (is.na(x)) stop("unknown choice of 'x'")
    x <- list(as.double(x))
  } else {
    PL <- getRFoptions(getoptions_="basic")$printlevel
    .Call(C_setlocalRFutils, NULL, 0)
    x <- C_UnifyXT(x, y, z, T, grid)
    .Call(C_setlocalRFutils, NULL, PL)
  }
  Call(gamma=gamma, x = x, a=a, scale=scale, Aniso=Aniso, proj=proj, var=var)
}
RMcov <- copyProp(RMcov, iRMcov)


RMpolynome <- function(degree, dim, value=NA,
		       coordnames = c("x", "y", "z", "T"),
                       proj=1:4) {
  if (hasArg("proj")) {
    if (!is.numeric(proj)) stop("'proj' must be a vector of integers")
    if (hasArg("dim")) stopifnot(dim == length(proj))
    else dim <- length(proj)
  } else proj <- 1:dim
  if (!hasArg("coordnames")) coordnames <- coordnames[proj]
  if (degree < 0 || degree > 5) stop("the degree is out of range")
  if (dim < 0  || dim > 4) stop("the dimension is out of range")
  x <- as.matrix(do.call("expand.grid", rep(list(0:degree), dim)))
  sums <- rowSums(x)
  y <- NULL
  for (i in 0:degree) {
    idx <- sums == i
    y <- rbind(y, x[idx, ])
  }
  n <- nrow(y)
#  y <- as.vector(y)
  z <- paste(rep(paste(" ", coordnames, sep=""), each=n),
             ifelse(y>1, "^", ""),
             ifelse(y>1, y, ""), sep="")
  z[ y == 0] <- ""
  dim(z) <- dim(y)
  z <- apply(z, 1, paste, collapse="", sep="")
  m <- length(z) - length(value)
  if (m > 0) value <- c(value, rep(NA, m)) else
  if (m < 0) value <- value[1:length(z)]
  cat( paste( value, z, collapse = " + ", sep=""), "\n" )
  
  z <- paste(rep(paste("R.p(", proj, ")", sep=""), each=n),
             ifelse(y>1, "^", ""),
             ifelse(y>1, y, ""), sep="")
  z[ y == 0 ] <- ""
  if (length(z) > 100) stop("maximum is ", MAXSUB, "^2 terms")
  dim(z) <- dim(y)
  z <- apply(z, 1, function(x)  paste(x[x!=""] , collapse="*", sep=""))
  value <- paste(R_CONST, "(", value, ")", sep="")
  z <- paste(value , z, sep="*")
  z[1] <- value[1]
  idx <- as.integer(length(z) / MAXSUB) * MAXSUB  
  if (idx > 0) {
    zz <- z[ 1:idx ]
    dim(zz) <- c(MAXSUB, length(zz) / MAXSUB)
    zz <- apply(zz, 2, function(x) paste("RMplus(", paste(x, collapse=", "), ")" ))
  } else zz <- NULL
  if (idx < length(z)) {
#    Print(zz, idx, idx + 1 == length(z))
    zz[length(zz) + 1] <-  if (idx + 1 == length(z)) z[length(z)] else 
       paste("RMplus(", paste(z[(idx+1) : length(z)], collapse=", "), ")" )
      
#    Print("A", idx, zz, (idx + 1) : length(z))
  }

  if (length(zz) > 1)
    zz <- paste("RMplus(", paste(zz, collapse=", "), ")")

 #  Print( zz)
  ##invisible
  return(eval(parse(text = zz)))
}
RMpolynome <- copyProp(RMpolynome, iR.c)


R.c <- function(a, b, c, d, e, f, g, h, i, j, l, m, n, o, p, q, ncol, factor) {
  if (length(a) == 1)
    iR.c(a=a, b=b, c=c, d=d, e=e, f=f, g=f, h=h, i=i, j=j, l=l,
         m=m, n=n, o=o, p=p, q=q, ncol=ncol, factor=factor)
  else {
    if(hasArg("b"))
      stop("several arguments may be given only when all of them are scalars")
    x <- as.list(a)
    names(x) <- letters[1:length(x)]
    if (hasArg("ncol")) x$ncol <- ncol
    if (hasArg("factor")) x$factor <- factor
    do.call("iR.c", x)
  }
}
RMpolynome <- copyProp(RMpolynome, iR.c)


xRMranef <- function(formula=NULL, data, x, y=NULL, z=NULL, T=NULL, grid,
                    var, scale, Aniso, proj, raw, norm) {
  if (hasArg("data") || is.numeric(formula) || is.data.frame(data) || 
      is(x, "RFsp") || isSpObj(x)) {
   formula <- RMcovariate(formula=formula, data=data, x=x, y=y, z=z, T=T,
                           grid=grid, raw=raw, norm=norm, addNA = TRUE)
  } else {
    if (!isS4(formula) || !is(formula, CLASS_CLIST))
      stop("'formula' is not a 'RMmodel' as expected")
     if (!missing(data) || !missing(x) || !is.null(y) ||
         !is.null(z) || !missing(grid) || !missing(raw) || !missing(norm))
       stop("If 'formula' is an 'RMmodel' then only 'var', 'scale', 'Aniso', and 'proj' might be given") 
  }
  if (missing(var)) {
    if (getRFoptions(getoptions_="basic")$printlevel > 0)
      message("Note that if 'var' is not given in 'RMranef', 'var' is set to 'NA' i.e., the variance is estimated'.")
    var <- NA
  }
  #RMraneffct(formula, var, scale, Aniso, proj)
}

  
XXXRMprod <- function(phi, var, scale, Aniso, proj) {
  #RMraneffct(phi, var, scale, Aniso, proj)
}

RMtrendplus <- function(C0, C1, C2, C3, C4, C5, C6, C7, C8, C9,
                        add.na=FALSE, var, scale, Aniso, proj) {
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('C0') && !is.null(subst <- substitute(C0))) 
    par.model[['C0']] <- CheckArg(C0, subst, TRUE)
  if (hasArg('C1') && !is.null(subst <- substitute(C1))) 
    par.model[['C1']] <- CheckArg(C1, subst, TRUE)
  if (hasArg('C2') && !is.null(subst <- substitute(C2))) 
    par.model[['C2']] <- CheckArg(C2, subst, TRUE)
  if (hasArg('C3') && !is.null(subst <- substitute(C3))) 
    par.model[['C3']] <- CheckArg(C3, subst, TRUE)
  if (hasArg('C4') && !is.null(subst <- substitute(C4))) 
    par.model[['C4']] <- CheckArg(C4, subst, TRUE)
  if (hasArg('C5') && !is.null(subst <- substitute(C5))) 
    par.model[['C5']] <- CheckArg(C5, subst, TRUE)
  if (hasArg('C6') && !is.null(subst <- substitute(C6))) 
    par.model[['C6']] <- CheckArg(C6, subst, TRUE)
  if (hasArg('C7') && !is.null(subst <- substitute(C7))) 
    par.model[['C7']] <- CheckArg(C7, subst, TRUE)
  if (hasArg('C8') && !is.null(subst <- substitute(C8))) 
    par.model[['C8']] <- CheckArg(C8, subst, TRUE)
   if (hasArg('C9') && !is.null(subst <- substitute(C9))) 
     par.model[['C9']] <- CheckArg(C9, subst, TRUE)
  par.model[[ADD_NA]] <- CheckArg(add.na, NULL, TRUE)
  par.model$trend <- TRUE
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
    par.general[['var']] <- CheckArg(var, subst, TRUE)
  if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
    par.general[['scale']] <- CheckArg(scale, subst, TRUE)
  if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
    par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
  if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
    par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new(CLASS_CLIST,
                        name = 'RMplus',  ## !!
                        submodels = submodels, 
                        par.model = par.model,
                        par.general = par.general)
  return(model)
}


RMtrendplus <- new(CLASS_RM, 
	.Data = RMtrendplus,
	type = c('trend'),
	isotropy = c('submodel dependent'),
	domain = c('submodel dependent'),
	operator = TRUE,
	monotone = 'submodel dependent monotonicity',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = SUBMODEL_DEP,
	vdim = SUBMODEL_DEP
	)


RMmodelplus <- function(C0, C1, C2, C3, C4, C5, C6, C7, C8, C9,
                        trend, var, scale, Aniso, proj) {
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('C0') && !is.null(subst <- substitute(C0))) 
    par.model[['C0']] <- CheckArg(C0, subst, TRUE)
  if (hasArg('C1') && !is.null(subst <- substitute(C1))) 
    par.model[['C1']] <- CheckArg(C1, subst, TRUE)
  if (hasArg('C2') && !is.null(subst <- substitute(C2))) 
    par.model[['C2']] <- CheckArg(C2, subst, TRUE)
  if (hasArg('C3') && !is.null(subst <- substitute(C3))) 
    par.model[['C3']] <- CheckArg(C3, subst, TRUE)
  if (hasArg('C4') && !is.null(subst <- substitute(C4))) 
    par.model[['C4']] <- CheckArg(C4, subst, TRUE)
  if (hasArg('C5') && !is.null(subst <- substitute(C5))) 
    par.model[['C5']] <- CheckArg(C5, subst, TRUE)
  if (hasArg('C6') && !is.null(subst <- substitute(C6))) 
    par.model[['C6']] <- CheckArg(C6, subst, TRUE)
  if (hasArg('C7') && !is.null(subst <- substitute(C7))) 
    par.model[['C7']] <- CheckArg(C7, subst, TRUE)
  if (hasArg('C8') && !is.null(subst <- substitute(C8))) 
    par.model[['C8']] <- CheckArg(C8, subst, TRUE)
   if (hasArg('C9') && !is.null(subst <- substitute(C9))) 
     par.model[['C9']] <- CheckArg(C9, subst, TRUE)
  if (hasArg('trend') && !is.null(subst <- substitute(trend))) 
     par.model[['trend']] <- CheckArg(trend, subst, TRUE)
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
    par.general[['var']] <- CheckArg(var, subst, TRUE)
  if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
    par.general[['scale']] <- CheckArg(scale, subst, TRUE)
  if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
    par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
  if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
    par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new(CLASS_CLIST, 
                        name = 'RMplus',  ## !!
                        submodels = submodels, 
                        par.model = par.model,
                        par.general = par.general)
  return(model)
}


RMmodelplus <- new(CLASS_RM, 
	.Data = RMmodelplus,
	type = c('of manifold type'),
	isotropy = c('submodel dependent'),
	domain = c('submodel dependent'),
	operator = TRUE,
	monotone = 'submodel dependent monotonicity',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = SUBMODEL_DEP,
	vdim = SUBMODEL_DEP
	)

