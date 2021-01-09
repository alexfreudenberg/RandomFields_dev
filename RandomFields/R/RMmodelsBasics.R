

CheckError <- function(name) {
  rarg <- paste0("r", name)
  if (exists(name)) {
    if (is(name, CLASS_CLIST))
      stop("'", name, "' may not be a submodel")
    else stop("'", name, "' is not of correct type")
  }
  if (exists(rarg) && is.function(get(rarg)))
    stop("'", name, "' may not be a random parameter")
  stop("Neither '", name, "' nor the distribution family 'p", name, "', 'd", name, "', 'q", name, "', 'r", name, "' is undefined.")
}

CheckMixed <- function(arg, subst, names) {
  ## never allows for distributions
  if (is.character(arg)) {
    if (length(arg) != 1) stop("'", deparse(substitute(arg)),
			       "' must be a single string")      
    arg <- -pmatch(toupper(arg), toupper(names))
    if (is.na(arg)) stop("'", deparse(substitute(arg)),
			 "': unknown character sequence")
    return(arg)
  }
  u <- rawTry(is.numeric(arg) || is.logical(arg) || is.language(arg)
              || is.list(arg) || (isS4(arg) && is(arg, CLASS_CLIST)))
  if (is.logical(u) && u) arg
  else CheckError(deparse(substitute(arg)))
}



CheckProj <- function(arg, subst) {
  if (is.character(arg)) {
    if (length(arg) != 1) stop("'proj' must be a single string")
    ## anyway most of the following part is written as if p was vector
    p <- -pmatch(arg, PROJECTION_NAMES)
    if (any(is.na(p))) {
      p <- pmatch(arg, RFoptions()$coord$coordnames)
      if (length(p) != length(arg))
	stop("projection definition '", arg,  "' is odd")
      if (any(is.na(p))) { 
	p <- pmatch(arg, c("x", "y", "z", "t", "X", "Y", "Z", "T"))
	if (!any(is.na(p))) { 
	  p <- (p - 1) %% 4 + 1
	  p[p == 4] <-
	  pmatch("time", PROJECTION_NAMES) -1 -length(PROJECTION_NAMES)
	} else {
	  p <- arg ## must be a single character string
	  while (!(substr(p, 1, 1) %in% 0:9)) p <- substring(p, 2)
	  warn <- options()$warn
	  options(warn = 0)            
	  p <- try(as.vector(p), silent=TRUE) ## to fool findall: tkEntry(
	  options(warn=warn)
	  if (is(p, "try-error"))
	    stop("'\"", arg,
		 "\"' is interpretated as a projection defintion but the character string is not recognized")
	}
        }
    }
    return(p)
  }
  u <- rawTry(is.numeric(arg) || is.logical(arg) || is.language(arg)
               || is.list(arg) || (isS4(arg) && is(arg, CLASS_CLIST)))
  if (is.logical(u) && u) return(arg)
  else if (substr(deparse(subst), 1, 1)=='R') arg
  else  do.call('RRdistr', list(subst))
}



CheckMaths <- function(arg, subst, distr) {
  u <- rawTry(is.numeric(arg) || is.logical(arg) || is.language(arg)
              || is.list(arg) || (isS4(arg) && is(arg, CLASS_CLIST)))
  if (is.logical(u) && u) arg
  else if (is.character(arg)) do.call('R.p', list(arg))
  else {
    name <- deparse(subst)
    if (substr(name, 1, 1)=='R') arg
    else if (distr) do.call('RRdistr', list(subst))
    else CheckError(name)
  }
}


CheckArg <- function(arg, subst, distr) {
  ## See also CheckArg in PrepareModel
  ##  print(deparse(subst))
  # str(arg)
  # Print("CheckArg", deparse(subst), is(arg, "formula"))
   u <- rawTry(is.numeric(arg) || is.logical(arg) || is.language(arg)
               || is.list(arg) || (isS4(arg) && is(arg, CLASS_CLIST)))
  if (is.logical(u) && u) arg
  else {
    name <- deparse(subst)
    if (substr(name, 1, 1)=='R') arg
    else if (distr) do.call('RRdistr', list(subst))
    else CheckError(name)
  }
}

CheckChar <- function(arg, subst, names=NULL, distr) {
  if (is.character(arg)) {
    if (length(arg) != 1) stop("'", deparse(substitute(arg)),
			       "' must be a single string")
    if (length(names) > 0) {
      arg <- pmatch(arg, names) - 1
      if (any(is.na(arg))) stop("'", deparse(substitute(arg)),
                                "': unknown string value")
    }
    return(arg)
  }
  u <- rawTry(is.numeric(arg) || is.logical(arg) || is.language(arg)
              || is.list(arg) || (isS4(arg) && is(arg, CLASS_CLIST)))
  if (is.logical(u) && u) arg
  else {
    name <- deparse(subst)
    if (substr(name, 1, 1)=='R') arg
    else if (distr) do.call('RRdistr', list(subst))
    else CheckError(name)
  }
}


copyProp <- function(what, from) {
  return(new(CLASS_RM,
             .Data = what,
             type = from["type"],
             isotropy = from["isotropy"],
             domain = from["domain"],
             operator = from["operator"],
             monotone = from["monotone"],
             finiterange = from["finiterange"],
             simpleArguments = from["simpleArguments"],
             maxdim = from["maxdim"],
             vdim = from["vdim"]
             ))
}
