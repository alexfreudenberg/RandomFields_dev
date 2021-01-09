
## Authors 
## Martin Schlather, schlather@math.uni-mannheim.de
##
### Diese Datei wandelt RMmodel in eine Liste um
##
## Copyright (C) 2015 -- 2016 Alexander Malinowski, Martin Schlather,
##                            ctr=Sebastian Gross
##               2017 -- 2019 Martin Schlather
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

#Zaehler <- 0


## NOTE:
## "definitions" in the man page is called here "symbol"


AddTransform <- function(M, model, x, params, orig.params, Names, nNA,
                         revOrdering,
                         symbols, formulae, genuine.formulae, 
                         general.coord, factor,
                         xdim, data, MODEL_AUX, add.na, Env) {
  ## Print("Entering AddTrafo")
  
  ## not that the number of formulae could have been reduced, due to
  ## constants given through formulae by the user      
  paramsNA <- params ## must be at the very beginning -- params will change
  n.params <- length(params)
  cumNA <- cumsum(c(0, nNA))
  totalNA <- cumNA[length(cumNA)]
  values <- as.double(100001 * MAXACCEPTED * (1:MAX_NA))
  for (i in 1:n.params) { ## replace NAs/NaNs by dummies to
    ## identify them on the C level
    ## fctn values all set to NaN
    params[[i]][if (any(is.nan(params[[i]]))) TRUE
                else is.na(params[[i]])] <- values[(cumNA[i] +1) : cumNA[i+1]]
    assign(Names[i], params[[i]], envir=Env)
  }
  
  M300 <- parseModel(model=model, x=x, Env=Env,
                     EnvDummies = new.env(parent=Env), add.na=add.na,
                     idx.coord = rep(FALSE, length(M$data.coordnames)),
                     coord.names=M$data.coordnames, data=data,
                     ##                         return_transform = FALSE,
                     general.coord=general.coord, factor=factor)$model
  
  ## prepare model to determine the indexing between R and C
  ## of the parameters that are NA
  ## rf_interface.cc, includeparam, these values cannot appear
  ## as parameter values
  
  coords <- if (!is.null(M$C_coords)) M$C_coords
            else {
              if (xdim==0 && !missing(data) && !is.null(data))
                xdim <- max(1, sum(M$idx.coord))
              if (xdim > 0) {                    
                xx <- matrix(1, nrow=3, ncol=xdim)
                  xx[3] = nrow(data)
                C_UnifyXT(xx, grid=TRUE)
              }
              else stop("neither 'dim' not 'x' given")
            }
  
  model.vdim <- integer(1)
  NAs.in.model <- ## user might have given NA directly in the model
    .Call(C_GetNAPositions, MODEL_AUX, list("Covariance", M300), coords,
          as.double(NA),
          if (totalNA == 0) values else values[-1:-totalNA],
          FALSE, model.vdim, 0L)
  tot.NAs.in.model <- sum(NAs.in.model)
  cumNA <- cumsum(c(0, c(nNA, NAs.in.model)))
  
  
  ## legend for transforms below:
  ## C : sequence in C running through all parameter 1..length(CtoP)<=totalNA
  ## R : 1..totalNA (sequence of all unknown values given by param)
  ## P : 1..length(param) (numbering of the params)
  ## F : 1..#formulae  (numbering of the formulae within params)
  ## G : 1..#genuine.formulae
  ## S : 1..#symbols; numbering of !formulae; #!formu+#symb=length(param)
  ## AmidB : A restricted to B, numbering starting from 1    
  ## note: input parameter are in order of C for the non formulae
  
  values <- values[OneTo(totalNA + tot.NAs.in.model)]
  CtoR <- .Call(C_GetNAPositions, MODEL_AUX, NULL, coords,
                values, double(0), FALSE, model.vdim, # side effect !!
                0L)
  n.CtoR <- length(CtoR)
  if (n.CtoR > 0) {
    Check.CtoR(CtoR, MODEL_AUX, values)
    RtoC <- rep(NA, (totalNA + tot.NAs.in.model))
    RtoC[CtoR] <- 1:n.CtoR

    genuine.symbols <- symbols & (nNA > 0)
    S2P <- which(genuine.symbols)
    S2Pxx <- which(c(genuine.symbols, rep(TRUE, length(NAs.in.model))))
    S2C <- cbind(cumNA[-length(cumNA)][S2Pxx] + 1, cumNA[-1][S2Pxx])
    S2C <- unlist(apply(S2C, 1, function(x) list(RtoC[x[1]:x[2]])),
                  recursive=FALSE)
    im.S2C <- unlist(S2C) # im = (flattened) image
    if (any(is.na(im.S2C))) {
      ## Print(RtoC, CtoR, S2C, im.S2C,NAs.in.model,RFgetModelInfo(MODEL_AUX))
      stop("incorrect 'params'. Probably not all newly introduced variables have been defined in 'RMdeclare'.")
    }
    CmidS <- integer(n.CtoR)

    
    CmidS[im.S2C] <- rank(im.S2C)
    ## S2CmidS is the sequence, by which input variables come in:
    ## (input variables have the same ordering than the one
    ## inside the ordering of all NAs and NaNs)
    S2CmidS <- lapply(S2C[OneTo(length(S2P))], function(x) CmidS[x])

##    Print(S2CmidS, cumNA, S2P, S2Pxx, S2C, CmidS, im.S2C, rank(im.S2C), CtoR, RtoC, M300, RFgetModelInfo(MODEL_AUX))
        
    G2P <- which(genuine.formulae)
    n.genuine.formulae <- length(G2P)
    G2C <- cbind(cumNA[-length(cumNA)][G2P] + 1, cumNA[-1][G2P])
    G2C <- unlist(apply(G2C, 1, function(x) list(RtoC[x[1]:x[2]])),
                  recursive=FALSE) ## to be save that it always a list
    used <- sapply(G2C, function(x) { finite <- !any(is.na(x));
      ##      Print(G2C, x, finite)
      if (!finite && !all(is.na(x))) stop("programming error in convert.R");
      finite})
    im.G2C <- unlist(G2C[used])
    
    AnsName <- "..bca.."
    if (AnsName %in% Names)
      stop("'", AnsName, "' is not allowed as parameter name")
    
    p <- A <- character(length(S2P))
    paramS2P <- paramsNA[S2P]
    namesS2P <- Names[S2P]
    subset <- which(sapply(paramS2P, length) > nNA[S2P])
    paramSubset <- paramS2P[subset]
    ## 'paste' macht bei Liste von Aufzaehlungen das richtige!!
    p[subset] <- paste0("<- ",  paramSubset, "\n\t", 
                        namesS2P[subset], "[",
                        lapply(paramSubset, function(x) which(is.na(x))),
                        "]")

    subset <- which(!sapply(paramS2P, is.vector))
    A[subset] <- paste0("\n\tdim(", namesS2P[subset], ") <-",
                        lapply(paramS2P[subset], dim))
    
    mock <- which(formulae & !genuine.formulae)
    
    G2text <- character(n.genuine.formulae) # for each formula a text
    for (i in OneTo(n.genuine.formulae))
      G2text[i] <- as.character(as.expression(orig.params[[ G2P[i] ]][[2]] ))
    ordering <- order(revOrdering[G2P])
    
    text <- paste0(sep="", "function(variab) {\n\t",
                   if (length(S2P) > 0) ## otherwise no variable is given
                     paste(namesS2P, p, "<- variab[", S2CmidS, "]",
                           A, "\n\t", collapse="")
                  ,
                   if (length(mock) > 0)
                     paste(Names[mock], "<-", params[mock],"\n\t",collapse="")
                  ,
                   if (length(G2P) > 0)
                     paste(Names[G2P][ordering], "<-", G2text[ordering],
                           "\n\t", collapse="")
                  ,
                   AnsName, " <- double(", n.CtoR, ")\n\t"
                  ,
                   if (length(im.S2C) > 0) 
                     paste0(AnsName, "[",list(sort(im.S2C)),"] <- variab\n\t")
                  ,
                   if (any(used))#list(Names) hier nicht moeglich; paste notw!
                     paste0(AnsName, "[", list(im.G2C), "] <- c(",
                            paste(Names[G2P[used]], collapse=","), ")\n\t")
                  ,
                  AnsName, "\n}\n")
    
    ##     cat(text)
##    Print(Names[G2P[used]], Names, G2P[used], G2P, used)
    
    fctn <- eval(parse(text=text), envir=NULL)
    fctn.env <- new.env(parent=baseenv())
    environment(fctn) <- fctn.env
    dots <- get("..dots..", env=Env)
    assign("..dots..", dots, env=fctn.env)
    if (length(dots) > 0)
      eval(parse(text = paste(names(dots),
                              "<- ..dots..[[", 1:length(dots), "]]",
                              collapse=";")),
           envir=fctn.env)

    ## Print(ls(envir = fctn.env), "OK")
    
    ## isNA: logical, Laenge (CtoR); true falls symbol
    isNA <- rep(FALSE, n.CtoR)
    isNA[im.S2C] <- TRUE
    
    CmidN2RmidN <- .Call(C_GetNAPositions, MODEL_AUX, NULL, coords,
                         values[-which(rep(formulae, nNA))], double(0),
                         FALSE, model.vdim, 0L)
    if (length(CmidN2RmidN) != sum(nNA[S2P]) + tot.NAs.in.model)
      stop("number of parameters declared to be estimated differs from the number of non-linear parameters to be estimated. In other words, for any declaration of the form 'xy = NA' in 'params' needs exactly one counter part within in the (RM)model defintion of the form 'arg = xy'")
  } # at least 1 variable to be filled on C level
  
  M$transform <- list(isNA=isNA, fctn)
  M$vdim <- model.vdim
  return(M)
}


IsFactor <- function(data, data.names, select) {
  if (is(data, "factor")) TRUE
  else if (is.data.frame(data)) {
    f <- sapply(data, is.factor)
    if (any(f & !select))
      stop("dependent variable (",
           paste(data.names[f & !select], collapse=","),
           ") seems to be a factor")
    f[select]
  } else if (is.matrix(data) || is(data,"RFsp"))
    rep(FALSE, sum(select))
  else if (is.vector(data)) FALSE
  else {
    Print(data) ## ok
    stop("unknown data format")
                ## simulation from randomFields -- unclear, what / how should
    ## any relevant information be extracted  TO DO???
  }
}


RemoveBlanks <- function(u) sub("^ *", "", sub(" *$", "", u))
Check.defn <- function(rightSide, Names, nNA, symbols) {
  d <- strsplit(rightSide, "RMdeclare")[[1]]
  ## RFdeclarations must be treated separately:
  if (length(d) == 2) {        
    cs <- strsplit(d[2], "")[[1]] ## split into single characters
    cs <- cumsum((cs == "(") - (cs == ")"))
    end <- which(cs == 0)[1] - 1
    s <- substr(d[2], which(cs == 1)[1]+1, ## next character after first "("
                end) ## last character before closing ")"
###         of RMdeclare (and not anything that has the user defined inbetween)
    zd <- strsplit(strsplit(s, ",")[[1]], "=")
    single <- sapply(zd, length) == 1  ## no "="
    nd <- RemoveBlanks(zd[single])
    neq <- matrix(RemoveBlanks(unlist(zd[!single])), nrow = 2)
    if (any(idx <- neq[1, ] != neq[2, ]))
      stop(paste0("In 'RMdeclare', in the declaration of ",
                  paste0("'", neq[1, idx], "'", collapse=","),
                  ", the left hand side does",
                  " not exactly match the right hand side."))
##    Print(cs, s, zd, nd, zd, d[2]) 
    rightSide <- paste(d[1], substring(d[2], end + 2)) # end + ")".Then next
  } else nd <- neq <- NULL
  
  s <- strsplit(rightSide, "(", fixed=TRUE)[[1]] #split at "(", ",", ")" 
  s <- unlist(strsplit(unlist(strsplit(s, ",")), ")")) 
  z <- strsplit(s, "=")
  nm <- sapply(z[sapply(z, length) == 2], ## count frequ. of right side
               ## of "=", ignoring blanks at beginning and of the words
               function(u) RemoveBlanks(u[2]))
  
  t <- table(c(nm, neq[1,], nd))
  nn <- Names[symbols & nNA > 0]
  f <- t[nn]
  f[is.na(f)] <- 0
  names(f) <- nn
  nn2 <- nn[f >= 2]             
   ## 
  if (any(f != 1)) {
    # Print(nm, neq, nd, nn2, nn, t, f,symbols & nNA > 0, rightSide)
    stop("We call an expression in 'params' that includes NAs on the right hand side a definition. (In particular, a definition is not given by a formula.) The name of each definition must appear exactly once in the RMmodel definition as a \"value\" of an argument. This is not the case for ",
         paste0("'", nn[f != 1], "'", collapse=","),
         ". Likely, you should do the following:\n", 
         if (any(f == 0)) {
           if (length(d) == 2)
             paste("* include", paste0("'", nn[f == 0], "'", collapse =","),
                   "in 'RMdeclare'.\n")
           else paste0("* add 'RMdeclare(", paste0(nn[f==0], collapse =","),
                       ")' to the model.\n")
         },
         if (length(nn2) > 0) {
           zd <- sapply(zd, function(x) x[1])
           if (any(idx <- nn2 %in% zd))
             paste0("* remove ",
                    paste0("'", nn2[idx], "'", collapse=", "),
                    " from 'RMdeclare'.\n")
           else {
             paste0("* introduce additional formulae, e.g.,",
                    paste0("'new.", nn2, " = ~", nn2, "'",
                           collapse=", "),
                    ", in 'params'. Replace the second appearance of ",
                    paste0("'", nn2, "'",collapse=","), " by ",
                    paste0("'new.", nn2, "'",collapse=","),
                    if (length(nn2) > 1) ", respectively", ".",
                    if (any(f > 2))
                      " Replace further appearances also by further formulae.",
                    "\n")
           }
         },
         "\n")
  }      
}


Check.CtoR <- function(indices, MODEL_AUX, values) {
  if (length(indices) <= 1)
    stop("Got a wrong index set.", CONTACT)  
  if (anyDuplicated(indices)) {
        
    GetParam <- function(p, name) {
  #    Print(p, name, v)
      N <- names(p)
      ##         Print(p, N, name, v,length(p))         
      if (is.numeric(p)) return(if (any(p == v)) list(name) else NULL)
      if (!is.list(p)) return(NULL)
      model <- if (p$name %in% c(RM_PLUS, RM_MULT)) NULL else p$name
      d <- NULL
      if (length(p$submodel) > 0) {
        n0 <- c(name, model)
        for (j in 1:length(p$submodels))
          d <- c(d, GetParam(p$submodels[[j]], n0))                
      }
      if (length(p$param) > 0) {
        if (p$name %in% DOLLAR) name <- c(name, p$submodel[[1]]$name)
        p.name <- names(p$param)
        for (j in 1:length(p$param)) {
          d <- c(d, GetParam(p$param[[j]], c(name, p.name[j]))) 
        }
      }
      return(d)
    }

    v <- values[indices[duplicated(indices)]]
    link <- GetParam(RFgetModelInfo(MODEL_AUX), NULL)
    if (length(link) < 2) stop("Severe error.", CONTACT)
    stop("There is a forbidden link between ",            
         paste("'", sapply(link[1:2], function(x) paste(x,collapse="..")),
               "'", collapse=" and ", sep=""),
         ". Are all relevant parameters given by formulae in 'params'?")
  }
}

PrepareModel2 <- function(model, ...,
                          params=NULL, ## in case of formulae
                          further.model = NULL, ##  may depend on the same
                          ##  parameter list and be a torsus (e.g. lower/upper)
                          arg=NULL, # takes the values "lower" and "upper"
                          ##   and is used to check parameter validity between
                          ## further.model and model
                          x=NULL, # list used within buildFactorList to
                          ## formulate RM_COVARIATE model and to prepare the
                          ## transform
                          trend=NULL, ## a model indicating a trend -- a further
                          ## option to define a trend; other possiblities are
                          ## trend=TRUE in '+' or to hope for automatic
                          ## detection
                          data = NULL, # data, needed here to detect on the fly
                          ## which columns are used are variables or coordinates
                          coord.opt = NULL, # these are the RFoptions for
                          ## the coordinate names, also for column detection
                          add.na=FALSE, # for trends: should NAs be added?
                          ## THis is needed for estimation in linear models
                          return_transform = TRUE,
                          xdim = 0,
                          Env = new.env(parent=.GlobalEnv), ## nur im
                          ##                       rekursiven Aufruf verwendet
                          EnvDummies = new.env(parent=Env)
                          ) {


 # stopifnot(Zaehler <- Zaehler + 1 < 2)
  ## (1) rffit->unifydata
  ## (2)
  ## (3)
  ## (4)
  
  default.max.number.data <- 99
  ##  Print("Prep")#, data)

  
  dots <- list(...)
  assign("..dots..", dots, envir=Env)
  if (length(dots)>0) {
    if (any(0 != sapply(dots, function(x) length(oldClass(x)))))
      stop("arguments in '...' may not have a class")
    eval(parse(text = paste(names(dots), "<- ..dots..[[", 1:length(dots), "]]",
                            collapse=";")),
         envir=Env)
  }
  
  if (!missing(model)) {
    RFopt <- RFoptions(GETOPTIONS=c("basic", "internal"))
    if (is(model, CLASS_FITLIST)) {
      if (RFopt$basic$helpinfo && RFopt$internal$note_mle) {
        RFoptions(internal.note_mle = FALSE)
        message("The MLE is extracted from a list of results. Use argument 'method' in any function of the package 'RandomFieldsLight' to extract different results, see ?RFfit for possible values of 'method'. Set RFoptions(helpinfo=FALSE) to avoid this message.")
      }
      model <- if (isS4(model)) model["ml"] else model[["ml"]]
    }
    if (is(model, CLASS_SINGLEFIT))
      model <- if (isS4(model)) model@model else model$model
    if (is(model, CLASS_RM))
      model <- do.call(model, list())
    if (is(model, CLASS_CLIST) && length(model@name) == 0)
      stop("model is empty")
  }
 
  if (length(trend) > 0) {
    ## trend might be given separately from the rest of the model
    ## this is a semi-implemented option for the user and particularly
    ## used by RandomFieldsLight
    if (isFormula  <- is(trend, "formula")) {
      Trend <- as.character(trend)
      Tv <-if (length(Trend) == 3) Trend[2]
      trend <- paste("RMtrendplus(add.na=", add.na, ", ~",
                     Trend[length(Trend)], ")")
    }

    if (missing(model)) {
      model <- if (isFormula) eval(parse(text=paste(Tv, "~", trend)))
               else RMtrendplus(trend, add.na=add.na)
    } else {
      if (is(model, "formula")) {
        m <- as.character(model)
        v <- if (length(m) == 3) m[2] else if (isFormula) Tv
        model <- paste(v,"~RMplus(RMmodelplus(trend=FALSE, ",
                       m[length(m)],"), ")
        if (isFormula) {
          if (length(Trend) == 3 && length(m) == 3 && !all(v == Tv))
            stop("dependent variables in the covariance model (",
                 paste(v, collapse=", "),
                 ") differ from those in the mean model (",
                 paste(Tv, collapse=", "), ".")
          model <- eval(parse(text=paste(model, trend, ")")))
        } else {
          assign("..trend..", trend, envir=Env)
          assign("..trend..", trend, envir=EnvDummies)
          model <- eval(parse(text=paste(model,"RMtrendplus(add.na=",
                                         add.na, ", ..trend..))")))
        }
      }
      
      else if (isFormula) {
        assign("..model..", model, envir=Env)
        assign("..model..", model, envir=EnvDummies)  
        model  <- eval(parse(text=paste(Tv, "~RMplus(..model.., ", trend, ")")))
        model  <- eval(parse(text=paste(Tv, "~RMplus(", trend, ",..model..)")))
        
      }
      
      else model <- RMplus(model, 
                           RMtrendplus(trend, add.na=add.na),
                           RMmodelplus(trend = FALSE, model))
    }
    
    return(PrepareModel2(model, params=params, 
                         further.model = further.model, arg=arg, x=x,
                         data=data,
                         coord.opt = coord.opt,
                         add.na = FALSE,
                         return_transform = return_transform,
                         xdim = xdim,
                         Env = Env,
                         EnvDummies = EnvDummies
                         ))
  }
  
  if (missing(model) || is.null(model)) stop("'model' must be given.")

  factor <-  NULL
  genuine.formulae <- formulae <- FALSE
  values <- double(0)  
 # warn <- options(warn=2)
#  options(warn 
  if (!missing(params) && length(params) > 0) {
    ## simple general stuff for params; reset params in case of further.model(s)
    ## latter are the lower bound model or the uppper bound model
    orig.params <- params
    n.params <- length(params)
    Names <- names(params)
    if (!all(nchar(Names) > 0)) stop("parameters are not all given by name")
    if (anyDuplicated(Names))
      stop("names on the left sides in 'params' must be unique.")
    formulae <- logical(n.params) ## FALSE
    for (i in 1:n.params) formulae[i] <- is(params[[i]], "formula")
    symbols <- !formulae

    if (!missing(further.model) && !is.null(further.model)) {
      if (is.list(further.model) && is.character(further.model[[1]]))
        stop("Do not mix different styles of model definition.\n  See ?RFformula and ?RFformulaAdvanced for details.")
      FurtherModels <- "In bounds, initial values and parscales, "
      for (n in names(further.model)) {
        idx <- which(n == Names)
        #Print(idx, n, Names)
        if (length(idx) > 0) {
          if (length(idx) != 1)
            stop("more than one formala given for '",  Names[idx], "'.")
          if (formulae[idx])
            stop(FurtherModels,
                 "formulae and values may not be given for\n",
                 "  dependent variables, here '",
                 Names[idx], "', since they do not make sense logically.")
          if (sum(is.na(params[[n]])) == length(further.model[[n]])) {
            params[[n]][is.na(params[[n]])] <- further.model[[n]]
          #  Print(params)
          } else {
            if (length(params[[n]]) != length(further.model[[n]]))
             stop("Specification of '", n, "' (length) in '", arg,
                   "' does not match the model definition")

            idx <- is.na(params[[n]])
            params[[n]][idx] <- further.model[[n]][idx]
            idx <- !idx & !is.na(further.model[[n]]) &
              !is.na(further.model[[n]])
            ok <-
              (if (arg == "lower")
               all(params[[n]][idx] >= further.model[[n]][idx])
              else if (arg == "upper")
               all(params[[n]][idx] <= further.model[[n]][idx])
              else 
               all(params[[n]][idx] == further.model[[n]][idx])
               )
            if (!ok) stop("Specification of '", n, "' in '", arg,
                          "' does not match the model definition.")
          }
        } else {
          assign(n, further.model[[n]], envir=Env)
          params[[n]] <- further.model[[n]]
        }
        assign(n, params[[n]], envir=Env)
      }
      if (n.params != length(formulae))
        stop(FurtherModels, "additional variables may not be defined.\n  Here ",
             paste("'", names(params)[!(names(params) %in% Names)], "'",
                   sep="", collapse=","), ".")
                                        #
      ##     Print(params); dddddd
    }

    nNA <- rep(0, n.params)
    for (i in 1:n.params) {
      if (symbols[i]) {
        assign(Names[i], params[[i]], envir=Env)
        nNA[i] <- sum(is.na(params[[i]]))
        if (any(is.nan(params[[i]])))
          stop("only NAs are allowed to indicate parameters to be estimated")
      }
    }
            
   
  } else params <- NULL
  
  n.formulae <- sum(formulae)
  if (n.formulae > 0) {
    ## identify genuine.formulae and #variables to be estimated; prepare Env
    ##
    ## params ist unveraendert orig.params, ausser !missing(further.model)
    ## letzteres schliesst wechselseitig return_transform aus!
    revOrdering <- rep(0, n.formulae) ## reverse ordering of formulae
    genuine.formulae <- formulae
    zaehler <- 0
    for (last in c(FALSE, TRUE)) {
      repeat {
        repet <- FALSE
        for (i in 1:n.params) {
          if (is(params[[i]], "formula")) {
            value  <- rawTry(eval(as.expression(params[[i]][[2]]), envir=Env))
                                        #Print(params[[i]], value)
            if (rawError(value)) { # error
              if (last) {
                ## Print(n.params, params, value, repet, formulae, genuine.formulae, symbols, nNA, ls(envir=Env))
                stop("This situation should not appear. Please check whether your definitions are sound and/or contact ", AUTHOR)
                value <- NaN 
                params[[i]] <- value
                assign(Names[i], value, envir=Env)
                nNA[i] <- 1
                revOrdering[i] <- zaehler <- zaehler + 1 ## andere Konstruktion
                ## geht nicht, da unten formeln zu nicht-Forelm erklaert werden
                }
            } else {
 #             Print(i, "", Names[i])
              revOrdering[i] <- zaehler <- zaehler + 1 ## andere Konstruktion
              ## geht nicht, da nachfolgend formeln u.U. zu nicht-Forelm
              ## erklaert werden
              genuine.formulae[i] <- any(is.na(value)) # non-trivial formulae,
              ##                      i.e. not of the form abc = ~123
              if (genuine.formulae[i]) {
                value[TRUE] <- NaN  ## there could be constants inside!
                ## We overwrite them by NaN. Otherwise there might be a mess
                ## between C-level identification of NaN positions and
                ## NaN on R level. In particular there is not check, whether
                ## the "constant" is indeed a constant!
                nNA[i] <- length(value)
              }
              params[[i]] <- value
              assign(Names[i], value, envir=Env)              
              repet <- !last
            } ## not try-error
          } ## is.formula
        } ## for i
        if (!repet) break
      } ## repeat
    } ## for last

    ## Print(formulae, genuine.formulae, params, revOrdering, ls(envir=Env), Names[formulae], Names[genuine.formulae])
  } # if any(formulae)
  ##
 
  
  data.varnames <-  data.coordnames <- data.names <-  character(0)
  idx.coord <- NULL
  ## determine names of variables on the left and on the right of ~
  ## replace ~x etc by R.p etc
  opt.given <- !missing(coord.opt) && !is.null(coord.opt)
  if (opt.given) {
    data.coordnames <- coord.opt$coordnames ## could be empty
    data.varnames <- coord.opt$varnames ## could be empty
  }

  if (length(data) > 0) {  
    ## try identify data.names, coord.names from data point of view
    if (is.list(data) && !is.data.frame(data)) data <- data[[1]]
    if (is.data.frame(data) || is.matrix(data)) data.names <- colnames(data)
    else if (is(data, "RFsp")) data.names <- colnames(data@data)
    nc <- length(data.names)
    general.coord <- paste0("X", 1:nc)

##    Print(data.names, coord.opt)
    if (nc == 0) {
      nc <- if (is(data, "RFsp")) ncol(data@data) else ncol(data)
      if (length(nc) == 0) stop("size of data cannot be determined")
      data.names <- paste0("D", 1:nc)
    } else if (any(data.names %in% general.coord)) general.coord <- NULL

    idx.coord <- rep(NA, length(data.names))
    if (length(data.coordnames) > 0) {
      p <- pmatch(data.coordnames, data.names)
      idx.coord[p[!is.na(p)]] <- TRUE
    } else if (opt.given && !is.na(coord.opt$coordidx[1])) {
      ## coord.opt$coordnames and coord.opt$coordidx are exclusive
      ende <- if (is.na(coord.opt$coordidx[2])) length(idx.coord)
              else coord.opt$coordidx[2]
      idx <- coord.opt$coordidx[1] : ende
      idx.coord[idx] <- TRUE
      data.coordnames <- data.names[idx]
    }
    if (length(data.varnames) > 0) {
      p <- pmatch(data.varnames, data.names)
      idx.coord[p[!is.na(p)]] <- FALSE
    } else if (opt.given && !is.na(coord.opt$varidx[1])) {
      ## coord.opt$coordnames and coord.opt$coordidx are exclusive
      ende <- if (is.na(coord.opt$varidx[2])) length(idx.coord)
              else coord.opt$varidx[2]
      idx <- coord.opt$varidx[1] : ende
             
      idx.coord[idx] <- FALSE
      data.varnames <- data.names[idx]
    }
  } else {
    if (length(data.coordnames) == 0) {
      data.coordnames <- paste0("X", 1:default.max.number.data) ###
      general.coord <- c("x", "y", "z", "T")
      if (length(intersect(data.varnames,
                           union(general.coord, data.coordnames))) > 0)
        stop("'coordnames' must be given.")      
    } else general.coord <- character(0)
    if (length(data.varnames) == 0) {
      data.varnames <- paste0("D", 1:default.max.number.data)
      if (length(intersect(data.varnames, data.coordnames)) > 0)
        stop("'varnames' must be given.")
    }
  }

  for (i in OneTo(length(general.coord)))
    assign(general.coord[i], general.coord[i], envir=Env) # character

  if (is(model, "formula")) {
    ## try identify data.names, coord.names from model point of view
    CM <- as.character(model)
#    Print(CM)
    if (length(CM) == 3) {
      leftSide <- CM[2]
      names <- c(data.names, data.varnames) ## important: data.names first!
      for (i in OneTo(length(names)))
        assign(data.names[i], data.names[i], envir=Env)
      ##Print(ls(envir=Env), leftSide)

      used.varnames <- eval(parse(text=leftSide), envir=Env) # indeed used ones
##      Print(data.varnames,  used.varnames)
      if (length(data) > 0)
        data.varnames <- if (length(data.varnames)==0) used.varnames
                         else intersect(used.varnames, data.varnames)
      var.idx <- match(used.varnames, names)
      if (any(is.na(var.idx))) {
        ## 30.11.20: falls !data sind neue Variablen akzeptiert??
        if (length(data.names) > 0)
          stop("unknown variable(s) on the left hand side of the formula")
        if (!all(is.na(var.idx))) ## dann muss der ganze Satz neu sein...
          stop("inconsistent use of variable names on the left side")
        if (exists("c", envir=Env))
          stop("please use another variable name instead of 'c'")
        assign("c", function(...) as.character(as.list(match.call())[-1]),
               envir=Env) ## Umdefinition von c(...)
        used.varnames <- rawTry(eval(parse(text=leftSide), envir=Env))
        if (rawError(used.varnames)) used.varnames <- leftSide
        rm("c", envir=Env)
        ## rm(list = ls(envir=Env), envir=Env) ## 30.11.20 auskommentiert
        ## dann sonst nicht symmetrisch zu else-Teil
        RFoptions(coords.varnames = used.varnames)        
      } else if (length(data) > 0) {        
        var.idx <- var.idx[var.idx <= length(data.names)]
        save.coord <- idx.coord[var.idx]
        save.coord[is.na(save.coord)] <- FALSE
        if (any(save.coord))
          stop("same names (",
               paste(data.names[var.idx[save.coord]], collapse = ","),
               ") indicate both value and coordinate. See RFoptions()$coord")
        idx.coord[var.idx] <- FALSE
        rm(list=data.names[var.idx], envir=Env)
        if (length(general.coord) > 0) ## for general.coord itself see below
         rm(list=general.coord[var.idx], envir=Env)
      }
    }

    if (length(params) > 0 && sum(nNA[symbols]) > 0)
      Check.defn(CM[length(CM)], Names, nNA, symbols)
  }

##  Print(data.names, idx.coord)
##  Print("AB",data.coordnames, data.varnames, data.names, idx.coord, data)
        
  if (length(data) > 0) {
    select <- idx.coord ## up to here NA=unclear
    select[is.na(select)] <- TRUE
    factor <- IsFactor(data, data.names, select)    
    data.names <- data.names[select]
    general.coord <- general.coord[select]
    idx.coord <- idx.coord[select]
    idx.coord[is.na(idx.coord)] <- FALSE ## from here on: FALSE=unclear
  } else {
    data.names <- data.coordnames
    idx.coord <- factor <- rep(FALSE, length(data.coordnames))
  }
    
  M <- parseModel(model=model, x=x, Env=Env,
                  EnvDummies = EnvDummies, add.na=add.na,
                  idx.coord = idx.coord,
                  coord.names=data.names,
                  general.coord=general.coord, factor=factor, data=data)
  class(M$model) <- CLASS_CLIST

  ## Print(idx.coord, data.names, factor, M$idx.coord)
  
  M$data.coordnames <- data.names[M$idx.coord & !factor]
  M$data.varnames <- data.varnames
  M$data.factornames <- data.names[M$idx.coord & factor]
  if (length(x) > 0) M$C_coords <- trafo.to.C_UnifyXT(x)
  
  if (return_transform && sum(genuine.formulae) > 0) {
    stopifnot(missing(further.model) || is.null(further.model))
    M <- AddTransform(M=M, model=model, x=x,
                      params=params, orig.params=orig.params,
                      Names=Names, 
                      nNA=nNA, revOrdering=revOrdering,
                      symbols=symbols,formulae=formulae,
                      genuine.formulae=genuine.formulae,
                      general.coord=general.coord, factor=factor,
                      xdim=xdim, data=data,
                      MODEL_AUX=MODEL_AUX,
                      add.na=add.na,
                      Env=Env)
    ##    Print(M)
 
  } else {
    model.vdim <- integer(1)
    if (length(M$C_coords) > 0) {
      .Call(C_GetNAPositions,
            MODEL_AUX, list("Covariance", M$model),
            M$C_coords, as.double(NA), as.double(NA),
            FALSE,
            model.vdim, # side effect !!
          0L)
    }
    M$vdim <- model.vdim
  }
  

  ##  print(M)

  M
}


parseModel <- function(model, x=NULL, Env=new.env(parent=.GlobalEnv),
                       EnvDummies = Env, add.na=NULL, idx.coord=NULL,
                       coord.names = character(0),
                       general.coord=character(0), factor=logical(0), data=NULL) {
                                        #
  ##  Print("parse", model, ls(envir=Env))
  
  ## check whether $model is already in list syntax
  if (is.list(model)) return(list(model=model, idx.coord=idx.coord))

  ## check whether $model has RMmodel syntax
  ##  Print(isRMmodel(model))

#  Print( idx.coord = idx.coord, coord.names=coord.names,
#        general.coord=general.coord, factor=factor)


  if (isRMmodel(model)) 
    return(buildCovList(model, x=x, Env=Env, EnvDummies=EnvDummies,
                        add.na=add.na,
                        idx.coord = idx.coord, coord.names=coord.names,
                        general.coord=general.coord, factor=factor, data=data
                        ))

  
   
 ## check whether $model has correct formula syntax
  if (!is(model, "formula")) stop(syntaxError) # 1 x isRMmodel aufgerufen
  
  Proj <- function(n, p=n) {
    val <- do.call(iR_P, list(proj=p, name=coord.names[n]))
    assign(coord.names[n], val, envir=Env)    
  }
  ProjG <- function(n, p=n) {
    val <- do.call(iR_P, list(proj=p, name=general.coord[n]))
    assign(general.coord[n], val, envir=Env)
  }
  Eval <- function(parse.txt) eval(parse.txt, envir=Env)

  Remove <- function(x) rm(list=x, envir=Env)

  assign("I", envir=Env, new(CLASS_RM, # 1 x isRMmodel aufgerufen
                .Data = function(x) {
                  eval(parse(text=deparse(substitute(x))),
                       envir = Env)
                },
                type = c('of manifold type'),
                isotropy = c('submodel dependent'),
                domain = c('submodel dependent'),
                operator = TRUE,
                monotone = 'submodel dependent monotonicity',
                finiterange = NA,
                simpleArguments = TRUE,
                maxdim = SUBMODEL_DEP,
                vdim = SUBMODEL_DEP
                ))

 # str(get("I", envir=Env))
#  print(get("I", envir=Env))
#
  
  rightSide <- as.character(model)
  rightSide <- rightSide[length(rightSide)]
  parseRS <- parse(text=rightSide)

  L <- length(coord.names)    
  Lg <- length(general.coord)
  notwithout <- logical(L)


##  Print("entering", coord.names, ":", general.coord, add.na)
  if (L > 0 ) {
    for (i in 1:L) { ## otherwise not clear which are there or not
      ##                and warnings will show up
      Proj(i)
      if (i <= Lg) ProjG(i)
    }
    rek <- function(idx) {
      for (i in idx) {
        Remove(coord.names[i])
        if (i <= Lg) Remove(general.coord[i])
      }     
      if (NotWithOut <- is(Try(Eval(parseRS)), "try-error")) {## 2 x isRMmodel aufgerufen
        for (i in idx) {
          Proj(i)
          if (i <= Lg) ProjG(i)
        }
        if (length(idx) == 1) notwithout[i] <<- TRUE
        else {
          half <- 1:as.integer(length(idx) / 2)
          rek(idx[-half])
          rek(idx[half])
        }
      }
    }
    rek(1:L)
  }
  Eval(parseRS)

##  Print(idx.coord, notwithout, L, coord.names, general.coord)
  
  idx.coord <- idx.coord | notwithout
  which.idx <- which(notwithout)
  which.factor <- factor[which.idx]  
#  Print(which.factor, CM, CM[1], CM[2])
#  print(model)
#  str(model)
  for (i in OneTo(length(which.idx))) {
     j <- which.idx[i]
    N <- coord.names[j]
#    Print(data, N)
    if (which.factor[i]) assign(N, data[[N]], envir=Env) else Proj(j, i)
    assign(N, N, envir=EnvDummies) # names of data columns
    ## may not appear in covariance models, except within formula
    ## so replacing numerical values by stop() should lead to errors
    ## if used within formula as they appear in the trend part
    if (j <= Lg) {
      G <- general.coord[j]
      if (which.factor[i]) assign(G, get(N, envir=Env), envir=Env)
      else ProjG(if (nchar(G) == 1) i else j, i)
      assign(G, G, envir=EnvDummies)
    }
  }
  
  ## extract tokens/summands
  tmpList <- list()
  chars <- strsplit(rightSide, "")[[1]]
  parToggle <- 0  
  token <- ""
  for (char in chars) {
    if (char == SYMBOL_PLUS && parToggle == 0) {
      tmpList <- c(tmpList, token)
      token <- ""			
    } else {
      if (char != " ") token <- paste(token, char, sep="")		
      if (char == SYMBOL_L_PAR) parToggle <- parToggle+ 1
      if (char == SYMBOL_R_PAR) parToggle <- parToggle- 1
    }
  }
  tmpList <- c(tmpList, token)

  summands <- vector("list", length(tmpList))
  for (i in 1:length(tmpList)) summands[[i]] <- removeParenthesis(tmpList[[i]])

  ## if (length(summands)==1) return(summands)  ## sonst steht unten paste(NULL), was "" gibt
  
  trendfct <- paste0(RM_TREND[1], "(")  
  isGenuineCovModel <- sapply(summands, FUN=function(s) {
    x <- Try(eval(parse(text=s),envir=EnvDummies)) 
    is(x, CLASS_CLIST) && isS4(x)## kostet Aufrufe isRMmodel
  })
 
  isFormalCovModel <- sapply(summands, FUN=function(s) {
    x <- Try(eval(parse(text=s),envir=Env))
    is(x, CLASS_CLIST) && isS4(x)## kostet Aufrufe isRMmodel
  })
  isRMTrend <- sapply(summands, FUN=function(s)
    substr(s, 1, length(trendfct)) == trendfct)

  ##  Print(add.na, chars, isGenuineCovModel, isRMTrend, summands, ls(envir=EnvDummies), get("a", envir=EnvDummies), ls(envir=Env), eval(parse(text=model),envir=EnvDummies))
  
  add.na <- add.na * (!isGenuineCovModel | add.na < 0) * !isRMTrend

##  Print(add.na, summands)

  listModel <- vector("list", length(summands) + 1)
  for (i in 1:length(summands)) {
    ## ksotest Aufrufe isRMmodel
    P <- buildFactorList(summands[[i]], Env=Env, EnvDummies = EnvDummies,
                         x=x, add.na=add.na[i],
                         isFormalCovModel = isFormalCovModel[i],
                         idx.coord = idx.coord, coord.names=coord.names,
                         general.coord=general.coord, factor=factor, data=data)
    idx.coord <- idx.coord | P$idx.coord
    listModel[[i+1]] <- P$model
  }
 ## Print(listModel)

  if (length(summands) == 1) listModel <- listModel[[2]]
  else listModel[[1]] <- SYMBOL_PLUS

  ##  Print(listModel, add.na)
  ## Print("ende parse")

  return(list(model=listModel, idx.coord=idx.coord))
}


buildCovList <- function(model, x, Env, EnvDummies, add.na,
                  idx.coord, coord.names, general.coord, factor, data) {
  ## Print("bCL", model)
  ## solves user defined model to a list by using
  ## in particular the models in RMmodels.R and the definition in
  ## RMmodelsBasics.R
  
  if (is.atomic(model) || is.list(model) ||
      is.language(model) || is.environment(model))
    return(list(model=model, idx.coord=idx.coord)) ## for recursive calling

  if (!isRMmodel(model)) stop('model must be of class ', CLASS_CLIST)
  
  name <- model@name
  if (name == RM_PLUS[1]) name <- SYMBOL_PLUS
  else if (name == RM_MULT[1]) name  <- SYMBOL_MULT

  P <- model@par.model
  if (length(P[[ADD_NA]]) > 0) {
    add.na <- P[[ADD_NA]]
    P[[ADD_NA]] <- NULL
##    Print(ADD_NA, add.na)
  }
  
  PD <- list(P=c(P, model@submodels),
             D=model@par.general[model@par.general != NO_DOLLAR_VALUE])
 
  for (j in 1:length(PD)) {
    P <- PD[[j]] 
    for (i in OneTo(length(P))) {
      X <- if (is(P[[i]], "formula")) 
             parseModel(P[[i]], x=x, Env=Env, EnvDummies=EnvDummies,
                        add.na = add.na,
                        idx.coord = idx.coord, coord.names=coord.names,
                        general.coord=general.coord, factor=factor, data=data)
           else buildCovList(P[[i]], x=x, Env=Env, EnvDummies=EnvDummies,
                             add.na=add.na,
                             idx.coord = idx.coord, coord.names=coord.names,
                             general.coord=general.coord, factor=factor, data=data)
      idx.coord <- idx.coord | X$idx.coord
      PD[[j]][[i]] <- X$model
    }
  }

  li <- c(list(name), PD$P)
  if (length(PD$D) > 0) li <- c(DOLLAR[1], PD$D, list(li))

  return(list(model=li, idx.coord=idx.coord))
}


removeParenthesis <- function(string) {
  n <- nchar(string)
  while (substring(string, 1, 1) == SYMBOL_L_PAR &&
         substring(string, n) == SYMBOL_R_PAR) {
           i <- 2; while(i <= n && substring(string, i ,i) == " ") i <- i + 1
           j <- n-1; while(j >= 1 && substring(string, j ,j) == " ") j <- j - 1
           string <- substring(string, i, j)
           n <- nchar(string)
         }
  return(string)
}


buildFactorList <- function(summand, Env, EnvDummies,
                            x, add.na, isFormalCovModel,
                            idx.coord, coord.names,
                            general.coord, factor, data)
                             {
  summand <- removeParenthesis(summand)
  C <- eval(parse(text=summand), envir=Env)
 
#  Print(summand)
 
  if (is.factor(C) || isFormalCovModel) {
##    Print(C, is.factor(C))
    if (is.factor(C)) {
      lev <- levels(C)
      if (length(lev) > MAXSUB^2 + 1)
        stop("max number of factors limited to ", MAXSUB^2)
      L <- list(RM_PLUS[1])
      i <- 2
      while (i <= length(lev)) {
        last <- min(length(lev), i + MAXSUB -1)
        plusList <- list(RM_PLUS[1])        
        while (i <= last) {
          ##          Print("addna")
          tmpList <- list(RM_COVARIATE)
          tmpList[[COVARIATE_NAME_NAME]] <- paste0(summand, i)
          tmpList[[COVARIATE_C_NAME]] <- as.numeric(C == lev[i]) ## 'data'
          tmpList[[COVARIATE_ADDNA_NAME]] <- add.na
          ##          tmpList[[COVARIATE_X_NAME]] <- x   ##
          if (add.na) tmpList[[COVARIATE_RAW_NAME]] <- TRUE 
          else tmpList[[COVARIATE_X_NAME]] <- x          
          plusList[[length(plusList) + 1]] <- tmpList
          i <- i + 1;
        }
        L[[length(L) + 1]] <- plusList
      }
      return(list(model = if (length(lev) == 2) tmpList
                          else if (length(lev) <= MAXSUB + 1) plusList
                          else  L,
                  idx.coord=idx.coord))
    } else {  ## (( && last))
      P  <- buildCovList(C, x=x, Env=Env, EnvDummies=EnvDummies,
                         add.na=add.na,
                         idx.coord = idx.coord, coord.names=coord.names,
                         general.coord=general.coord, factor=factor, data=data)
      tmpList <- P$model
      idx.coord <- idx.coord | P$idx.coord
    }

 ##   Print(add.na, P, tmpList, idx.coord)
    
    if (abs(add.na) > 1) {
      tL <- if ((tmpList[[1]] %in% RM_MULT)) tmpList[-1] else list(tmpList)    
      names <- sapply(tL, function(x) x[[1]])
      if (all(names != R_CONST)) {## user given model does not start with
        ## R_CONST (or ~1 had not been transformed to R_CONST(NA)), add "NA *"
        const <- list(R_CONST)
        const[[CONST_A_NAME]] <- NA
        const[[COVARIATE_NAME_NAME]] <- summand
        tmpList <- c(SYMBOL_MULT, tL, list(const))
      }
    }

  } else { ## pure numbers
    if (length(C) > 1) {      
      tmpList <- list(RM_COVARIATE)         
      tmpList[[COVARIATE_C_NAME]] <- C
      tmpList[[COVARIATE_ADDNA_NAME]] <- add.na
      if (add.na) tmpList[[COVARIATE_RAW_NAME]] <- TRUE 
      else tmpList[[COVARIATE_X_NAME]] <- x   
    } else {
#      Print(add.na, C); stopifnot(add.na != 0)
      tmpList <- list(R_CONST)
      tmpList[[CONST_A_NAME]] <-
        if (is.finite(C) && C == 1 && add.na) NA else as.numeric(C)
    }
    tmpList[[COVARIATE_NAME_NAME]] <- summand
  }

  return(list(model=tmpList, idx.coord=idx.coord))
}

