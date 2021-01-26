
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
## You should have received a cop of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  

#Zaehler <- 0


## NOTE:
## "definitions" in the man page is called here "symbol"

CopyDotsTo <- function(dots, Env, first=FALSE) {
  assign("..dots..", dots, envir=Env)
  if (length(dots)>0) {
    if (first && any(0 != sapply(dots, function(x) length(oldClass(x)))))
      stop("arguments in '...' may not have a class")
    txt <-  paste(names(dots), "<- ",
                  "..dots..[[", 1:length(dots), "]]", collapse=";")
    eval(parse(text = txt), envir=Env)
  }
}
  
CheckDots <- function(dots, Names, txt) {
  n <- names(dots)
  if (any(n %in% Names))
    stop(txt, " may not be named after a data column name. ",
         "Here, coincidence is for '",
         paste(names(dots)[names(dots) %in% Names], collapse="', '"), "'.")
  reserved <- c("I", "c", "cbind")
  if (any(n %in% reserved))
      stop("the reserved word", S(reserved), " '",
           paste(reserved, collapse="', '"),
           "' may not be used in the dots (...)")
}

  

ExtractNames <- function(Names, data, dont.add.data, model, RFopt, dots, Env) {
  coord.opt <- RFopt$coords
  link.coord <-link.data <- data.names <- character(0)
  data.trafo <- NULL
  is.var <- is.x <-  integer(0)
  factor <- logical(0)
  data.initial <- coord.opt$dataframe_initial
  repet <- vdim <- 0

  ## determine names of variables on the left and on the right of ~
  ## replace ~x etc by R.p etc
  if ((coord.initial <- coord.opt$coord_initial) != "") {
    default.link.coord <- paste0(coord.initial, 1:coord.opt$max_coord)
    nice.link.coord <- c(coord.opt$cartesian_names, coord.opt$earth_coord_names)
    link.coord <- c(nice.link.coord, default.link.coord)
    names(link.coord) <- c(rep("coord", length(coord.opt$cartesian_names)),
                           rep("earth", length(coord.opt$earth_coord_names)),
                           rep("default", coord.opt$max_coord))
  }
  
#  if (length(data) == 0) {
#    if (length(coord.opt$varnames) == 0 && data.initial != "") {
#      link.data <- paste0(data.initial, 1:coord.opt$max_columns)
                                        #    }                                        #  } else { ##  data given
  is.var.eFN <- NULL
  if (length(data) != 0) {
    ## first extract data column names
    if (is.list(data) && !is.data.frame(data)) data <- data[[1]]
    if (is(data, "RFsp")) {
      data.names <- colnames(data@data)
      if (data@.RFparams$n > 1)
        data.names <-
          sapply(data.names, ## AUTHOR: Sebastian Gross, 26.08.2011
                 FUN= function(x) strsplit(x, "\\.n[[:digit:]]+$")[[1]][1])
      repet <- 1:data@.RFparams$n
      vdim <- data@.RFparams$vdim
      is.var <- 1:ncol(data@data)
    } else {
      if (is.data.frame(data) || is.matrix(data)) data.names <- colnames(data)
      data.names[data.names == ""] <- VOID # ansonsten fehler in assign(i,i,...)
      if (dont.add.data)
        CheckDots(dots=dots, Names=data.names, txt="additional arguments")
      nc <- length(data.names)
      if (nc == 0) {
        nc <- if (is(data, "RFsp")) ncol(data@data) else ncol(data)
        if (length(nc) == 0) stop("size of data cannot be determined")
      }
      if (data.initial != "") link.data <- paste0(data.initial, 1:nc)

      if (length(data.names) > 0) {   
        is.var  <- extractFromNames(RFopt=RFopt, cn=data.names)
        if (is.matrix(is.var)) {
          is.var.eFN <- is.var
          is.var <- c(is.var)
          vdim <- nrow(is.var.eFN)
          repet <- ncol(is.var.eFN)
          data.names[is.var] <- rownames(is.var.eFN) ## unify names
        }
        is.x <- extractFromNames("coord", RFopt=RFopt, cn=data.names)
        if (any(data.names %in% default.link.coord)) {
          L <- min(length(data.names), length(default.link.coord))
          if (!all(data.names[1:L] == default.link.data[1:L])) {
            link.coord <- nice.link.coord
            string <- char <- if (substr(coord.initial, 1, 1)=="X")  "C" else "X"
            setRFoptions(coords.dataframe_initial = string)
            while (any(string == substr(Names, 1, length(string))))
              string <- paste0(string, char)
            warning("column names of the data should not be the same as the default names for coordinates. Better change the default names by RFoptions((coord.initial='", string, "')")
        }
        }
        if (any(data.names %in% link.data)) {
          L <- min(length(data.names), length(link.data))
          if (!all(data.names[1:L] == link.data[1:L])) {
          data.coord <- NULL
          string <- char <- if (substr(coord.initial, 1, 1)=="V") "D" else "V"
          setRFoptions(coords.coord_initial = string)
          while (any(string == substr(Names, 1, length(string)))) string <- paste0(string, char)
          warning("column names of the data should not be the same as the default names for data names. Better change the default names by RFoptions((data.initial='", string, "')")
          }
        }      
      }
      
    }

    ## !! folllowing assignment also needed in PrepareModel !!
    DN <- data.names
    ##  Print(data, DN)
    for (i in DN) assign(i, i, envir=Env) ## will allow formalae on left hand side
  }
  
  if (any(Names %in% c(data.names, link.data, link.coord)))
    stop("'coordnames' and 'varnames' of RFoptions may not be used as variable names; neither colnames(data) nor the names created with dataframe_initial and coord_initial. I.e. forbidden are\n", paste(c(data.names, link.data, link.coord), collapse=", "))
  
  if (is(model, "formula")) { ## ueberschreibt alle vorigen defaults!!
    ## try identify data.names, coord.names from model point of view (from left side!)
    CM <- as.character(model)
    
    if (length(CM) == 3) { ## data.varnames given by left hand side of formula
      leftSide <- parse(text=CM[2])
      
      used.varnames <- eval(leftSide, envir=Env) # indeed used ones
        is.var <- match(used.varnames, DN)## geordnet !!!
      if (length(is.var) == 0 && length(link.data) > 0) {
        rm(list=data.names, envir=Env)
        DN <- link.data
        ## entweder ueber colnames oder D*, nicht gemischt
        for (i in DN) assign(i, i, envir=Env)
        used.varnames <- eval(leftSide, envir=Env)
        is.var <- match(used.varnames, DN)
        if (length(is.var) > 0) data.names <- DN
      }
      
      if (any(is.na(is.var))) { ## fiktive Namen vom User
        assign("c", function(...) as.character(as.list(match.call())[-1]),
               envir=Env) ## Umdefinition von c(...) --- warum?
        used.varnames <- rawTry(leftSide, envir=Env)
        used.varnames <- if (rawError(used.varnames)) NULL else CM[2]
        is.var <- match(used.varnames, data.names)
        if (!all(is.na(is.var))) {
          is.var <- is.var[!is.na(is.var)] ## NAs entweder aufgrund von Fehler
          ## von User oder wegen Formeln von User, kann erst bei UnifyData abgefangen werden
          data.trafo <- leftSide
        }
        rm("c", envir=Env)
        setRFoptions(coords.varnames = used.varnames) ## 16.12.20 noetig?
      }
      vdim <- length(is.var)
      if (!is.null(is.var.eFN)) { ## then it is a
        ## matrix an names have been doubled. now the model proably refers to
        ## a subset of the namesn
        if (!all(is.var  %in% c(is.var.eFN)))
          stop("model variables do not match the given data")
        ## achtung! is.var haben durch Nutzer vorgegebener Ordnung
        is.var <- as.vector(is.var.eFN[pmatch(is.var, is.var.eFN[,1]), ])
      }
    }
    if (length(is.var) > 0) names(is.var) <- DN[is.var]
  }

  if (length(data.names) > 0 && !is(data, "RFsp")) {
    factor <- IsFactor(data, data.names, is.var=is.var) ## cannot be x-coord
    is.x <- setdiff(is.x, which(factor))
  }
  
  return(list(factor=factor,
              data.names=data.names,
              is.var=is.var,
              is.x=is.x,
              link.coord=link.coord,
              link.data=link.data,
              data.trafo = data.trafo,
              vdim = vdim,
              repet = repet
              ))
}



AddTransform <- function(M, model, params, orig.params, Names, nNA,
                         revOrdering,
                         symbols, formulae, genuine.formulae, 
                         factor,
                         xdim,  MODEL_AUX, add.na, Env, EnvDummies) {
  ##  Print("Entering AddTrafo")
  
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
  
  M300 <- parseModel(model=model, Env=Env, EnvDummies = EnvDummies,
                     add.na=add.na)$model
  
  ## prepare model to determine the indexing between R and C
  ## of the parameters that are NA
  ## rf_interface.cc, includeparam, these values cannot appear
  ## as parameter values
 
  model.vdim <- 0L
  NAs.in.model <- ## user might have given NA directly in the model
    .Call(C_GetNAPositions, MODEL_AUX, list("Dummy", M300), M$C_coords,
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
  CtoR <- .Call(C_GetNAPositions, MODEL_AUX, NULL, M$C_coords,
                values, double(0), FALSE, model.vdim, # side effect !!
                0L)
  n.CtoR <- length(CtoR)
  if (n.CtoR > 0) {
    Check.CtoR(CtoR, MODEL_AUX, values)
    RtoC <- rep(NA, (totalNA + tot.NAs.in.model))
    RtoC[CtoR] <- 1:n.CtoR

    genuine.symbols <- symbols & (nNA > 0)
    constants <- which(symbols & !genuine.symbols)
    S2P <- which(genuine.symbols)
    S2Pxx <- which(c(genuine.symbols, rep(TRUE, length(NAs.in.model))))
    S2C0 <- cbind(cumNA[-length(cumNA)][S2Pxx] + 1, cumNA[-1][S2Pxx])
    S2C <- unlist(apply(S2C0, 1, function(x) list(RtoC[x[1]:x[2]])),
                  recursive=FALSE)
    im.S2C <- unlist(S2C) # im = (flattened) image

    ##    Print(M300, NAs.in.model, tot.NAs.in.model, cumNA, values, CtoR, RtoC, genuine.symbols, S2P, S2Pxx, S2C0, S2C, im.S2C)
    
    if (any(is.na(im.S2C))) {
      stop("incorrect 'params'. Probably not all newly introduced variables have been defined in 'RMdeclare'.")
    }
    CmidS <- integer(n.CtoR)

    
    CmidS[im.S2C] <- rank(im.S2C)
    ## S2CmidS is the sequence, by which input variables come in:
    ## (input variables have the same ordering than the one
    ## inside the ordering of all NAs and NaNs)
    S2CmidS <- lapply(S2C[OneTo(length(S2P))], function(x) CmidS[x])
        
    G2P <- which(genuine.formulae)
    n.genuine.formulae <- length(G2P)
    G2C <- cbind(cumNA[-length(cumNA)][G2P] + 1, cumNA[-1][G2P])
    G2C <- unlist(apply(G2C, 1, function(x) list(RtoC[x[1]:x[2]])),
                  recursive=FALSE) ## to be save that it always a list
    used <- sapply(G2C, function(x) { finite <- !any(is.na(x));
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

#    Print("NACH MOCK", ls(envir = Env))
    
    G2text <- character(n.genuine.formulae) # for each formula a text
    for (i in OneTo(n.genuine.formulae))
      G2text[i] <- as.character(as.expression(orig.params[[ G2P[i] ]][[2]] ))
    ordering <- order(revOrdering[G2P])
    
    def.txt <- paste0(
        sep="", "function(variab) {\n\t",
        if (length(constants) > 0)
          paste(Names[constants], "<-", params[constants], "\n\t", collapse="")
       ,
        if (length(S2P) > 0) ## otherwise no variable is given
          paste(namesS2P, p, "<- variab[", S2CmidS, "]", A, "\n\t", collapse="")
       ,
        if (length(mock) > 0)
          paste(Names[mock], "<-", params[mock],"\n\t",collapse="")
       ,
        if (length(G2P) > 0)
          paste(Names[G2P][ordering], "<-", G2text[ordering],"\n\t",collapse="")
    )
    estim.txt <- paste0(
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
    all <- c(namesS2P, Names[mock], Names[G2P])
    return.txt <- paste0("list(", paste(all, sep=" = ", all, collapse=","),
                         ")\n}\n")
    
    fctn <- eval(parse(text=paste0(def.txt, estim.txt)), envir=NULL)
    fctn.env <- new.env(parent=baseenv())
    environment(fctn) <- fctn.env
    CopyDotsTo(get("..dots..", env=Env), fctn.env)

    params.fctn <- eval(parse(text=paste0(def.txt, return.txt)), envir=NULL)
    environment(params.fctn) <- fctn.env
  
    
    ## isNA: logical, Laenge (CtoR); true falls symbol
    isNA <- rep(FALSE, n.CtoR)
    isNA[im.S2C] <- TRUE
    
    CmidN2RmidN <- .Call(C_GetNAPositions, MODEL_AUX, NULL, M$C_coords,
                         values[-which(rep(formulae, nNA))], double(0),
                         FALSE, model.vdim, 0L)
    if (length(CmidN2RmidN) != sum(nNA[S2P]) + tot.NAs.in.model)
      stop("number of parameters declared to be estimated differs from the number of non-linear parameters to be estimated. In other words, for any declaration of the form 'xy = NA' in 'params' needs exactly one counter part within in the (RM)model defintion of the form 'arg = xy'")
  } # at least 1 variable to be filled on C level
  
  M$transform <- list(isNA=isNA, fctn=fctn, params.fctn=params.fctn)
  M$vdim <- model.vdim
  return(M)
}


IsFactor <- function(data, data.names, is.var) {
  if (is(data, "factor")) TRUE
  else if (is.data.frame(data)) {
    f <- sapply(data, is.factor)
    if (any(f[is.var])) {
      f[-is.var] <- FALSE
      stop("dependent variable (", paste(data.names[f], collapse=","),
           ") seems to be a factor")
    }
    f
  } else if (is.matrix(data) || is(data,"RFsp")) rep(FALSE, ncol(data))
  else if (is.vector(data)) FALSE
  else {
    stop("unknown data format")
                ## simulation from randomFields -- unclear, what / how should
    ## any relevant information be extracted  TO DO???
  }
}



RemoveBlanks <- function(u) sub("^ *", "", sub(" *$", "", u))
Check.declare <- function(rightSide, Names, nNA, symbols) {
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
  if (length(f) > 0) names(f) <- nn
  nn2 <- nn[f >= 2]             
  ##
  
  if (any(f != 1)) {
    stop("We call an expression in 'params' that includes NAs on the right hand side of a definition. (In particular, a definition is not given by a formula.) The name of each definition must appear exactly once in the RMmodel definition as a \"value\" of an argument. This is not the case for ",
         paste0("'", nn[f != 1], "'", collapse=","),
         ". Likely, you should do one of the following:\n", 
         if (any(f == 0)) {
           if (length(d) == 2)
             paste("* include", paste0("'", nn[f == 0], "'", collapse =","),
                   "in 'RMdeclare'.\n")
           else paste0("* add 'RMdeclare(", paste0(nn[f==0], collapse =","),
                       ")' to the model.\n")
         },
         if (any(nn[f == 0] %in% unlist(z[sapply(z, length) == 1])))
           paste("* name the argument", S(which(f==0)),"'",
                 paste(nn[f == 0], collapse="', '"), "'."),
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
      N <- names(p)
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
    link <- GetParam(internRFgetModelInfo(MODEL_AUX), NULL)
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
                          ## transform; if not given, RMcovariate will fail,
                          ## in general (for call from
                          ## internRFgetModelInfo_model, OK)
                          trend=NULL, ## a model indicating a trend -- a further
                          ## option to define a trend; other possiblities are
                          ## trend=TRUE in '+' or to hope for automatic
                          ## detection
                          data = NULL, # data, needed here to detect on the fly
                          ## which columns are used are variables or coordinates
                          add.na=FALSE, # for trends: should NAs be added?
                          ## THis is needed for estimation in linear models
                          ## add.na == TRUE: only RMcovaiate gets additional NA 
                          ## add.na == 2: all linear functions get additional NA
                          return_transform = FALSE,
                          xdim = 0,
                          Env = new.env(parent=.GlobalEnv), ## nur im
                          ##                       rekursiven Aufruf verwendet
                          EnvDummies = new.env(parent=Env),
                          dont.add.data = TRUE ## false needed only in of covariates in kriging
                          ) {

  ##  Print("enering prep", model, missing( params), data, missing(x), if (!missing(x)) x)
  
 # stopifnot(Zaehler <- Zaehler + 1 < 2)
  ## (1) rffit->unifydata
  ## (2)
  ## (3)
  ## (4)
  
  RFopt <- getRFoptions(getoptions_=c("basic", "general", "messages", "coords"))
  coord.opt <- RFopt$coords
  default.max.number.data <- RFopt$coords$max_columns
  default.max.number.coord <- RFopt$coords$max_coord
  dummy.matrix.fill <- Inf
  if (length(xdim) == 0) xdim <- 0

  ## before parseModel is evaluated since columns might be added
  ## in case of kriging and "..." give

  dots <- list(...)
  CopyDotsTo(dots, Env, first = TRUE)
  if (!missing(model)) {
    if (is(model, CLASS_FITLIST)) {
      Help("mle")
      model <- if (isS4(model)) model["ml"] else model[["ml"]]
    }
    if (is(model, CLASS_SINGLEFIT))
      model <- if (isS4(model)) model@model else model$model
    if (is(model, CLASS_RM))
      model <- do.call(model, list())
    if (is(model, CLASS_CLIST) &&
        ((isS4(model)) && length(model@name) == 0) || (is.list(model) && length(model) == 0))
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
                     RMshape(Trend[length(Trend)]), ")")
    }

    if (missing(model)) {
      model <- if (isFormula) eval(parse(text=paste(Tv, "~", trend)))
               else RMtrendplus(RMshape(trend), add.na=add.na)
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
          assign("..trend..", RMshape(trend), envir=Env)
          assign("..trend..", RMshape(trend), envir=EnvDummies)
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
      
      else model <- RMplus(RMtrendplus(RMshape(trend), add.na=add.na),
                           RMmodelplus(trend = FALSE, model))
    }
    
    return(PrepareModel2(model, params=params, 
                         further.model = further.model, arg=arg, x=x,
                         data=data,
                          add.na = FALSE,
                         return_transform = return_transform,
                         xdim = xdim,
                         Env = Env,
                         EnvDummies = EnvDummies
                         ))
  }
  
  if (missing(model) || is.null(model)) stop("'model' must be given.")
  if (missing(params)) params <-NULL
    
  link.names <- NULL
  genuine.formulae <- formulae <- FALSE  
  values <- double(0)  
  ## warn <- options(warn=2)
  ##  options(warn
  Names <- if ( length(params) > 0) names(params) else NULL
  DataNames <- ExtractNames(Names=Names, data=data, dont.add.data=dont.add.data,
                            model=model, RFopt=RFopt, dots=dots, Env=Env)
  data.names <- DataNames$data.names

  if (!missing(data) && length(data) > 0) {
    data1 <- data[[1]]
    if (length(colnames(data1)) == 0 && length(data.names) > 0 &&
        is.matrix(data1)) colnames(data1) <- data.names
  } else data1 <- NULL

 
  is.x <- DataNames$is.x
  link.coord <- DataNames$link.coord
  nlc <- names(link.coord)
  wechsel <- which(nlc[-1] != nlc[-length(nlc)])

##  Print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!", params)

  if (length(params) > 0) {
   
    ## simple general stuff for params; reset params in case of further.model(s)
    ## latter are the lower bound model or the uppper bound model
    orig.params <- params
    n.params <- length(params)
    Names <- names(params)
    if (!all(nchar(Names) > 0)) stop("parameters are not all given by name")
    CheckDots(dots, Names=Names, txt="The left sides within 'params'")
    if (anyDuplicated(Names))
      stop("names on the left sides in 'params' must be unique.")
    if (any("I" == Names)) stop("Do not use the reserved word 'I' in 'params'")
    
    formulae <- logical(n.params) ## FALSE
    for (i in 1:n.params) formulae[i] <- is(params[[i]], "formula")
    symbols <- !formulae

    if (!missing(further.model) && !is.null(further.model)) {
      if (is.list(further.model) && is.character(further.model[[1]]))
        stop("Do not mix different styles of model definition.\n  See ?RFformula and ?RFformulaAdvanced for details.")
      FurtherModels <- "In bounds, initial values and parscales, "
      for (n in names(further.model)) {
        idx <- which(n == Names)
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

    ##  Print(nNA, formulae, symbols, ls(envir=Env)) 
  
    if ((n.formulae <- sum(formulae)) > 0) {
      ## identify genuine.formulae and #variables to be estimated; prepare Env
      ## 18.12.20: third cathegory: 
      ##
      ## params ist unveraendert orig.params, ausser !missing(further.model)
      ## letzteres schliesst wechselseitig return_transform aus!
      revOrdering <- rep(0, n.formulae) ## reverse ordering of formulae
      genuine.formulae <- formulae
      zaehler <- 0
      for (last in c(FALSE, TRUE)) {
        repeat {
          repetition <- FALSE 
          for (i in 1:n.params) {
            if (is(params[[i]], "formula")) {
              value  <- rawTry(eval(as.expression(params[[i]][[2]]),envir=Env))
              if (rawError(value)) { # error
                if (last) {
                  stop("This situation should not appear. Please check whether your definitions are sound and/or contact ", AUTHOR)
                  value <- NaN 
                  params[[i]] <- value
                  assign(Names[i], value, envir=Env)
                  nNA[i] <- 1
                  revOrdering[i] <- zaehler <- zaehler + 1 # andere Konstruktion
                  ##geht nicht, da unten formeln zu nicht-Forelm erklaert werden
                }
              } else {
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
                repetition <- !last
              } ## not try-error
            } ## is.formula
          } ## for i
          if (!repetition) break
        } ## repeat
      } ## for last

    } # if any(formulae)
    ##

    if (is(model, "formula") && length(params) > 0 && sum(nNA[symbols]) > 0) {
      CM <- as.character(model)
      Check.declare(CM[length(CM)], Names, nNA, symbols)
    }

    ## search for constructions of the form 'X2=D4' in params
    ## all are previously found indication of coords are ignored!
    i <- 1
    len <- length(link.coord)
    while(i <= len && !exists(link.coord[i], envir = Env, inherits=FALSE))
      i <- i + 1
    if (i <= len) {
      is.x <- NULL ## delete previous knowledge about coordinates
      add <- wechsel - c(0, wechsel[-length(wechsel)])
      restlaenge  <- len - wechsel[length(wechsel)]
      step <- c(rep(add, add), restlaenge:1)
      ##      sockel <- rep(c(0, wechsel), c(add, restlaenge)) ## nicht loeschen
      k <- i
      while(k <= len) { ## all numbers in a sequence must be givem
        if (!exists(link.coord[k], inherits=FALSE)) { ## should be end of list
          end <- k + step[k] ## different coord system name start from one again
          k <- k + 1 ## but is there a gap ?
          while(k<end && !exists(link.coord[k], envir = Env, inherits=FALSE))
            k <- k+1
          if (k<end) stop("The definitions of all coordinates must be given.")
          next
        }
        value <- get(link.coord[k], envir=Env)
##        Print(ls(envir=Env), value, link.coord, k, link.coord[k])
        if (!is.character(value) || length(value) != 1)
          stop("value of '", link.coord[k],
               "' is neither a column name nor an abstract reference",
               "to a column: ",
               paste0("'",
                      internalRFoptions(getoptions_="coords")$cord_initial,
                      1:2, "'",
                      collapse=", "), "...")
        idx <- which(data.names == value)
        if (DataNames$factor[idx])
          stop("coordinate '", link.coord[k], "' may not be a factor")
        if (length(idx) != 1 || is.na(idx))
          stop("'", link.coord[k], "' does not have an apropriate value.")
        is.x <- c(is.x, idx) ## Konstruktion erspart den sockel oben
        link.names <- c(link.names, link.coord[k])
        k <- k + 1
      }
    }
  }

  ## prepare environment for the use of x,y,z,T (or other column names
  ## for coordinates) as values to be replaced in the linear model part
  if (length(is.x) > 0) {    
    for (k in 1:length(is.x)) {
      n <- data.names[is.x[k]]
      value <- do.call(iR_P, list(proj=k, name=n))
      assign(n, value, envir=Env)
      if (length(link.names) > 0) assign(link.names[k], value, envir=Env)
    }
  } else if (length(link.coord) > 0){ ## is.x == NULL, just take the link.coord
    wechsel <- c(0, wechsel) ## change in coordinates system names
    for (w in 2:length(wechsel)) {
      sockel <- wechsel[w-1]
      for (k in (sockel+1):wechsel[w])
        assign(link.coord[k], envir=Env, ## noch nicht
               ## sauber, da eigentlich noch new=EARTH_COORD gesetzt werden muesste
               ## fuer link.coord = lon/lat
               do.call(iR_P, list(proj=k-sockel, name=link.coord[k])))
    }
  }

  factor <- DataNames$factor
  is.covariate <- EnvTest <- NULL
  assign("c", cbind, envir=Env) ## cbind multivariate covariates
  is.undetectable <- which(data.names == VOID)
  is.unclear <- (1:length(data.names))[-c(is.x, DataNames$is.var,
                                          is.undetectable)]
  dn <- data.names[is.unclear] ## in multivariate case,
  ##     names might be used for more than 1 column; only feasible when
  ##     if data are passed by matrix. data.frames don't allow it
  if (length(is.unclear) > 0) names(is.unclear) <- dn
  
  if (length(dn) > 0) {
    unidn <- unique(dn)
    simple <-  length(dn) == length(unidn)
    CopyDotsTo(dots, Env, first = TRUE)
    
    covariates <- data1[, is.unclear, drop=FALSE]
    for (i in unidn) {## covariate
      d <- data1[ , if (simple) i else which(i == dn), drop=FALSE]
      d <- if (is.factor(d[[1]])) d[[1]] else as.matrix(d)

      if (dont.add.data) assign(i, envir=Env, d)

      else if (exists(i, envir=Env, inherits=FALSE)) { ## koennte sein,
        ##                                  Nutzer hat mehr
        is.covariate <- c(is.covariate, i) ## string zeilen eingegeben
        dot.covariate <- get(i, envir=Env)      
        ## in kriging the given and the new coordinates may contain
        ## covariates, even in the covariance model (non-stat variance!)
        ## so the data must be jointly passed to RMcovariate
        ## RMcovariate must be able to distinguish.

        if (is.factor(d) && !is.character(dot.covariate)) {
          if (!is.factor(dot.covariate))
            stop("a component cannot be a factor in one part and not a factor in the other")
          assign(paste0(i, "..factor.."), d, envir=Env)
        } else if (is.numeric(dot.covariate)) {
          dot.covariate <- as.matrix(dot.covariate)
          
          L <- nrow(dot.covariate) - nrow(d)
          
          if ((L <- length(dot.covariate) - nrow(data1)) > 0)
            d <-  do.call(rbind, c(list(d), as.list(rep(dummy.matrix.fill, L))))
          else if (L < 0)
            dot.covariate <- do.call(rbind, c(list(dot.covariate),
                                          as.list(rep(dummy.matrix.fill, -L))))
          colnames(d) <- rep("..", ncol(d))
          assign(i, envir=Env, cbind(dot.covariate, d))
        } else if (!is.character(dot.covariate))
          stop("unexpected value in creating the model", CONTACT)
      }
    }
    for (i in data.names) assign(i, i, envir=EnvDummies)
    for (i in link.coord) assign(i, i, envir=EnvDummies)

    ## erst jetzt !! factor muss vorher noch gespeichert werden
    is.unclear <- setdiff(is.unclear, which(factor))#da factor sicher covariate
    ##                                     und komplett extra behandelt wird  
    if (!dont.add.data) is.unclear <- setdiff(is.unclear, is.covariate)## mues-
    ##                                                     sen alle erkannt sein
    if (length(is.unclear) > 0) { ## jetzt ohne factors und
      ##                                       nur bei dont.add.data
      EnvTest <- new.env(parent=.GlobalEnv)
      for (i in is.unclear) {
        n <-data.names[i]
        m <- get(n, envir=Env)
        assign(n, envir=EnvTest, m[min(2, nrow(m)), , drop=FALSE]) # to get
        ## testing fast
      }
    }
  }

#  if (length(ls(envir=Env)) >= 15)
#    Print(ls(envir=Env), get("Z1", envir=Env), dots, data1, data, dn, is.unclear, data.names, is.x, DataNames$is.var)
  ##stopifnot(length(ls(envir=Env)) < 15)
 
  ## Prepare the identity operator (might have been used by the user)
  assign("I", envir=Env, new(CLASS_RM, # 1 x isRM model aufgerufen
                .Data = function(x)
                  eval(parse(text=deparse(substitute(x))), envir = Env),
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

  
  unclear.names <- data.names[is.unclear]

  
  M <- parseModel(model=model, Env=Env, EnvDummies = EnvDummies, add.na=add.na,
                  unclear=unclear.names, EnvTest=EnvTest)
  class(M$model) <- CLASS_CLIST

  M$is.x <- is.x
  M$is.factor <- which(DataNames$factor)
  M$is.var <- DataNames$is.var
  M$data.trafo <- DataNames$data.trafo
  M$data.names <- data.names
  idx <- unclear.names %in% M$unclear
  M$is.unclear <- c(is.unclear[idx], is.undetectable)
  M$is.covariate <- union(which(data.names %in% is.covariate),
                          is.unclear[!idx]) ## nur die, die
  ## ueber data erfasst wurden, nicht diejenigen, die direkt eingegeben wurden
  M$vdim <- DataNames$vdim
  M$repet <- DataNames$repet

  
  if (xdim != 0 && length(M$is.x) != 0 && length(M$is.x) != xdim) {
    stop("dimension mismatch", CONTACT)
  }
  
  M$C_coords <- if (length(x) > 0) trafo.to.C_UnifyXT(x)
                else if (length(M$is.x) > 0) {
                  C_UnifyXT(lapply(data, function(d) d[, M$is.x, drop = FALSE]),
                            grid = FALSE)
                } else if (xdim > 0) {                    
                  xx <- matrix(1, nrow=3, ncol=xdim)
                  if (length(data1) > 0)
                    xx[3, 1] <- nrow(data1) ## nur erste Dim echt
                  C_UnifyXT(xx, grid=TRUE)
                } else stop("'x' not given (nor 'dim')")

                                         #
  if (return_transform && sum(genuine.formulae) > 0) {
    stopifnot(missing(further.model) || is.null(further.model))
    M <- AddTransform(M=M, model=model, 
                      params=params, orig.params=orig.params,
                      Names=Names, 
                      nNA=nNA, revOrdering=revOrdering,
                      symbols=symbols,formulae=formulae,
                      genuine.formulae=genuine.formulae,
                      factor= factor,
                      xdim=xdim, 
                      MODEL_AUX=MODEL_AUX,
                      add.na=add.na,
                      Env=Env, EnvDummies = EnvDummies)
  } else if (DataNames$vdim != 0) M$vdim <- DataNames$vdim
  else {
    model.vdim <- integer(1)
    if (length(M$C_coords) > 0) {
      .Call(C_GetNAPositions,
            MODEL_AUX, list("Dummy", M$model),
            M$C_coords, as.double(NA), as.double(NA),
            FALSE,
            model.vdim, # side effect !!
          0L)
    }
    M$vdim <- model.vdim
    if (DataNames$vdim > 0 && model.vdim != DataNames$vdim)
      stop("the multidimensionality of the data and the model do not seem to be consistent")
  }

  if (any(M$add.na)) Help("addNA")

##  Print(M, data,DataNames, return_transform && sum(genuine.formulae) )
  
  M
}

detect.covariates <- function(S, unclear, EnvTest) {
  idx <- 1:((length(unclear) + 1) / 2)
  new.unclear <- NULL
  for (j in 0:1) { #
    m <- unclear[idx]
    if (length(m) == 0) next
    for (i in m) assign(i, list(get(i, envir=EnvTest))) ## disable values
    try <- Try(eval(parse(text=s),envir=Env))
    for (i in m) assign(i, get(i, envir=EnvTest)[[1]]) ## restore
    new.unclear <- c(if (is.numeric(try)) unclear[idx]
                     else if (length(unclear) > 1) 
                       detect.covariates(S=S, unclear=unclear[idx],
                                         EnvTest=EnvTest),
                     new.unclear)
    idx <- -idx
  }
  return(new.unclear)
}

covariate.names <- function(m) {
  if (m[[1]]==RM_COVARIATE || m[[1]]==R_CONST) return(m[[COVARIATE_NAME_NAME]])
  unique(unlist(lapply(m, function(x) if (is.list(x)) covariate.names(x))))
}

parseModel <- function(model, Env, EnvDummies, add.na=NULL, 
                       unclear=NULL, EnvTest=NULL) {
                                        #
  ##  Print("parse", model, ls(envir=Env), is.list(model),  is(model, CLASS_CLIST), isS4(model))
 
  if (!isS4(model) && is(model, CLASS_CLIST)) {
    ans <- covariate.names(model)
    ##    Print(model, ans, unclear, setdiff(unclear, ans))
    return(list(model=model, unclear=setdiff(unclear, ans)))
  }

 
  
  if (isS4(model) && is(model, CLASS_CLIST))
    return(buildCovList(model, Env=Env, EnvDummies=EnvDummies,
                        add.na=add.na, unclear=unclear, EnvTest=EnvTest))
  
  if (!is(model, "formula")) stop("Malformed model expression -- maybe you have used a wrong or obsolete definition, or just used an incorrect option name or incorrect ordering of the arguments. See ?RMmodel for the model definition. Check manual for further information (RMmodel, RFsimulate)")

  rightSide <- as.character(model)
  rightSide <- rightSide[length(rightSide)]
  parseRS <- parse(text=rightSide)

  
  ## extract tokens/summands
  summands <- list()
  chars <- strsplit(rightSide, "")[[1]]
  parToggle <- 0  
  token <- ""
  for (char in chars) {
    if (char == SYMBOL_PLUS && parToggle == 0) {
      summands <- c(summands, token)
      token <- ""			
    } else {
      if (char != " ") token <- paste(token, char, sep="")		
      if (char == SYMBOL_L_PAR) parToggle <- parToggle+ 1
      if (char == SYMBOL_R_PAR) parToggle <- parToggle- 1
    }
  }
  summands <- c(summands, token)
  summands <- lapply(summands, removeParenthesis)
  
  isGenuineCovModel <- sapply(summands, FUN=function(s) { # pos def
    x <- Try(eval(parse(text=s),envir=EnvDummies))    
    isS4(x) && is(x, CLASS_CLIST) ## kostet Aufrufe 
  })
 
  isFormalCovModel <- sapply(summands, FUN=function(s) { # parts using data or
    ##                                                     reference to coord
    x <- Try(eval(parse(text=s),envir=Env))
    isS4(x) && is(x, CLASS_CLIST) ## kostet Aufrufe 
  })
  
  trendfct <- paste0(RM_TREND[1], "(")
  isRMTrend <- sapply(summands, FUN=function(s)
    substr(s, 1, length(trendfct)) == trendfct)

  add.na <- add.na * (!isGenuineCovModel | add.na < 0) * !isRMTrend
  
  listModel <- vector("list", length(summands) + 1)
  m <- 1
  for (k in 1:length(summands)) {   ## ksotest Aufrufe 
    m <- m + 1;
    S <- summands[[k]]
    ##    Print(k, S, ls(envir=Env), Env, .GlobalEnv)
    C <- eval(parse(text=S), envir=Env) ## could be a function of the data col

    ## Print(S, C, is.factor(C))
    
    if (is.factor(C)) {
      lev <- levels(C)     
      if (length(lev) > MAXSUB^2 + 1)
        stop("max number of factors limited to ", MAXSUB^2)
      if (length(lev) == 1)
        stop("number of factors is >=2 -- a constant is modelled through '~1'")
      S2 <- paste0(S, "..factor..")
      if (extra <- exists(S2, envir=Env, inherits=FALSE)) {#factor may not
        ##                                         be part of a formula!
        C2 <- eval(parse(text=S2), envir=Env)
        lev2 <- levels(C2)
        if (length(lev2) != length(lev))
          ## it is assumed that the factors match. No control!!
          stop("a factor must have the same number of levels in both, the conditioning data set and the one for prediction")
      }      

      L <- list(RM_PLUS[1])
      i <- 2
      while (i <= length(lev)) {
        last <- min(length(lev), i + MAXSUB -1)
        plusList <- list(RM_PLUS[1])        
        while (i <= last) {
          model <- list(RM_COVARIATE)
          model[[COVARIATE_NAME_NAME]] <- paste0(S, i)
          model[[COVARIATE_EXTRA_DATA_NAME]] <- extra
          L1 <- as.numeric(C == lev[i]) ## 'data'
          if (extra) {
            dummy.matrix.fill <- Inf
            L2 <- as.numeric(C2 == lev2[i]) ## 'data'
            delta <- length(L1 - L2)
            L1 <- cbind(c(L1, rep(max(0, -delta))), c(L2, rep(max(0, delta))))
          }
          model[[COVARIATE_C_NAME]] <- L1  ## 'data'
          model[[COVARIATE_ADDNA_NAME]] <- add.na[k] > 1
          model[[COVARIATE_RAW_NAME]] <- TRUE 
          plusList[[length(plusList) + 1]] <- model
          i <- i + 1;
        }
        L[[length(L) + 1]] <- plusList
      }
      listModel[[k+1]] = if (length(lev) == 2) model
                         else if (length(lev) <= MAXSUB + 1) plusList
                         else  L
                  
    } else {
      if (isFormalCovModel[k]) { 
        P  <- buildCovList(C, Env=Env, EnvDummies=EnvDummies, add.na=add.na[k], 
                           unclear=unclear, EnvTest=EnvTest)
        model <- P$model
        if (abs(add.na[k]) > 1) {
          tL <- if ((model[[1]] %in% RM_MULT)) model[-1] else list(model) 
          names <- sapply(tL, function(x) x[[1]])
          if (all(names != R_CONST)) {## user given model does not start with
            ## R_CONST (or ~1 had not been transformed to R_CONST(NA)),add "NA*"
            const <- list(R_CONST)
            const[[CONST_A_NAME]] <- NA
            const[[COVARIATE_NAME_NAME]] <- S
            model <- c(SYMBOL_MULT, tL, list(const))
          }
        }
      } else { ## pure numbers (might be within matrix)
        if (length(C) > 1) {
          swaps <- 0 
          swap.position <- integer(1 + if (is.matrix(C)) ncol(C) / 2 else 0)
          swap.position[swaps + 1] <- 0
          if (extra <- is.matrix(C) && ncol(C)>1) { ## potentielle c, cbind
            if (length(n <- colnames(C)) > 0) {     ## abfangen
              w <- n == ".."
              step <- 2 * w - 1
              sum <- 0;
              for (i in 1:length(w)) {
                ## Ex: ---+++-+--++
##                Print(C, w, i, sum)
                if ((w[i] && sum >= 0) || (!w[i] && (sum != 0 && w[i-1])) )
                  stop("Data size mismatch. If you cannot find the error: ",
                       CONTACT)
                sum <- sum + step[i]
                if (sum == 0) {
                  swaps <- swaps + 1 
                  swap.position[swaps + 1] <- i
                }
              }
              if (sum != -length(w) && ## sum == -length(w) iff ".." unused
                  sum != 0 ## sum == 0 iff ".." used and correctly used
                  )
                stop("Data size mismatch. If you cannot find the error: ",
                     CONTACT)
              extra <- sum == 0 ## or swaps > 0
            }
          }
          for (L in 1:max(1, swaps)) {
            model <- list(RM_COVARIATE)         
            model[[COVARIATE_C_NAME]] <-
              if (swaps <= 1) C ## fuer L=1 && swaps=1 : if & else das Gleiche
              else C[, (swap.position[L] + 1):swap.position[L+1], drop=FALSE]
            model[[COVARIATE_ADDNA_NAME]] <- add.na[k] > 1
            model[[COVARIATE_EXTRA_DATA_NAME]] <- extra
            model[[COVARIATE_RAW_NAME]] <- TRUE 
            if (L < swaps) {
              model[[COVARIATE_NAME_NAME]] <- paste(S, L)
              listModel[[m]] <- model
              m <- m + 1
            } else if (swaps > 1) S <- paste0(S, L)
          }

          if (length(unclear) > 0)
            unclear <- detect.covariates(S=S, EnvTest=EnvTest,
                                         unclear=setdiff(unclear, rownames(C)))
        } else { ## a simple constant, could be a spatial constant mean
          model <- list(R_CONST)
          model[[CONST_A_NAME]] <-
            if (is.finite(C) && C == 1 && add.na[k]) NA else as.numeric(C)
        }
        model[[COVARIATE_NAME_NAME]] <- S
      } ## length(C) > 1
      listModel[[m]] <- model
    } ## not a factor
  } # for summand

  if (length(summands) == 1) listModel <- listModel[[2]]
  else listModel[[1]] <- SYMBOL_PLUS

  return(list(model=listModel, add.na=add.na, unclear=unclear))
}


buildCovList <- function(model, Env, EnvDummies, add.na, unclear, EnvTest) {
  ## Print("bCL", model)
  ## solves user defined model to a list by using
  ## in particular the models in RMmodels.R and the definition in
  ## RMmodelsBasics.R
  
  if (is.atomic(model) || is.list(model) ||
      is.language(model) || is.environment(model))
    return(list(model=model)) ## for recursive calling

  if (!isS4(model) || !is(model, CLASS_CLIST))
    stop('model must be of class ', CLASS_CLIST)

  P <- model@par.model
  if (length(P[[ADD_NA]]) > 0) {
    add.na <- P[[ADD_NA]]
    P[[ADD_NA]] <- NULL
  }

  name <- model@name
  if (name == RM_PLUS[1]) name <- SYMBOL_PLUS
  else if (name == RM_MULT[1]) name  <- SYMBOL_MULT
  else if (name == RM_COVARIATE) 
    unclear <- setdiff(unclear, P[[COVARIATE_NAME_NAME]])
  
  PD <- list(P=c(P, model@submodels),
             D=model@par.general[model@par.general != NO_DOLLAR_VALUE])
 
  for (j in 1:length(PD)) {
    P <- PD[[j]] 
    for (i in OneTo(length(P))) {
      X <- if (is(P[[i]], "formula")) 
             parseModel(P[[i]], Env=Env, EnvDummies=EnvDummies,
                        add.na = add.na, unclear=unclear,
                        EnvTest=EnvTest)
           else buildCovList(P[[i]], Env=Env, EnvDummies=EnvDummies,
                             add.na=add.na,unclear=unclear,
                             EnvTest=EnvTest)
      PD[[j]][[i]] <- X$model
      unclear <- X$unclear
    }
  }

  li <- c(list(name), PD$P)
  if (length(PD$D) > 0) li <- c(DOLLAR[1], PD$D, list(li))

  return(list(model=li, unclear=unclear))
}


removeParenthesis <- function(string) { ## Warum noch verwendet? 15.12.21??
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

##   library(RandomFields); c(RMexp()); c(RMexp(), RMexp())
