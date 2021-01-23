
## Authors 
## Martin Schlather, schlather@math.uni-mannheim.de
##
##
## Copyright (C) 2017 -- 2017 Martin Schlather
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 3
## of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied$x warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  



#   stop("big data sets currently not allowed") 
#    printlevel <- mle.methods <- lsq.methods <- recall <- sdvar <- general <- TRUE




## to do: grid
GetNeighbourhoods <- function(model, Z, X,
                              splitfactor, maxn, split_vec, shared=FALSE,
                              consider.correlations = TRUE ## ok for kriging
                              ## not ok for bigdata as parameters are NA!!
                             ) {

  model.nr <- MODEL_AUX
  rfInit(model= list("RFfctn", model), x = Z$coords, reg=model.nr,
         RFopt=getRFoptions(SAVEOPTIONS=NULL))  
  lc <- length(Z$coords)
  locfactor <- as.integer(splitfactor * 2 + 1) ## total number of
  ##    neighbouring boxes in each direction, including the current box itself
  maxn <- as.integer(maxn)        ## max number of points, including neighbours
  minimum <- as.integer(split_vec[NEIGHB_MIN+1])# min. nmbr of pts in neighbhood
  splitn <- as.integer(split_vec[NEIGHB_SPLIT+1]) # nmbr of locations when split
  maximum <- as.integer(split_vec[NEIGHB_MAX+1])## max nmbr of points when still
  ##                                    neighbours of neighbours are included.
  ##                         Note that, mostly, an additional box is included.
  xdimOZ <- Z$xdimOZ
  tsdim <- Z$tsdim
  newZ <- list()
 
  if (Z$dist.given) { ## to do
    ndata <- nrow(Z$data[[1]])
    stop("space splitting for 'distances' not programmed yet")
    if (is.vector(dist)) {
      j <- 1:ndata
      composite <- list()
      li <- 0
      maxi <-  (splitfactor[2] - 1) / Z$vdim 
        # mini <-  splitfactor[3] / Z$vdim
      while (length(j) >= maxi) {
        distpos <- (j[1]-1) * (ndata - j[1] / 2) + 1
          distidx <- distpos : (distpos + (ndata-j[1]-1))
        
        locs <- (j[1]+1):ndata
        locidx <- which(!is.na(pmatch(locs, j)))
        locidx <- c(locs[locidx], j[1])            
        li <- li + 1
        composites[[li]] <- locidx
        j <- j[is.na(pmatch(j, locidx))]
      }
      if (length(j) > 0) {
        li <- li + 1
        composite[[li]] <- j
      }
      
      ## kann jetzt noch optimiert werden hinsichtlich schaetzqualitaet durch
      ## zusammenfassen
    } else {
      stop("not programmed yet")
    }
    
    ## return(result)
    return(newZ)
  }


  newZ <- list()
  totpts <- sapply(Z$coords, function(x) x$totpts) ## passt das ? 2.2.19
  for (set in 1:lc) {
     ## n <- as.integer(dim(Z$coords[[set]])[2])
    coords <- Z$coords[[set]]
    data <- Z$data[[set]]

    oldlen <- length(newZ)
    nOT <- nT <- totpts[set]
    if (Z$has.time.comp) nOT <- nOT / coords$T[3]
    
    if (consider.correlations) { ## to rto      
      natsc <- .C(C_MultiDimRange, as.integer(model.nr),
                  as.integer(set),
                  natsc = double(tsdim)
                  )$natsc
    } else natsc <- rep(1, tsdim)  
    
     
    if (coords$grid) { ## given
      stop("space splitting for 'grid' not programmed yet")
      XX <- cbind(coords$x, coords$T)
      dim(data) <- c(XX[3, ], dim(data)[2])
      pts <- natsc / XX[2, ]
      total <- prod(pts)
      factor <-(total / splitn) ^ (1/tsdim)
      pts <- pts / factor
      blocks <- round(XX[3,] / pts)
      blocks[blocks < 1] <- 1
      old <- FALSE
      while(prod(ceiling(XX[3, ] / blocks)) > maximum) {
        stopifnot(any(!old))
        idx <- blocks == max(blocks[!old])        
        old <- old | idx
        blocks[idx] <- blocks[idx] + 1
      }
      minpts <- trunc(XX[3, ] / blocks)
      remaining <- XX[3, ] - minpts * blocks
      blocknpts <- (blocks - remaining) * minpts
      listseq <- combi <- list()
      for (i in 1:tsdim) {
        combi[[i]] <- c(-1, 1)
        listseq[[i]] <- 1:blocknpts[i]
      }      
      combi <- if (tsdim > 1) do.call(expand.grid, combi)
               else list(as.matrix(combi))
      for (j in 1:nrow(combi)) {
        L <- list(data)
        idx.combi <- combi[j, ] 
        for (i in 1:tsdim) L[[i+1]] <- idx.combi[i] * listseq[[i]]
        L[[length(L) + 1]] <- TRUE ## vdim + repet
     
        if (all(idx.combi > 0 | remaining[idx.combi < 0] > 0)) {          
          newZ[[length(newZ) + 1]] <- Z
          pts <- minpts + (idx.combi < 0)
          newZ[[length(newZ)]]$coords$x[3, ] <- pts
          newZ[[length(newZ)]]$data <- matrix(nrow=prod(pts), do.call("[", L))
        }
      }
    } else { # !coords$grid
      maxnOT <- maxn
      nDsplitn <- nT / splitn
      x <- t(coords$x)

      u <- apply(coords$x, 2, function(x) length(unique(x)))
      Range <- apply(coords$x, 2, range)
      
      if (Z$has.time.comp) {
        T <- coords$T
        u <- c(u, T[XLENGTH + 1])
        Range <- cbind(Range, c(0,(T[XLENGTH + 1]-1)) * T[XSTEP+1] +
                              T[XSTART + 1])
        Tseq <- seq(0:(T[XLENGTH + 1]-1)) * T[XSTEP + 1] + T[XSTART + 1]
        x <- rbind(x, 0)  ## dummy to get the dimensins right!
      }
      
      rd <- apply(Range, 2, diff)  ### incl T 
      len <- pmax(1e-10 * max(rd), rd * (1 + 1e-10))### incl T 
      units <- pmax(1, len * natsc) ### incl T

#      print(Range)
 #     Print(rd, len, units, Range); kkkk
            
      ## * "gerechte" Aufteilung in alle Richtungen waere nDsplitn
      ## * die Richtung in die viele units sind, soll eher aufgespalten werden
      ## * ebenso : wo viele Werte sind eher aufspalten
      blockidx <- (nDsplitn / prod(units * u))^{1/tsdim} * units * u > 0.5
      reddim <- sum(blockidx)
      units <- units[blockidx]
      zaehler <- 1
      blocks <- rep(1, tsdim)
      OK <- integer(1)
      
      repeat {
        blocks[blockidx] <- (nDsplitn / prod(units))^{1/reddim} *
          locfactor * zaehler * units * Z$vdim

#        Print(blocks, blockidx, nDsplitn , prod(units), {1/reddim} ,
#          locfactor * zaehler * units * Z$vdim, x, Range[1,], len)
        
        blocks <- as.integer(ceiling(blocks))
        coord.idx <- floor((x - Range[1,]) / (len / blocks))
        if (Z$has.time.comp) coord.idx[nrow(coord.idx), ] <- 0
        cumblocks <- cumprod(blocks) ### alles bis hier inkl. T        
        totblocksOT <- as.integer(cumblocks[xdimOZ])

        Ccumblocks <- as.integer(c(1, cumblocks))
        cumblocks <- Ccumblocks[-length(Ccumblocks)]
 
        if (Z$has.time.comp && blocks[tsdim] > 1) {
          maxnOT <- as.integer(maxn / ceiling(T[3] / blocks[tsdim]))
        }
        
        ## zuordnung der coordinaten_Werte zu den jeweiligen "blocks"
        ## given ist liegend
        cumidx <- as.integer(colSums(coord.idx * cumblocks))

#        Print(cumidx, coord.idx, Range, x, x - Range[1,], cumblocks, len, blocks)
 #       str(cumidx)
  #      xxx
   #     str(nOT)
    #    str(totblocksOT)
        nOT <- as.integer(nOT)
        elms.in.boxes <- .Call(C_countelements, cumidx, nOT, totblocksOT)

        neighbours <- .Call(C_countneighbours, xdimOZ, blocks, locfactor,
                            Ccumblocks, elms.in.boxes, maxnOT)
        ## if there too many points within all the neighbours, then split
        ## into smaller boxes
        zaehler <- zaehler * 2
      
        ## image(neighbours, zlim=c(0:(prod(blocks)-1)))
        if (!is.null(neighbours)) break;
      } # repeat

      l <- list()

  ###    Print(cumidx, xdimOZ, nOT, Ccumblocks, elms.in.boxes)
      
      l[[1]] <- .Call(C_getelements, cumidx, xdimOZ, nOT, Ccumblocks,
                      elms.in.boxes)


      
      l1len <- sapply(l[[1]], length)

      if (!hasArg("X")) {
        ## Print(cumidx, xdimOZ, nOT, Ccumblocks, elms.in.boxes)
        
	L <- .Call(C_getelements, cumidx, xdimOZ, nOT, Ccumblocks,
		   elms.in.boxes)
        for (idx in L) {
          ##Print(coords, idx, data)
          new.x <-coords$x[idx, , drop=FALSE]
          
          if (Z$has.time.comp && blocks[tsdim] > 1) {
            minpts <- trunc(T[3] / blocks[tsdim])
            remaining <- T[3] - minpts * blocks
            blocknpts <- (blocks - remaining) * minpts
            combi <- c(-1, 1)

            ##Print("this is wrong")

            
            listseq <- 1:blocknpts

            
            
            new.data <- data
            ## Print(nrow(data) / T[3], T[3], nrow(new.data), new.data, Z$data[[set]], data, idx)
            dim(new.data) <- c(nrow(data) / T[XLENGTH + 1], T[XLENGTH + 1],
                               ncol(new.data))
            for (j in 1:length(combi)) {
              idx.combi <- combi[j] 
              if (idx.combi >0 || remaining[idx.combi < 0] > 0) {
                ## newz wird z.Zt nicht weiter verwendet
                newZ[[length(newZ) + 1]] <- Z
                newZ[[length(newZ)]]$coords$x <- new.x
                pts <- minpts + (idx.combi < 0)
                newZ[[length(newZ)]]$coords$T[3] <- pts
                newZ[[length(newZ)]]$data <-
                  data.matrix(new.data[ ,idx.combi * listseq, ])
                n <- length(newZ[[length(newZ)]]$data)
                nr <- length(idx) * pts
                dim(newZ[[length(newZ)]]$data) <- c(nr, n / nr)
              }
            }
          } else {
            newZ[[length(newZ) + 1]] <- Z
            newZ[[length(newZ)]]$coords$x <- new.x
            newZ[[length(newZ)]]$data <- data[idx, , drop=FALSE]  ## Z$data[idx, ]
          }
        } # for idx        
      }
      
    }# !coords$grid
    
    if (hasArg(X)) {
 ##      l1len <- sapply(l[[1]], length)
      if (X$grid) {
        stop("not programmed yet")
      } else { 
        ## now calculate the boxes for the locations where we will interpolate
        nX <- nrow(X$x) ## how many new locations?
        tX <- t(X$x)
        i <- pmax(0, pmin(blocks-1, floor((tX - Range[1,]) / (len / blocks))))
        dim(i) <- dim(tX)
        i <- as.integer(colSums(i * cumblocks))
      
        res.in.boxes <- .Call(C_countelements, i, nX, totblocksOT)
        
        notzeros <- res.in.boxes > 0

   #     Print(i, xdimOZ, nX, Ccumblocks, res.in.boxes)
        
        l[[3]] <-
          .Call(C_getelements, i, tsdim, ## tsdim, da Einzelpunkte!
                nX, Ccumblocks, res.in.boxes)[notzeros]
        ## TO DO : idx[[3]] passt nicht, da sowohl fuer Daten
        ##          als auch coordinaten verwendet wird. Bei repet > 1
        ##         ist da ein Problem -- ueberpruefen ob repet=1
      }
    } else {
      notzeros <- TRUE
    }

        
    ll <- .Call(C_getneighbours, xdimOZ, ## ts.xdim
                blocks, locfactor, Ccumblocks,
                neighbours)[notzeros]
    less <-
      sapply(ll, function(x) sum(elms.in.boxes[x]) < minimum) | !shared
    ##                  if !shared then all(less)==TRUE
      
    if (any(less)) {
      not.considered.yet <- sapply(l[[1]], length) > 0   
      newll <- ll
      for (i in which(less)) {
        current <- ll[[i]]
        elements <- sum(elms.in.boxes[current] *
                        (shared | not.considered.yet[current]))# number of pts in a neighbourhood
        while (elements < minimum) {
          new <- unique(unlist(ll[current])) # neighbours of neighbours, but not
          new <- new[which(is.na(pmatch(new, current)))]# neighbours themselves
          nn <- elms.in.boxes[new] * (shared | not.considered.yet[new]) # how many pts are in each of these boxes?
          ordr <- order(nn)
          new <- new[ordr]
          nn <- nn[ordr]
          cs <- elements + cumsum(nn)
          smaller <- sum(cs <= maximum) ## now, check which neighbours of
          ## the neigbours can be included in the list of neighbours of i
          ## to increase the number of points in the kriging neighbourhood
          if (smaller == 0) break; ## none
          if (smaller == length(cs) || cs[smaller] >= minimum ||
              cs[smaller+1] > maxn) {
            if ( (elements <- cs[length(cs)]) <= maxn ) {            
              current <- c(current, new)            
            } else {
              current <- c(current, new[1:smaller])
              elements <- cs[smaller]
            }
            if (smaller != length(cs)) break
          } else {
            ## smaller < length(cs) && cs[smaller] < minimum &&
            ## cs[smaller+1]<=maxn
            ## i.e., include the next one, but there is no chance to include
            ## more within the rules.
            elements <- cs[smaller+1]
            current <- c(current, new[1:(smaller+1)])
            break;
          }
        }
        current <- current[l1len[current] > 0]
        if (!shared) current <- current[not.considered.yet[current]]
        newll[[i]] <- current
        not.considered.yet[current] <- FALSE                            
      }
      newll <- newll[sapply(newll, length) > 0]
      l[[2]] <- newll
    } else l[[2]] <- ll
#  } ## locations to be estimated not on grid
#} ## locations to be estimated
  } ## for

  return(if (shared) l else lapply(l[[2]], function(x) unlist(l[[1]][x])))


  if (hasArg("X")) {
    return(if (shared) l else lapply(l[[2]], function(x) unlist(l[[1]][x])))
  } else { ## to do
    return(if (shared) l else lapply(l[[2]], function(x) unlist(l[[1]][x])))
    return(newZ)
  }
  
}

GetComposites <- function(Z, cliquesize) {
  stopifnot(cliquesize == 2)
  return(Z)
}


BigDataSplit <- function(Z, RFopt) {
  fit <- RFopt$fit
  method <- fit$likelihood   
  if (is.na(method <- pmatch(method, RC_LIKELIHOOD_NAMES)))
    stop("unknown value for 'likelihood'.")
  method <- RC_LIKELIHOOD_NAMES[method] # kein + 1 notwendig

  totpts <- sapply(Z$coords, function(x) x$totpts)
  
  if (method == "full" ||
      (method %in% c("auto", "tesselation") && all(totpts<=fit$max_neighbours)
       )) return(Z)
  if (method == "auto") { ## default
    method <- "tesselation"
    if (RFopt$basic$printlevel>=PL_IMPORTANT)
      message("Too many locations to use standard estimation methods.\n",
              "Hence an approximative methods is used. However, it is *not* ",
              "ensured\n",
              "that they will do well in detecting long memory dependencies.")
  }

  if (method == "tesselation") {
    if (any(diff(fit$cliquesize) <= 0)) stop("in case of 'tesselation', 'cliquesize' must contain three increasing numbers")
    
    return(GetNeighbourhoods(Z$model, Z=Z, 
                             splitfactor=fit$splitfactor_neighbours,
                             maxn=fit$max_neighbours,
                             split_vec=fit$cliquesize,
                             consider.correlations = FALSE,# should be improved!
                             shared=FALSE)
           )
  } else if (method == "composite") {
    if (any(diff(fit$cliquesize) != 0)) stop("in case of 'composite', 'cliquesize' must be a single value")
    return(GetComposites(Z=Z, cliquesize = fit$cliquesize))
  } else stop("unknown 'likelihood' value")  
}
