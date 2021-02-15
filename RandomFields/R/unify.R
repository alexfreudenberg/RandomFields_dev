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


GetRanges <- function(Z, mindist_pts) {
  ## mindist_pts = RFopt$fit$smalldataset / 2
  if (length(Z$coords) > 0) Z <- Z$coords## could be just a list of distances,
  ## see UnifyXT
  else Z <- list(Z)
  ans  <- lapply(Z, function(z) {
    if (length(z$dist.given) == 0 || z$dist.given) {
      if (is.matrix(z$x)){
        r <- apply(z$x, 2, range)
        rS <- rowSums(z$x^2)
        min <- min(rS)
        max <- max(rS)
      } else {
        r <- range(z$x)
        min <- r[1]^2
        max <- r[2]^2
      }
      list(min=min, max=max, r=r)
    } else {
      if (z$grid) {
        d <- z$x[XSTEP+1, ] * (z$x[XLENGTH+1, ] - 1)
        min <- sum(z$x[XSTEP+1, ]^2)
        max <- sum(d^2)
        r <- rbind(z$x[XSTART + 1, ], z$x[XSTART+1,] + d)
      } else {
        if (nrow(z$x) >= 2)
          d <- if (nrow(z$x) <= mindist_pts) dist(z$x)
               else dist(z$x[sample(nrow(z$x), mindist_pts), ])
        min <- if (nrow(z$x) < 2) Inf else min(d)^2
        max <- if (nrow(z$x) < 2) 0 else max(d)^2
        r <- apply(z$x, 2, range)
      }
      if (z$has.time.comp) {
        d <- z$T[XSTEP+1] * (z$T[XLENGTH+1]-1)
        min <- min + z$T[XSTEP+1]^2
        max <- max + d^2
        r <- cbind(r, c(z$T[XSTART + 1], z$T[XSTART+1] + d))
      }
      list(min=min, max=max, r=r)
    }
  })

  rangex <- sapply(ans, function(x) x$r)
  dim <- length(rangex) / length(Z) / 2
  base::dim(rangex) <- c(2, dim, length(Z))
  
  list(rangex=apply(rangex, 2, range),
       mindist=sqrt(min(sapply(ans, function(x) x$min))),
       maxdist=sqrt(max(sapply(ans, function(x) x$max))) )
}
  
 

seq2grid <- function(x, name, grid, warn_ambiguous, gridtolerance) {
  xx <- matrix(nrow=3, ncol=length(x))
  step0 <- rep(FALSE, length(x))
  gridnotgiven <- missing(grid) || length(grid) == 0
  
  for (i in 1:length(x)) {
    if (length(x[[i]]) == 1) {
      xx[,i] <- c(x[[i]], 1, 1)
      next
    }
    step <- diff(x[[i]])
    if (step[1] == 0.0) {
      
      ok <- step0[i] <- all(step == 0.0)      
    } else {
      ok <- max(abs(step / step[1] - 1.0)) <= gridtolerance
    }

    if (!ok) {
      if (gridnotgiven) return(FALSE)
      if (!TRUE)
        Print(i, x[[i]][1:min(100, length(x[[i]]))], #
              step[1:min(100,length(step))],
              range(diff(step[1:min(100,length(step))])))
      stop("Different grid distances detected, but the grid must ",
           "have equal distances in each direction -- if sure that ",
           "it is a grid, increase the value of 'gridtolerance' which equals ",
           gridtolerance,".\n")
    }

    xx[,i] <- c(x[[i]][1], step[1], if (step0[i]) 1 else length(x[[i]]))
  }

#  if (FALSE && gridnotgiven &&  length(x) > 1) Warning("ambiguous","Ambiguous interpretation of coordinates. Better give 'grid=TRUE' explicitly. (This message appears only once per session.)" )
 
  if (any(step0)) {
    if (all(step0)) {
      if (gridnotgiven) return(FALSE)
      else stop("Within a grid, the coordinates must be distinguishable")
    } else if (gridnotgiven) Note("ambiguous")
  }

  return(xx)
}

RFearth2cartesian <- function(coords, units=NULL, system = "cartesian",
                              grid=FALSE) {
  ## may not be called internally !!
  RFopt <- internalRFoptions(getoptions_="coords")
  if (is.character(system)) system <- pmatch(system, ISO_NAMES) - 1
  stopifnot(system %in%
            c(CARTESIAN_COORD, GNOMONIC_PROJ, ORTHOGRAPHIC_PROJ))
  if (missing(units) || is.null(units)) {
    global.units <- RFopt$new_coordunits[1]
    units <-if (global.units[1] == "") "km" else global.units
  }
  if (!is.matrix(coords)) coords <- t(coords)
  res <- RFfctn(COPY=FALSE, # OK
                RMtrafo(new=system), coords, grid=grid,
                coords.new_coord_system = "keep",
                coords.new_coordunits=units,
                coords.coord_system="earth")
  dimnames(res) <- list(NULL, c("X", "Y", "Z", "T")[1:ncol(res)])
  return(res)
}

RFearth2dist <- function(coords, units=NULL, system="cartesian",
                         grid=FALSE, ...) {
  ## may not be called internally !!
  RFopt <- internalRFoptions(getoptions_="coords")
 if (is.character(system)) system <- pmatch(system, ISO_NAMES) - 1
  stopifnot(system %in%
            c(CARTESIAN_COORD, GNOMONIC_PROJ, ORTHOGRAPHIC_PROJ))
  if (missing(units) || is.null(units)) {
    global.units <- RFopt$new_coordunits[1]
    units <- if (global.units[1] == "") "km" else global.units 
  }
  if (!is.matrix(coords)) coords <- t(coords)
  z <- RFfctn(COPY=FALSE, # OK
              RMtrafo(new=system), coords, grid=grid,
                coords.new_coord_system = "keep",
                coords.new_coordunits=units,
                coords.coord_system="earth")
  return(dist(z, ...))
}

UnifyXT <- function(x=NULL, y=NULL, z=NULL, T=NULL, grid,
                    distances=NULL,
                    dim=NULL, # == spatialdim!
                    length.data,
                    y.ok = FALSE, 
                    Ty=T, gridy){
#  Print("unify", x)
  ## do not pass anything on "..." ! --- only used for internal calls
  ## when lists are re-passed

  ## converts the given coordinates into standard formats
  ## (one for arbitrarily given locations and one for grid points)
 
  RFopt <- getRFoptions()
  PL <- RFopt$basic$printlevel
  duplicated <- which(RFopt$general$duplicated_loc == DUPLICATEDLOC_NAMES) - 1

  if (!missing(x) && !is.null(x) ) {
    if (is(x, "UnifyXT")) {
      if (length(y) > 0 || length(z) > 0 || length(T) >0)
        stop("'y', 'z', 'T', may not begiven if class(x) is 'Unify")
      return(x)
    }
    if (is.list(x) && !is.data.frame(x)) {
      ##      Print(!is.list(x[[1]]), x,length(x), is.data.frame(x[[1]]))
      if (!is.list(x[[1]])) return(do.call("UnifyXT", x)) # simple list
      ## nested list
      L <- lapply(x, function(xx)
        if (is(xx, "UnifyXT")) xx
        else if (is.data.frame(xx)) UnifyXT(xx)
        else do.call("UnifyXT", xx))
      if (length(L) > 1) {    
        if (!all(diff(sapply(L, function(x) x$has.time.comp)) == 0) ||
          !all(diff(sapply(L, function(x) x$spatialdim)) == 0))
          stop("all sets must have the same dimension")
        if (!all(diff(sapply(L, function(x) x$dist.given)) == 0))
          stop("either all the sets must be based on distances or none")        
      } ## else L[[1]]
      class(L) <- "UnifyXT"
      return(L)
    }
  }

  if (!missing(gridy) && length(y) > 0 && y.ok) {
    L <- UnifyXT(x=x, T=T, grid=grid)
    y <- UnifyXT(x=y, T=Ty, grid=gridy)
    if (L$spatialdim != y$spatialdim)
      stop("'x' and 'y' do not indicate the same spatial dimension")
    if (xor(length(L$T) == 0, length(y$T) == 0))
      stop("Time must be given for both 'x' and 'y', or for none")
    L$y <- y$x
    L$Ty <- y$T
    L$gridy <- y$grid
    L$totptsy <- as.integer(y$totpts)
    return(L)
  } else if (!missing(gridy))
    stop("'gridy' is set although it should not", CONTACT)
  
  curunits <- RFopt$coords$coordunits
  newunits <-  RFopt$coords$new_coordunits
  coord_system <-  RFopt$coords$coord_system
  new_coord_system <-  RFopt$coords$new_coord_system
  ex.red <- RFopt$internal$examples_reduced
  change.of.units <- FALSE
  
  if (!missing(distances) && !is.null(distances)) { ## length==0 OK!    
    stopifnot(is.matrix(distances) || is(distances, "dist") ||
                                       is.vector(distances),
              (!missing(dim) && !is.null(dim)),
              (missing(grid) || length(grid) == 0),
              missing(x) || is.null(x),
              length(y)==0,
              length(z)==0,
              length(T)==0)
    
    if (coord_system != new_coord_system && new_coord_system != "keep")      
      stop("coordinate systems differ")
    
    if (is.list(distances)) {
      L <- list()
      for (i in 1:length(distances))
        L[[i]] <- do.call("UnifyXT", list(distances=distances[[i]], dim=dim))
      class(L) <- "UnifyXT"
      return(L)
    }
    
    if (is(distances, "dist")) {
      x <- as.vector(distances)
      len <- length(distances)
    } else if (is.matrix(distances) || is.vector(distances)) {
      if (is.matrix(distances)) {        
        len <- nrow(distances)
        if (missing(dim) || is.null(dim)) dim = ncol(distances)
        else if (dim != ncol(distances))
          stop("matrix of distances does not fit the given dimension")
        ## see below!
        if (duplicated != DUPLICATEDLOC_RISKERROR &&
            duplicated != DUPLICATEDLOC_REPEATED **
            any(rowSums(distances != 0) == 0))
          stop("distances must be all greater than zero")
     } else {
        len <- length(distances)
        if (missing(dim) || is.null(dim))
          stop("dim is not given although 'distances' are used")
        ## NOTE!
        ## SpatialNugget in nugget.cc uses that this condition is satisfied!
        if (!all(distances > 0)) stop("distances must be all greater than zero")
       }
      x <- distances
    } else {
      stop("'distances' not of required format.")
    }

    if (duplicated != DUPLICATEDLOC_RISKERROR &&
        duplicated != DUPLICATEDLOC_REPEATED) {
      mindist <- GetRanges(distances, RFopt$fit$smalldataset / 2)$mindist
      if (mindist <= RFopt$nugget$tol) {
        if (duplicated == DUPLICATEDLOC_ERROR)
          stop("Some locations are not distinguishable. If this is intentional, see ?RPnugget and ?RFoptions 'duplicated_locations' to deal with repeated measurements.")
        mindist <- (1e-6 * (RFopt$nugget$tol == 0) + 0.2 * RFopt$nugget$tol) /
          sqrt(spatialdim)
        if (is.vector(distances)) distances[distances == 0] <- mindist
        else distances[rowSums(abs(distances)) == 0, ] <- mindist
      }
    }

    if (ex.red && len > ex.red^2 / 2) {
      totpts <- as.integer(ex.red)
      len <- as.integer(totpts * (totpts - 1) / 2)
      x <- if (is.matrix(x)) x[1:len ,] else x[1:len]
    } else {
      totpts <- as.integer(1e-9 + 0.5 * (1 + sqrt(1 + 8 * len)))
      if (totpts * (totpts-1) / 2 != len)
        stop("number of distances not of the form 'n (n - 1) / 2'.")     
    }

    if (duplicated == DUPLICATEDLOC_ERROR && anyDuplicated(x))
         stop("Some locations are not distinguishable. If this is intentional, see ?RPnugget and ?RFoptions 'duplicated_locations' to deal with repeated measurements. First duplicated is ", anyDuplicated(x))


    ## keep exactly the sequence up to 'distances'
    if (storage.mode(x) != "double") storage.mode(x) <- "double"
    L <- list(x = as.matrix(x), #0
              T= double(0),  #1
              grid = FALSE, #2
              spatialdim=as.integer(dim),#3
              has.time.comp=FALSE, #4
              dist.given = TRUE, #5
              y = double(0),   #6
              Ty= double(0),  #7
              gridy = FALSE, #8
              totpts = as.integer(totpts), ## 9 number of points
            #  l = len, ## 10, ?? physical length??
              totptsy = 0L, ## 11 number of points
          #    ly= 0, ## 12,?? physical length??
              coordunits = rep(curunits, length=dim), #13
              new_coordunits = rep(newunits, length=dim) #14
              )
    class(L) <- "UnifyXT"
    return(L)
  }
  
                                        #stopifnot(!missing(x))
  if (is(x, "RFsp") || isSpObj(x)) {
    return(UnifyXT(x=coordinates(x),
                   y=y,
                   z=z, T=T, grid=grid,
                   distances=FALSE, length.data=length.data,
                   y.ok=y.ok))
  }    
  if (is.raster(x)) x <- as(x, 'GridTopology')
  
  if ((missing(grid) || length(grid) == 0) && !missing(length.data)) {
    new <-  Try(UnifyXT(x=x, y=y, z=z, T=T, grid=TRUE,
                        distances=FALSE, dim=if (!missing(dim)) dim,
                        length.data = length.data, y.ok = y.ok))
    if (grid <- !is(new, CLASS_TRYERROR)) {
      ratio <- length.data / new$totpts

      if (grid <- ratio == as.integer(ratio)) {
        if (PL>=PL_IMPORTANT && new$spatialdim > 1)
          message("Grid detected. If it is not a grid, set grid=FALSE. (To avoid this message set 'grid'.)\n")
      }
    }
    return(if (grid) new
           else UnifyXT(x=x, y=y, z=z, T=T, grid=FALSE, distances = FALSE,
                        dim = if (!missing(dim)) dim,
                        length.data = length.data, y.ok = y.ok))
  } # if (missing(grid) && !missing(length.data))


  gridtriple <- FALSE

  if (is.GridTopology <- is(x, "GridTopology")){
    x <- rbind(x@cellcentre.offset,
               x@cellsize,
               x@cells.dim)
    if ((missing(grid) || length(grid) == 0)) grid <- TRUE else stopifnot(grid)
    gridtriple <- TRUE
  }
  ##else {
  ##  is.GridTopology <- FALSE
  ##}

  
  if (is.data.frame(x)) {
    if (ncol(x)==1) x <- as.vector(x) else x <- as.matrix(x)
  }

  stopifnot(length(x) != 0)
                                        #  stopifnot(all(unlist(lapply(as.list(x), FUN=function(li) is.numeric(li))))) ## wann benoetigt???

  stopifnot(is.numeric(x))# um RFsimulte(model, data) statt data=data abzufangen
  
  
                                        #  stopifnot(all(is.finite(x)), all(is.finite(y)), all(is.finite(z))) ; s.u. unlist

  ## from here on it is assumed that it is a standard input,
  ## i.e. (x,y,z,T) define one set of coordinates
  ## OR: x and y give the coordinates for a kernel, so they have the
  ##     same structure. Especially, x and y must be matrices and
  ##     dim(x) == dim(y)
  if (is.matrix(x)) {
    if (!is.numeric(x)) stop("x is not numeric.")
    if (nrow(x) == 1) {
      grid <- gridtriple <- FALSE
    } else if ((d <- ncol(x)) > 1 &&
               x[1, ncol(x)] == x[2, ncol(x)])  ## fast check whether
      ##    expanded grid currentlyl, for y not programmed yet
    {

      dims <- integer(d)
      start <- step <- double(d)
      for (i in 1:d) {
        u <- sort(unique(x[, i]))
        start[i] <- u[1]
        dims[i] <- length(u)
        diffu <- diff(u)
        step[i] = if (all(diff(diffu) < 1e-14)) diffu[1] else NA
      }
      if (all(is.finite(step)) && prod(dims) == nrow(x)) {
        co <- vector("list", d)
        for (i in 1:d) co[[i]] <- seq(start[i], by=step[i], length.out=dims[i])
        coo <- as.matrix(do.call("expand.grid", co))
        if (all.equal(x, coo)) {
          x <- matrix(c(start, step, dims), byrow=TRUE, nrow=3)
          grid <- gridtriple <- TRUE
        }
      }
    }
   
    if (length(z)!=0) stop("If 'x' is a matrix, then 'z' may not be given")
    if (length(y)!=0) {
      if (!y.ok) stop("If 'x' is a matrix, then 'y' may not be given")
      if (length(T)!=0)
        stop("If 'x' is a matrix and 'y' is given, then 'T' may not be given")

      if (!is.matrix(y) || any(base::dim(y) != base::dim(x)))
        stop("'y' does not match 'x' (they must be matrices of equal size)")
    }
    
   change.of.units <- coord_system == COORD_SYS_NAMES[coord_auto + 1] &&
      ncol(x) >= 2 && ncol(x) <= 3 && !is.null(n <- dimnames(x)[[2]])

##    Print(RFopt, change.of.units, coord_system, x)
     
    if (change.of.units) {
      if (any(idx <- earth_coord_names(names=n, opt=RFopt$coords))) {
        if (length(idx) == 2 && !all(idx == 1:2))
          stop("earth coordinates not in order longitude/latitude")
        cur <- curunits[1]
        newunits <- RFopt$coords$new_coordunits
        curunits <- RFopt$coords$coordunits
        curunits[1:2] <- RFopt$coords$earth_coord_names[1:2]
        if (newunits[1] == "") newunits[1] <-  UNITS_NAMES[units_km + 1]
        newunits[2:3] <- newunits[1]
        Note("coordinates", cur=cur, units=newunits[1])
        coord_system <- COORD_SYS_NAMES[earth + 1]
        setRFoptions(coords.coord_system = coord_system,
                     coords.coordunits = curunits,
                     coords.new_coordunits = newunits,
                     messages.note_coordinates=FALSE)
      } else {
        setRFoptions(coords.coord_system =  COORD_SYS_NAMES[cartesian + 1])
      }
    }    
  
    spatialdim <- ncol(x)
    len <- nrow(x)
    if (spatialdim==1 && len != 3 && (missing(grid) || length(grid) == 0)) {
      if (length(x) <= 2) grid <- TRUE
      else {
        dx <- diff(x)
        grid <- max(abs(diff(dx))) < dx[1] * RFopt$general$gridtolerance
      }
    } 

    if ((missing(grid) || length(grid) == 0) &&
        any(apply(x, 2, function(z) (length(z) <= 2) || max(abs(diff(diff(z))))
                  > RFopt$general$gridtolerance))) {
      grid <- FALSE
    }

    if ((missing(grid) || length(grid) == 0) || !is.logical(grid)) {
      grid <- TRUE
      if (spatialdim > 1) Warning("ambiguous")
    }

    if (grid && !is.GridTopology) {
      if (gridtriple <- len==3) {
        if (PL >= PL_SUBIMPORTANT && RFopt$messages$warn_oldstyle) {
          message("Grid coordinates recognized of size ",
                  paste0("seq(", x[XSTART+1,], ", by=", x[XSTEP +1,],
                         ", len=", x[XLENGTH + 1,], ")",
                         collapse =" x "), ".")
        }
      } else len <- rep(len, times=spatialdim)   # Alex 8.10.2011
    }

    if (grid && !gridtriple) {
      ## list with columns as list elements -- easier way to
      ## do it??
      x <- lapply(apply(x, 2, list), function(r) r[[1]])   
      if (length(y) != 0) y <- lapply(apply(y, 2, list), function(r) r[[1]])
    }    
  } else { ## x, y, z given separately (x not a matrix)
    if (length(y)==0 && length(z)!=0) stop("y is not given, but z")
    xyzT <- list(x=if (!missing(x) && !is.null(x)) x, y=y, z=z, T=T)
    for (i in 1:4) {
      if (!is.null(xyzT[[i]]) && !is.numeric(xyzT[[i]])) {
        if (PL>PL_IMPORTANT) 
          message(names(xyzT)[i],
                  " not being numeric it is converted to numeric")
        assign(names(xyzT)[i], as.numeric(xyzT[[i]]))
      }
    }
    remove(xyzT)
    spatialdim <- 1 + (length(y)!=0) + (length(z)!=0)
    if (spatialdim==1 && ((missing(grid) || length(grid) == 0) || !grid)) {
      ## ueberschreibt Einstellung des Nutzers im Falle d=1
      if (length(x) <= 2) newgrid <- TRUE
      else {
        dx <- diff(x)
        newgrid <- max(abs(diff(dx))) < dx[1] * RFopt$general$gridtolerance
      }
      if ((missing(grid) || length(grid) == 0)) grid <- newgrid
      else if (xor(newgrid, grid)) Warning("on_grid", grid)
    }
    len <- c(length(x), length(y), length(z))[1:spatialdim]
    
    if (!(missing(grid) || length(grid) == 0) && !grid) { ## sicher nicht grid, ansonsten ausprobieren
      if (any(diff(len) != 0)) stop("some of x, y, z differ in length")
      x <- cbind(x, y, z)
      ## make a matrix out of the list
      len <- len[1]
    } else {
      if ((missing(grid) || length(grid) == 0) && any(len != len[1]))
        grid <- TRUE
      x <- list(x, y, z)[1:spatialdim]
    }
    y <- z <- NULL ## wichtig dass y = NULL ist, da unten die Abfrage
  }  ## end of x, y, z given separately 
  
  if (!all(is.finite(unlist(x)))) {
    stop("coordinates are not all finite")
  }


  if ((missing(grid) || length(grid) == 0) || grid) {
    if (gridtriple) { ## gridtriple implies grid
      if (len != 3)
        stop("In case of simulating a grid with option gridtriple, exactly 3 numbers are needed for each direction")
      totpts <- prod(x[XLENGTH + 1,])     
      if (length(y)!=0) {
        if (!all(y[XLENGTH + 1,] == x[XLENGTH + 1,]))
          stop("the grids of x and y do not match ")
      }
    } else {     
      xx <- seq2grid(x, "x",  grid,
                     RFopt$messages$warn_ambiguous, RFopt$general$gridtolerance)

      if (length(y)!=0) {
        yy <- seq2grid(y, "y", grid,
                       RFopt$messages$warn_ambiguous,
                       RFopt$general$gridtolerance)
        if (xor(is.logical(xx), is.logical(yy)) ||
            (!is.logical(xx) && !all(yy[XLENGTH + 1,] == xx[XLENGTH + 1,])))
          stop("the grids for x and y do not match")      
      }
      if (missing(grid) || length(grid) == 0) grid <- !is.logical(xx)       
      if (grid) {
        x <- xx
        totpts <- prod(len)
        len <- 3
        if (length(y) != 0) y <- yy
      } else {
        x <- sapply(x, function(z) z)
        if (length(y) != 0) y <- sapply(y, function(z) z)
      }
    }

    if (grid) {
      if (any(x[XLENGTH + 1, ] <= 0))
        stop(paste("length of the grid must be postive. Got as length",
                   paste(x[XLENGTH + 1,], collapse=",")))
      if (any(x[XSTEP + 1, ] <= 0)) {
        idx <- (x[XSTEP + 1, ] == 0) | duplicated <= DUPLICATEDLOC_LASTERROR
       if (any(x[XSTEP + 1, ] < 0) || (any(idx & x[XLENGTH + 1, idx] != 1)))
          stop(paste("step of the grid must be postive. Got as step",
                     paste(x[XSTEP + 1,], collapse=",")))
      }
    }
  }
  
  if (!grid) {
    totpts <- nrow(x)
    if (length(y)==0 && duplicated == DUPLICATEDLOC_ERROR) {
      if (totpts < 2000 && anyDuplicated(x)) 
         stop("Some locations are not distinguishable. If this is intentional, see ?RPnugget and ?RFoptions 'duplicated_locations' to deal with repeated measurements. First duplicated is ", anyDuplicated(x))
      }
      ## fuer hoehere Werte con total ist ueberpruefung nicht mehr praktikabel
  }
  
  if (coord_system == "earth") {
                                        # if (ncol(x) > 4) stop("earth coordinates have maximal 3 components")
    opt <- getRFoptions(getoptions_="coords") ## muss nochmals neu sein
    global.units <- opt$new_coordunits[1]
    if (global.units[1] == "") global.units <- "km"
    
    Raumdim <- ncol(x) #if (grid) ncol(x) else
    new_is_cartesian <- new_coord_system %in% CARTESIAN_SYS_NAMES
    if (new_is_cartesian) {
      if (sum(idx <- is.na(opt$zenit))) {
        zenit <- (if (grid) x[XSTART+1, 1:2] +
                            x[XSTEP +1, 1:2] * (x[XLENGTH + 1, 1:2] - 1)
                  else if (opt$zenit[!idx] == 1) colMeans(x[, 1:2])
                  else if (opt$zenit[!idx] == Inf) colMeans(apply(x[, 1:2], 2,
                                                                  range))
                  else stop("unknown value of zenit"))
        setRFoptions(zenit = zenit)
      }

      code <- switch(new_coord_system,
                     "cartesian" = CARTESIAN_COORD,
                     "gnomonic" = GNOMONIC_PROJ,
                     "orthographic" = ORTHOGRAPHIC_PROJ,
                     stop("unknown projection method")
                     )
      message("New coordinate system: ", new_coord_system, ".\n")
      x <- RFfctn(COPY=FALSE, # OK
                  RMtrafo(new=code), x, grid=grid, 
                  coords.new_coordunits=global.units,
                  coords.new_coord_system = "keep")
      
      if (length(y) != 0)         
        y <- RFfctn(COPY=FALSE, #OK
                    RMtrafo(new=code), y, grid=grid, 
                    coords.new_coordunits=global.units,
                    coords.new_coord_system = "keep")
      
      if (new_coord_system == "cartesian") {
        Raumdim <- max(3, Raumdim)
        spatialdim <- Raumdim
      }
      base::dim(x) <- c(length(x) /Raumdim, Raumdim)
                                        #x <- t(x)

      ## never try to set the following lines outside the 'if (new_coord_system'
      ## as in case of ..="keep" none of the following lines should be set
      setRFoptions(coords.coord_system = 
                     if (new_is_cartesian) "cartesian" else new_coord_system)
      grid <- FALSE
    } else if (!(new_coord_system %in% c("keep", "sphere", "earth"))) {
      warning("unknown new coordinate system")
    }
  }

  if (has.time.comp <- length(T)!=0) {
    Ttriple <- length(T) == 3;
    if (length(T) <= 2) Tgrid <- TRUE
    else {
      dT <- diff(T)
      Tgrid <- max(abs(diff(dT))) < dT[1] * RFopt$general$gridtolerance
    }
    if (is.na(RFopt$general$Ttriple)) {
      if (Ttriple && Tgrid)
        stop("ambiguous definition of 'T'. Set RFoptions(Ttriple=TRUE) or ",
             "RFoptions(Ttriple=FALSE)")
      if (!Ttriple && !Tgrid) stop("'T' does not have a valid format")
    } else if (RFopt$general$Ttriple) {
      if (!Ttriple)
        stop("'T' is not given in triple format 'c(start, step, length)'")
      Tgrid <- FALSE
    } else {
      if (!Tgrid) stop("'T' does not define a grid")
      Ttriple <- FALSE
    }
    if (Tgrid)
      T <- as.vector(seq2grid(list(T), "T", Tgrid,
                              RFopt$messages$warn_ambiguous,
                              RFopt$general$gridtolerance))
    totpts <- totpts * T[XLENGTH + 1]
  }

  if (!missing(dim) && !is.null(dim) && spatialdim != dim) {
    stop("'dim' should be given only when 'distances' are given. Here, 'dim' contradicts the given coordinates.")
  }

  if (ex.red) {
    if (grid) {
      x[XLENGTH + 1, ] <- pmin(x[XLENGTH + 1, ], ex.red)
      if (length(y) > 0) y[XLENGTH + 1, ] <- pmin(y[XLENGTH + 1, ], ex.red)
      totpts <- as.integer(prod(x[XLENGTH + 1, ]))
    } else {
      len <- totpts <- as.integer(min(nrow(x), ex.red^spatialdim))
      x <- x[1:len, , drop=FALSE]
      if (length(y) > 0) y <- y[1:len, , drop=FALSE]
    }
    
    if (has.time.comp) {
      T[XLENGTH + 1] <- min(T[XLENGTH + 1], 3)
      totpts <- as.integer(totpts * T[XLENGTH + 1])
    }
  }
  
  
  ## keep exactly the sequence up to 'grid'
  if (length(x) > 0) {
    if (storage.mode(x) != "double") storage.mode(x) <- "double"
  } else x <- double(0)
  if (length(y) > 0) {
    if (storage.mode(y) != "double") storage.mode(y) <- "double"
  } else y <- double(0)


  ## Print("UNIFY-XT ", totpts);

  L <- list(x=x, #0
            T=as.double(T), #1
            grid=as.logical(grid), #2
            spatialdim=as.integer(spatialdim), #3
            has.time.comp=has.time.comp, #4
            dist.given=FALSE, #5
            y=y, #6
            Ty=as.double(T), #7
            gridy=as.logical(grid), #8
            totpts=as.integer(totpts), ## 9, nr of locations
       #     l=as.integer(len),           ## 10, physical "length/cols" of input
            totptsy=if (length(y) == 0) 0L
                    else if (is.matrix(y)) nrow(y)
                    else length(y), ## 11, nr of locations
       #     ly=as.integer(len),          ## 12, physical "length/cols" of input
            coordunits = rep(curunits, spatialdim + has.time.comp),  #13
            new_coordunits =
              if (change.of.units) newunits
              else rep(newunits, length = spatialdim + has.time.comp) #14
            )
  class(L) <- "UnifyXT"

  return(L)  
}


splitUnifyXT <- function(x, split) {
  if (!is(x,"UnifyXT") || length(x) != 1 || !is.list(split) || x[[1]]$grid ||
      x[[1]]$has.time.comp || x$dist.given || length(y) > 0)
    stop(CONTACT)
  y <- x[[1]]
  y$x <- NULL
  ans <- lapply(split, function(s) {
    y$x <- x[[1]]$x[s, , drop=FALSE]
    y$totpts <- length(s)
    y
  })
  class(ans) <- "UnifyXT"
  ans
}

trafo.to.C_UnifyXT <- function(new) {
  if (is(new, "C_UnifyXT")) return(new)
  if (!is(new, "UnifyXT")) stop("wrong argument to C_UnifyXT")
  if (is.list(new[[1]])) {
    for(i in 1:length(new)) {
      if (length(new[[i]]$x)>0 && !new[[i]]$grid) new[[i]]$x = t(new[[i]]$x)
      if (length(new[[i]]$y)>0 && !new[[i]]$grid) new[[i]]$y = t(new[[i]]$y)
    }
  } else {
    if (length(new$x)>0 && !new$grid) new$x = t(new$x)
    if (length(new$y)>0 && !new$grid) new$y = t(new$y)
  }
  class(new) <- "C_UnifyXT"
  new
}


C_UnifyXT <- function(...) 
  return(if (...length() == 1 && is(...elt(1), "C_UnifyXT")) ...elt(1)
         else if (...length() == 1 && is(...elt(1), "UnifyXT"))
           trafo.to.C_UnifyXT(...)
         else trafo.to.C_UnifyXT(UnifyXT(...)))


## used by RFratiotest, fitgauss, Crossvalidation, likelihood-ratio,  RFempir
UnifyData <- function(model, x, y=NULL, z=NULL, T=NULL,  grid=NULL, data,
                      distances=NULL, RFopt,
                      dim=NULL, allowFirstCols=TRUE, vdim = NULL,
                      params=NULL,
                      further.models=NULL, ## != NULL iff RFfit
                      model.dependent.further.models = FALSE,
                      return_transform = FALSE,
                      ...) {

  ##  if (missing(x)) Print(data, T) else Print(data, T, x)

  ##  if (missing(x)) Print(data, T) else Print(data, T, x)
  ##  if (!missing(model)) print(model)
  ##if (!missing(further.models))print(further.models)

  
  RFopt <- internalRFoptions(COPY=FALSE, internal.examples_reduced=FALSE) # OK
  PL <- as.integer(RFopt$basic$printlevel)
  .Call(C_setlocalRFutils, NULL, min(PL, PL_IMPORTANT)) # fuer UnifyXT-Aufrufe
  duplicated <- which(RFopt$general$duplicated_loc == DUPLICATEDLOC_NAMES) - 1
#  Print(duplicated , RFopt$general$duplicated_loc, DUPLICATEDLOC_NAMES)
  
  add.na <- RFopt$fit$addNAlintrend ## * !is.null(further.models) 22.10.19 ## for fitting

  ## Achtung! Reicht nicht, oben auf NULL zu setzen 9.11.20!!
  if (missing(dim)) dim <- NULL
  if (missing(grid)) grid <- NULL
  if (missing(further.models))  further.models <- NULL
  orig.further.models <- further.models
  
  dist.given <- !missing(distances) && length(distances)>0
  ##  if (dist.given) {   printf("geht nicht");  bug; }

  prepRecall <- matrix.indep.of.x.assumed <- FALSE  
  neu <- RFsp.info <- C_coords <- NULL
  if (missing(data)) stop("missing data")
  missing.x <- missing(x) || is.null(x)
  missing.vdim <- missing(vdim) || length(vdim) == 0 || vdim == 0
  missing.further <- length(further.models) == 0
  
  ##  Print(is(data[[1]], "RFsp"), is(data, "RFsp"));
  ##  str(data)

  if (isSpObj(data)) data <- sp2RF(data)
  if (!is.list(data) || is.data.frame(data)) data <- list(data)
  orig.data <- data
  xdim <- 0

  if (isRFsp <- is(data[[1]], "RFsp")) {
    if ( (!missing.x && length(x)!=0) || length(y)!=0   || length(z) != 0 ||
         length(T) !=  0 || dist.given || length(dim)!=0 || length(grid) != 0)
      stop("data object already contains information about the locations. So, none of 'x' 'y', 'z', 'T', 'distance', 'dim', 'grid' should be given.")
    sets <- length(data)
    x <- RFsp.info <- vector("list", sets)
    
    if (!is.null(data[[1]]@.RFparams)) {
      if (length(vdim) > 0) stopifnot( vdim == data[[1]]@.RFparams$vdim)
      else vdim <- data[[1]]@.RFparams$vdim
    }
    
    repet <- numeric(length(data))
    for (i in 1:length(data)) {
      xi <- list()
      xi$grid <- isGridded(data[[i]])
      compareGridBooleans(grid, xi$grid)
      
      ## PUT data.columns here, if allow for different
      ## column labels in different sets 
      
      dimensions <-
        if (xi$grid) data[[i]]@grid@cells.dim else nrow(data[[i]]@data)
      totpts <- prod(dimensions)
      if (data[[i]]@.RFparams$vdim > 1) {
        dimensions <- c(dimensions, data[[i]]@.RFparams$vdim)      
        if (RFopt$general$vdim_close_together)
          dimensions <- dimensions[c(length(dimensions),
                                     1:(length(dimensions)-1))]
      }
      L <- nrow(data[[i]]@data) * ncol(data[[i]]@data)
      repet[i] <- L[i] / prod(dimensions)
      if (repet[i] != as.integer(repet[i]))
        stop("number of calculated repetitions does not match",
             " the length of the data.")
      if (repet[i] > 1) dimensions <- c(dimensions, repet[i])
      
      RFsp.info[[i]] <- list(data.params = data[[i]]@.RFparams,
                             dimensions = dimensions,
			     gridTopology=if (xi$grid) data[[i]]@grid else NULL,
			     coords = if (!xi$grid) data[[i]]@coords else NULL)
      
      tmp <- RFspDataFrame2conventional(data[[i]])
      xi$x <- tmp$x
      if (!is.null(tmp$T)) xi$T <- tmp$T
      data[[i]] <- as.matrix(tmp$data)
      base::dim(data[[i]]) <- c(totpts, L / totpts)
      x[[i]] <- xi
    }
    
    neu <- UnifyXT(x=x)

#    Print("A")
    
    info <- data.columns(data=orig.data, model=model, ## ultralangsam !!
                         x = neu,  RFopt=RFopt, vdim=vdim,
                         params=params, add.na=add.na,
                         return_transform = return_transform,
                         ...)
    if (!is.null(info$is.var)) {	
      sel <- info$is.var
      for (i in 1:length(data)) {
        nr <- nrow(data[[i]])
        r <- repet[i]
       # base::dim(data[[i]]) <- c(nr, ncol(data[[i]]) / r, r)
        if (!is.null(info$data.trafo))
          stop("currently formulae on the left hand side are not supported")
       	data[[i]] <- data[[i]][, sel, drop=FALSE]
        base::dim(data[[i]]) <-  c(nr, length(sel))
	if (!is.null(names(sel))) colnames(data[[i]]) <- names(sel)
      }
    }
    
  } else { # !isRFsp   
    sets <- length(data)
    for (i in 1:sets)
      if (is.vector(data[[i]])) data[[i]] <- as.matrix(data[[i]])
      else if(!is.numeric(data[[i]]) && !is.data.frame(data[[i]]))
        stop("'data' not numeric")

    if (dist.given) {  
      stopifnot(missing.x || length(x)==0, length(y)==0, length(z)==0)
      if (!is.list(distances)) distances <- list(distances)
      if (length(distances) != sets)
        stop("number of sets of distances does not match number of sets of data")
      for (i in 1:sets) 
        if (any(is.na(data[[i]])))
          stop("missing data are not allowed if distances are used.")

      stopifnot(missing(T) || length(T)==0)
      if (is.matrix(distances[[1]])) {
         spatialdim <- tsdim <- xdimOZ <- nrow(distances[[1]])
        if (length(dim) > 0 && dim != spatialdim)
          stop("unclear specification of the distances: either the distances is given as a vector or distance vectors should given, where the number of rows matches the spatial dimension")
      } else {
        xdimOZ <- 1L
        spatialdim <- tsdim <- as.integer(dim)      
      }
      
      coordunits <- RFopt$coords$coordunits
      has.time.comp <- FALSE      

      neu <- UnifyXT(distances=distances, dim = spatialdim)
      xdim <- spatialdim
#      Print("B")
     info <- data.columns(data, ## auch PM2
                          model=model,  ## auch PM2
                          xdim=xdim, ## auch PM2
                          x = neu, ## auch PM2
                          RFopt = RFopt, ## nur data.col
                          force=allowFirstCols, ## nur datalcol
                          halt=!allowFirstCols, ## nur datalcol
                           vdim=vdim, ## nur datacol
                           params=params, ## nur PM2
                           add.na=add.na, ## nur PM2
                          return_transform = return_transform,
                             ...) ## nur PM2
    } else { ## distances not given      
      ##      Print("distances not given")
      if (!missing.x) {
        if (is.data.frame(x)) x <- data.matrix(x)
        if (is.list(x)) {
          if (length(y)!=0 || length(z)!=0 || length(T)!=0)
            stop("if x is a list, then 'y', 'z', 'T' may not be given")       
          
          if (!is.list(x[[1]])) {
            if (length(data) == 1) x <- list(x)
            else stop("number of sets of 'x' and 'data' differ")
          }
          
        } else {
          x <- list(x=x)
          if (length(y)!=0) {
            stopifnot(!is.list(y))
            x$y <- y
          }
          if (length(z)!=0) {
            stopifnot(!is.list(z))
            x$z <- z
          }
          if (length(T)!=0) {
            stopifnot(!is.list(T))
            x$T <- T
          }
          if (!missing(grid) && !is.null(grid))
            x$grid <- grid
          ##if (!is.list(data)) data <- list(data.matrix(data))
          x <- list(x)          
        }
        neu <- UnifyXT(x=x)
      } else neu <- NULL
      
      ##
     if (missing.x) xdim <- dim      
      info <- data.columns(data, model=model,
                           xdim=xdim,
                           x = neu,
                           RFopt = RFopt,
                           only.data = !missing.x,  ## nur datacol
                           force=allowFirstCols, halt=!allowFirstCols,
                           vdim=vdim, params=params, add.na=add.na,
                        return_transform = return_transform,
                        ...)
##      Print(info)
    } # ! distance
    sets <- length(data)
    
  } # !isRFsp

##  Print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

  newmodel <- info$model$model
  
  varnames <- names(info$is.var)
  varnames[varnames == VOID] <- ""
  coordnames <- names(info$is.x)
  if (!is.null(info$model$C_coords)) C_coords <- info$model$C_coords

  
  ## Print(isRFsp, dist.given, missing.x)

  if (!isRFsp && !dist.given) { ## dritter, benoetigt def von varnames  
    if (missing.x) { ## dec 2012: matrix.indep.of.x.assumed        
      x <- list()
      if (length(info$is.x) == 0) {
        matrix.indep.of.x.assumed <- TRUE
        Note("detection", "Columns representing coordinates not found. So no spatial context is assumed.")
        for (i in 1:sets) {
          x[[i]] <- 1:nrow(data[[i]])
          storage.mode(x[[i]]) <- "numeric"
        }
      } else {
          for (i in 1:sets) {
            xx <- data.matrix(data[[i]][, info$is.x, drop=FALSE])#if data.frame
            colnames(xx) <- coordnames
            x[[i]] <- list(x=xx, grid=FALSE)
          }
      } ## length(info$is.x) != 0
      
      for (i in 1:sets) {
        data[[i]] <- data[[i]][ , info$is.var, drop=FALSE]
##Print(i, data[[i]], info$is.var)        
        if (!is.null(info$data.trafo)) {
          stop("currently formulae on the left hand side are not supported")
        }
        colnames(data[[i]]) <- varnames
        data[[i]] <- data.matrix(data[[i]])
      }
      
      neu <- UnifyXT(x=x) 
      C_coords <- trafo.to.C_UnifyXT(neu)
      
    } ## missing xgiven; KEIN ELSE, auch wenn nachfolgend z.T. gedoppelt wird
  }
  
  for (i in 1:sets) 
    if (is.data.frame(data[[i]])) data[[i]] <- data.matrix(data[[i]])


  attr(data, "Unified") <- TRUE
 
  if (!dist.given) { ##   x coordinates, not distances
 
    if (!is.list(neu[[1]])) neu <- list(neu)

    coordunits<- neu[[1]]$coordunits
    spatialdim <- as.integer(neu[[1]]$spatialdim)
    has.time.comp <- neu[[1]]$has.time.comp
    tsdim <- as.integer(spatialdim + has.time.comp)
    
    if (has.time.comp &&
        any(sapply(neu, function(x) x$T[2]) <= RFopt$nugget$tol))
      stop("the step of the temporal grid component ",             
           if (any(sapply(neu, function(x) x$T[2]) == 0))
             "equals zero." else "is smaller than nugget tolerance 'tol'.")
    
    if (any(sapply(neu, function(x) x$grid && any(x$x[2,]<=RFopt$nugget$tol))))
          stop("the step of some spatial grid component ",             
               if (any(sapply(neu, function(x) x$grid && any(x$x[2,] == 0))))
                 "equals zero." else "is smaller than nugget tolerance 'tol'.")

    
    if (duplicated == DUPLICATEDLOC_SCATTER){
      for (i in 1:length(neu))
        if (!neu[[i]]$grid) .Call(C_scatter, neu[[i]]$x)
        else if (any(neu[[i]]$x[XSTEP +1,] == 0))
          stop("step=0 detected and grid points cannot be scattered.")
    }
   
    xdimOZ <- ncol(neu[[1]]$x)
  }

  if (is.null(coordnames))
    coordnames <- SystemCoordNames(locinfo=neu[[1]], opt=RFopt$coords) 
  
  totpts <- sapply(neu, function(x) x$totpts)
  ldata <- sapply(data, length)

  if (missing(model)) {
    if (!missing.further)
      stop("'model' is not given, but 'further.models'")
   } else {
     model.vdim <- info$model$vdim
     if (model.vdim != 0) {
       if (length(vdim) == 0) vdim <- model.vdim
       else if (vdim != model.vdim)
         stop("given multivariate dimension differs from the dimensions expected by the model")
     }
     
     if (is.null(newmodel)) stop("new model not found.", CONTACT)
     
     if (!missing.further) {
       #Print(coordnames, varnames,extractRepeatedNames(varnames))
       setRFoptions(coords.coordnames = coordnames,
                    coords.varnames = extractRepeatedNames(varnames))
       for (m in 1:length(further.models))
         if (!is.null(further.models[[m]]) &&
             !is.numeric(further.models[[m]])) {           
           further.models[[m]] <-
             if (model.dependent.further.models) {

               str(model)
               
               nf <- names(further.models)
               PrepareModel2(model=model, params=params, xdim=xdim, x=neu,
                             add.na=add.na, data=orig.data, arg = nf[m],
                             further.model=further.models[[m]], ...)$model
             } else {
               PrepareModel2(further.models[[m]], params=params,
                             xdim = xdim, x=neu, add.na=add.na, 
                             data=orig.data, ...)$model
             }
          }
     }
   }

  if (length(vdim) == 0) {
    repetitions <- sapply(data, ncol)
    if (min(repetitions) == 1) vdim <- 1 ## ggT==1 wuerde reichen
  } else repetitions <- as.integer(ldata / (totpts * vdim))

  if (length(vdim) == 0 && !missing(model) && missing.x) {
    ## Verbesserung kann nur erreicht werden, wenn vorher
    ##            keine x-Koord fuer Modell bekannt waren.
    return(UnifyData(model=model, x=neu, data=data,
                     distances=distances, RFopt=RFopt,
                     dim=dim, allowFirstCols=allowFirstCols,
                     vdim = vdim,
                     params=params,
                     return_transform = return_transform,
                     further.models=orig.further.models, ...))
  }

  if (length(vdim) == 0) {
    if (length(info$is.x) + length(info$is.var) == ncol(orig.data[[1]])) {
      info$data.info <- "guess"
      vdim <- 1
    } else stop("Not detectable which columns are data columns.\nPlease set RFoptions(varnames=) or RFoptions(varidx=) accordingly.")
  }
  
  if (PL > PL_SUBIMPORTANT ||
      (PL >= PL_IMPORTANT && info$data.info != "safe")
      && .Call(C_areDataIdxDifferentFromFormer, info$is.var)
      ) {
    nn <- min(length(info$is.var), RFopt$messages$vec_len)
    reduced <- nn < length(info$is.var)
    nnRep <- min(length(repetitions), RFopt$messages$vec_len)
    reducedRep <- nn < length(repetitions)

    Note("detection", "Using columns ", paste(info$is.var[1:nn], collapse=","),
         if (reduced) "...",
         " as data columns",
         if (any(repetitions) > 1)
           paste("with", paste(repetitions[1:nnRep], collapse=",") ,
                  if (reducedRep) "...",
                 "repetitions"),
           
         ".")
  }
 
  if (any(repetitions)==0) stop("no or not sufficiently many data are given")
  if (!is.na(vdim) && any(ldata != repetitions * totpts * vdim))
    stop("mismatch of data dimensions")

  vrep <- repetitions * vdim

  for (i in 1:length(data)) {
    base::dim(data[[i]]) <- c(ldata[i] / vrep[i], vrep[i])
  }


  ##  Print(vdim, vrep, repetitions); stopifnot(vdim ==2)
                                        #  Print(repetitions)

  ##  Print("end unify") 

  .Call(C_setlocalRFutils, NULL, PL) # fuer UnifyXT-Aufrufe
  return(list(
    ## coords = expandiertes neu # #
    model = if (missing(model)) NULL else newmodel,
    orig.model = if (missing(model)) NULL else model,
    further.models = further.models,
    transform = info$model$transform,
 ##   dimdata=dimdata, isRFsp = isRFsp,  # 3.2.0 : not given anymore
## oben auch noch loeschen!!
 ##  len = len,    # 3.2.0 : not given anymore
    data=data,
    RFsp.info = RFsp.info,
    coords = neu,    
    dist.given=dist.given,
    spatialdim=spatialdim,
    tsdim=tsdim,
    coordunits=coordunits,
    has.time.comp = has.time.comp,
    matrix.indep.of.x.assumed = matrix.indep.of.x.assumed,
    xdimOZ = xdimOZ,
    vdim = vdim,
    coordnames=coordnames,
    varnames=varnames,
    is.var.orig = info$is.var,
    repetitions = repetitions,
    C_coords = C_coords,
    orig.trend = list(...)$trend,
    DataNames = info$M$data.coordnames ## info$M might be NULL already;
    ##                                       from model only!
  ))
}

