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


seq2grid <- function(x, name, grid, warn_ambiguous, gridtolerance) {
  xx <- matrix(nrow=3, ncol=length(x))
  step0 <- rep(FALSE, length(x))
  gridnotgiven <- missing(grid) || length(grid) == 0
  
  for (i in 1:length(x)) {
    if (length(x[[i]]) == 1) {
      xx[,i] <- c(x[[i]], 0, 1)
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

  if (FALSE && gridnotgiven && warn_ambiguous && length(x) > 1) {
    RFoptions(internal.warn_ambiguous = FALSE)
    message("Ambiguous interpretation of coordinates. Better give 'grid=TRUE' explicitly. (This message appears only once per session.)")
  }

  if (any(step0)) {
    if (all(step0)) {
      if (gridnotgiven) return(FALSE)
      else stop("Within a grid, the coordinates must be distinguishable")
    } else {
      if (gridnotgiven && warn_ambiguous) {
        RFoptions(internal.warn_ambiguous = FALSE)
        warning("Interpretation as degenerated grid. Better give 'grid' explicitely. (This warning appears only once per session.)")
      }
    }
  }

  return(xx)
}

RFearth2cartesian <- function(coords, units=NULL, system = "cartesian",
                              grid=FALSE) {
  if (is.character(system)) system <- pmatch(system, ISO_NAMES) - 1
  stopifnot(system %in%
            c(CARTESIAN_COORD, GNOMONIC_PROJ, ORTHOGRAPHIC_PROJ))
  if (missing(units) || is.null(units)) {
    global.units <- RFoptions()$coords$new_coordunits[1]
    units <-if (global.units[1] == "") "km" else global.units
  }
  if (!is.matrix(coords)) coords <- t(coords)
  res <- RFfctn(RMtrafo(new=system), coords, grid=grid,
                coords.new_coord_system = "keep",
                coords.new_coordunits=units,
                coords.coord_system="earth")
  dimnames(res) <- list(NULL, c("X", "Y", "Z", "T")[1:ncol(res)])
  return(res)
}

RFearth2dist <- function(coords, units=NULL, system="cartesian",
                         grid=FALSE, ...) {
  if (is.character(system)) system <- pmatch(system, ISO_NAMES) - 1
  stopifnot(system %in%
            c(CARTESIAN_COORD, GNOMONIC_PROJ, ORTHOGRAPHIC_PROJ))
  if (missing(units) || is.null(units)) {
    global.units <- RFoptions()$coords$new_coordunits[1]
    units <- if (global.units[1] == "") "km" else global.units 
  }
  if (!is.matrix(coords)) coords <- t(coords)
  z <- RFfctn(RMtrafo(new=system), coords, grid=grid,
                coords.new_coord_system = "keep",
                coords.new_coordunits=units,
                coords.coord_system="earth")
  return(dist(z, ...))
}

UnifyXT <- function(x=NULL, y=NULL, z=NULL, T=NULL, grid, distances=NULL,
                    dim=NULL, # == spatialdim!
                    length.data,
                    y.ok = FALSE, 
                    printlevel = RFopt$basic$printlevel,
                    allow_duplicated =
                      RFopt$internal$allow_duplicated_locations ){
  ## do not pass anything on "..." ! --- only used for internal calls
  ## when lists are re-passed

  ## converts the given coordinates into standard formats
  ## (one for arbitrarily given locations and one for grid points)
  #print("UnifyXT in convert.R")#Berreth
  
  RFopt <- RFoptions()
##Print("Enter")
  
  if (!missing(x) && !is.null(x) ) {
    if (is(x, "UnifyXT")) return(x)  
    if (is.list(x) && !is.data.frame(x)) {
      if (!is.list(x[[1]]))
        return(do.call("UnifyXT", c(x, list(printlevel=printlevel))))
      L <- list()
      for (i in 1:length(x)) {        
        L[[i]] <-
          if (is(x[[i]], "UnifyXT")) x[[i]]
          else do.call("UnifyXT", c(x[[i]], list(printlevel=printlevel)))
      }
      if (length(x) > 1) {    
        if (!all(diff(sapply(L, function(x) x$has.time.comp)) == 0) ||
          !all(diff(sapply(L, function(x) x$spatialdim)) == 0))
          stop("all sets must have the same dimension")
        if (!all(diff(sapply(L, function(x) x$dist.given)) == 0))
        stop("either all the sets must be based on distances or none")
      }
      class(L) <- "UnifyXT"
      return(L)
    }
  }

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
        L[[i]] <- do.call("UnifyXT", list(distances=distances[[i]], dim=dim,
                                          printlevel=printlevel))
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
        if (any(rowSums(distances != 0) == 0))
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
   
    if (ex.red && len > ex.red^2 / 2) {
      LEN <- as.integer(ex.red)
      len <- as.integer(LEN * (LEN - 1) / 2)
      x <- if (is.matrix(x)) x[1:len ,] else x[1:len]
    } else {
      LEN <- as.integer(1e-9 + 0.5 * (1 + sqrt(1 + 8 * len)))
      if (LEN * (LEN-1) / 2 != len) LEN <- NaN
    }

    ## keep exactly the sequence up to 'distances'
    if (storage.mode(x) != "double") storage.mode(x) <- "double"
    L <- list(x = as.matrix(x), #0
              y = double(0),   #1
              T= double(0),  #2
              grid = FALSE, #3
              spatialdim=as.integer(dim),#4
              has.time.comp=FALSE, #5
              dist.given = TRUE, #6
              restotal = LEN, ## number of points
              l = LEN, ## ?? physical length??
              coordunits = rep(curunits, length=dim),
              new_coordunits = rep(newunits, length=dim)
              )
    class(L) <- "UnifyXT"
    return(L)
  }

   
 
 #stopifnot(!missing(x))
  if (is(x, "RFsp") || isSpObj(x)) {
    return(UnifyXT(x=coordinates(x), y=y, z=z, T=T, grid=grid,
                   distances=distances, dim=dim, length.data=length.data,
                   y.ok=y.ok, printlevel=printlevel))
  }    
  if (is.raster(x)) x <- as(x, 'GridTopology')
 
  if ((missing(grid) || length(grid) == 0) && !missing(length.data)) {
    new <-  Try(UnifyXT(x=x, y=y, z=z, T=T, grid=TRUE, distances=distances,
                        dim=if (!missing(dim)) dim,
                        length.data = length.data, y.ok =y.ok,
                        printlevel = printlevel), silent=TRUE)
    if (grid <- !is(new, "try-error")) {
      ratio <- length.data / new$restotal

      if (grid <- ratio == as.integer(ratio)) {
        if (printlevel>=PL_IMPORTANT && new$spatialdim > 1)
          message("Grid detected. If it is not a grid, set grid=FALSE. (To avoid this message set 'grid'.)\n")
      }
    }
    return(if (grid) new else {
      UnifyXT(x, y, z, T, grid=FALSE, distances,
              if (length(distances) > 0) dim=1,
              length.data = length.data,
              printlevel = printlevel) }
           )
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
 
  if (is.matrix(x)) {
    if (!is.numeric(x)) stop("x is not numeric.")
    if ((d <- ncol(x)) > 1 &&
        x[1, ncol(x)] == x[2, ncol(x)]) ## fast check whether expanded grid
    { ## currentlyl, for y not programmed yet
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

    
    if (length(z)!=0) stop("If x is a matrix, then z may not be given")
    if (length(y)!=0) {
      if (!y.ok) stop("If x is a matrix, then y may not be given")
      if (length(T)!=0)
        stop("If x is a matrix and y is given, then T may not be given")
      if (!is.matrix(y) || ncol(y) != ncol(x) ||
          nrow(x)==3 && nrow(y)!=3 && ((missing(grid) || length(grid) == 0) ||
                                grid))
        stop("y does not match x (it must be a matrix)")
    }

    change.of.units <- coord_system == COORD_SYS_NAMES[coord_auto + 1] &&
      ncol(x) >= 2 && ncol(x) <= 3 && !is.null(n <- dimnames(x)[[2]])
    if (change.of.units) {
      if (any(idx <- earth_coordinate_names(n))) {
        if (length(idx) == 2 && !all(idx == 1:2))
          stop("earth coordinates not in order longitude/latitude")
        cur <- curunits[1]
        newunits <- RFopt$coords$new_coordunits
        curunits <- RFopt$coords$coordunits
        curunits[1:2] <- COORD_NAMES_EARTH[1:2]
        if (newunits[1] == "") newunits[1] <-  UNITS_NAMES[units_km + 1]
        newunits[2:3] <- newunits[1]                
        if (RFopt$internal$warn_coordinates) {
          message("\n\nNOTE: ",
                  "Earth coordinates are detected.\n",
                   "If this is not what you wish, change option 'coord_system'",
                  "\nto RFoptions(coord_system = \"cartesian\").\n",
                  "Note further that angles in R.cos, R.sin, R.tan, RMangle",
                  " are now expected\nin DEGREE and ",
                  "R.acos, R.asin, R.atan, R.atan2 return results in degree.\n",
                 "\nCurrent units are ",
                  if (cur=="") "not given." else paste0("'", cur, "'."),
                  "\nIn rare cases earth coordinates will be transformed ",
                  "within submodels.",
                  "\nThen it will be transformed into units of '",
                  newunits[1],
                  "'. In particular, the\nvalues of all scale parameters of ",
                  "any of these submodels, defined in R^3,\nare ",
                  "understood in units of '", newunits[1],
                  "'. Change option 'units' if ",
                  "necessary.\n" ,
                  "(This message appears in long form only once per session.)\n"
                  )
        } else message("Earth coordinates detected. In case it is necessary they will be transformed ",
                     " into units of ",  newunits[1], ".")
        coord_system <- COORD_SYS_NAMES[earth + 1]
        RFoptions(coords.coord_system = coord_system,
                  coords.coordunits = curunits,
                  coords.new_coordunits = newunits,
                  internal.warn_coordinates=FALSE)
       } else {
         RFoptions(coords.coord_system =  COORD_SYS_NAMES[cartesian + 1])
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
    } # else {

    if ((missing(grid) || length(grid) == 0) &&
        any(apply(x, 2, function(z) (length(z) <= 2) || max(abs(diff(diff(z))))
                  > RFopt$general$gridtolerance))) {
      grid <- FALSE
    }

    if ((missing(grid) || length(grid) == 0) || !is.logical(grid)) {
      grid <- TRUE
      if (spatialdim > 1 && RFopt$internal$warn_ambiguous) {
        RFoptions(internal.warn_ambiguous = FALSE)
        warning("Ambiguous interpretation of the coordinates. Better give the logical parameter 'grid=TRUE' explicitely. (This warning appears only once per session.)")
      }
    }

    if (grid && !is.GridTopology) {
      if (gridtriple <- len==3) {
        if (printlevel >= PL_SUBIMPORTANT && RFopt$internal$warn_oldstyle) {
         message("Grid coordinates recognized of size ",
                  paste0("seq(", x[1,], ", by=", x[2,], ", len=", x[3,], ")",
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
  } else { ## x, y, z given separately
    if (length(y)==0 && length(z)!=0) stop("y is not given, but z")
    xyzT <- list(x=if (!missing(x) && !is.null(x)) x, y=y, z=z, T=T)
    for (i in 1:4) {
      if (!is.null(xyzT[[i]]) && !is.numeric(xyzT[[i]])) {
        if (printlevel>PL_IMPORTANT) 
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
      else if (xor(newgrid, grid) && RFopt$internal$warn_on_grid) {
        RFoptions(internal.warn_on_grid = FALSE)
        message("coordinates", if (grid) " do not",
                " seem to be on a grid, but grid = ", grid)
      }
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
    if (gridtriple) {
      if (len != 3)
        stop("In case of simulating a grid with option gridtriple, exactly 3 numbers are needed for each direction")
      lr <- x[3,] # apply(x, 2, function(r) length(seq(r[1], r[2], r[3])))
      ##x[2,] <- x[1,] + (lr - 0.999) * x[3,] ## since own algorithm recalculates
      ##                               the sequence, this makes sure that
      ##                               I will certainly get the result of seq
      ##                               altough numerical errors may occurs
      restotal <- prod(x[3, ])
      if (length(y)!=0 && !all(y[3,] == x[3,]))
        stop("the grids of x and y do not match ")        
    } else {     
      xx <- seq2grid(x, "x",  grid,
                     RFopt$internal$warn_ambiguous, RFopt$general$gridtolerance)
     if (length(y)!=0) {
        yy <- seq2grid(y, "y", grid,
                       RFopt$internal$warn_ambiguous,
                       RFopt$general$gridtolerance)
        if (xor(is.logical(xx), is.logical(yy)) ||
            (!is.logical(xx) && !all(yy[3,] == xx[3,])))
          stop("the grids for x and y do not match")      
      }
      if (missing(grid) || length(grid) == 0) grid <- !is.logical(xx)       
      if (grid) {
        x <- xx
        if (length(y) != 0) y <- yy
        restotal <- prod(len)
        len <- 3
      } else {
        x <- sapply(x, function(z) z)
        if (length(y) != 0) y <- sapply(y, function(z) z)
      }
    }
    if (grid && any(x[3, ] <= 0))
      stop(paste("step must be postive. Got as steps",
                 paste(x[3,], collapse=",")))
    ##if (len == 1) stop("Use grid=FALSE if only a single point is simulated")
  }
 
  if (!grid) {
    restotal <- nrow(x)
    if (length(y)==0 && !allow_duplicated) {
      if (restotal < 200 && any(as.double(dist(x)) == 0)) {
        d <- as.matrix(dist(x))
        diag(d) <- 1
        idx <-  which(as.matrix(d) ==0)
        if (printlevel>PL_ERRORS)
          Print(x, dim(d), idx , cbind( 1 + ((idx-1)%% nrow(d)), #
                                       1 + as.integer((idx - 1)  / nrow(d))) ) 
        stop("Some locations are not distinguishable. If this is intentional, see ?RPnugget and ?RFoptions 'allow_duplicated_locations' to deal with repeated measurements.")
      }
      ## fuer hoehere Werte con total ist ueberpruefung nicht mehr praktikabel
    }
  }
 
  if (coord_system == "earth") {
    # if (ncol(x) > 4) stop("earth coordinates have maximal 3 components")
    opt <- RFoptions()$coords ## muss nochmals neu sein
    global.units <- opt$new_coordunits[1]
    if (global.units[1] == "") global.units <- "km"
   
    Raumdim <- ncol(x) #if (grid) ncol(x) else
    new_is_cartesian <- new_coord_system %in% CARTESIAN_SYS_NAMES
    if (new_is_cartesian) {
      if (sum(idx <- is.na(opt$zenit))) {
        zenit <- (if (grid) x[1, 1:2] + x[2, 1:2] * (x[3, 1:2] - 1)
                  else if (opt$zenit[!idx] == 1) colMeans(x[, 1:2])
                  else if (opt$zenit[!idx] == Inf) colMeans(apply(x[, 1:2], 2,
                                                                  range))
                  else stop("unknown value of zenit"))
         RFoptions(zenit = zenit)
      }

      code <- switch(new_coord_system,
                     "cartesian" = CARTESIAN_COORD,
                     "gnomonic" = GNOMONIC_PROJ,
                     "orthographic" = ORTHOGRAPHIC_PROJ,
                     stop("unknown projection method")
                     )
      message("New coordinate system: ", new_coord_system, ".\n")
      x <- RFfctn(RMtrafo(new=code), x, grid=grid, 
                   coords.new_coordunits=global.units,
                   coords.new_coord_system = "keep")
      
       if (length(y) != 0)         
         y <- RFfctn(RMtrafo(new=code), y, grid=grid, 
                   coords.new_coordunits=global.units,
                   coords.new_coord_system = "keep")
     
      if (new_coord_system == "cartesian") {
        Raumdim <- max(3, Raumdim)
        spatialdim <- Raumdim
      }
      dim(x) <- c(length(x) /Raumdim, Raumdim)
      #x <- t(x)

      ## never try to set the following lines outside the 'if (new_coord_system'
      ## as in case of ..="keep" none of the following lines should be set
      RFoptions(coords.coord_system = 
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
                              RFopt$internal$warn_ambiguous,
                              RFopt$general$gridtolerance))
    restotal <- restotal * T[3]
  }

   if (!missing(dim) && !is.null(dim) && spatialdim != dim) {
    stop("'dim' should be given only when 'distances' are given. Here, 'dim' contradicts the given coordinates.")
  }

  if (ex.red) {
    if (grid) {
      x[3, ] <- pmin(x[3, ], ex.red)
      if (length(y) > 0) y[3, ] <- pmin(y[3, ], ex.red)
      restotal <- as.integer(prod(x[3, ]))
    } else {
      len <- restotal <- as.integer(min(nrow(x), ex.red^spatialdim))
      x <- x[1:len, , drop=FALSE]
      if (length(y) > 0) y <- y[1:len, , drop=FALSE]
    }
    
    if (has.time.comp) {
      T[3] <- min(T[3], 3)
      restotal <- as.integer(restotal * T[3])
    }
  }
 
   
  ## keep exactly the sequence up to 'grid'
  if (length(x) > 0) {
    if (storage.mode(x) != "double") storage.mode(x) <- "double"
  } else x <- double(0)
  if (length(y) > 0) {
    if (storage.mode(y) != "double") storage.mode(y) <- "double"
  } else y <- double(0)

  L <- list(x=x, #0
            y=y, #1
            T=as.double(T), #2
            grid=as.logical(grid), #3
            spatialdim=as.integer(spatialdim), #4
            has.time.comp=has.time.comp, #5
            dist.given=FALSE, #6
            restotal=as.integer(restotal), ## 7, nr of locations
            l=as.integer(len),             ## 8, physical "length/rows" of input
            coordunits = rep(curunits, spatialdim + has.time.comp),  #9
            new_coordunits =
              if (change.of.units) newunits
              else rep(newunits, length = spatialdim + has.time.comp)) #10
  class(L) <- "UnifyXT"

 return(L)  
}


trafo.to.C_UnifyXT <- function(new) {
  if (is.list(new[[1]])) {
    for(i in 1:length(new)) {
      if (length(new[[i]]$x)>0 && !new[[i]]$grid) new[[i]]$x = t(new[[i]]$x)
      if (length(new[[i]]$y)>0 && !new[[i]]$grid) new[[i]]$y = t(new[[i]]$y)
    }
  } else {
    if (length(new$x)>0 && !new$grid) new$x = t(new$x)
    if (length(new$y)>0 && !new$grid) new$y = t(new$y)
  }
  new
}


C_UnifyXT <- function(...) return(trafo.to.C_UnifyXT(UnifyXT(...)))
   

## used by RFratiotest, fitgauss, Crossvalidation, likelihood-ratio,  RFempir
UnifyData <- function(model, x, y=NULL, z=NULL, T=NULL,  grid=NULL, data,
                      distances=NULL, RFopt, mindist_pts=2,
                      dim=NULL, allowFirstCols=TRUE, vdim = NULL,
                      params=NULL,
                      further.models=NULL, ## != NULL iff RFfit
                      ...) {

  ##  if (missing(x)) Print(data, T) else Print(data, T, x)
  ##  if (!missing(model)) print(model)
  ##if (!missing(further.models))print(further.models)

  
  RFoptOld <- internal.rfoptions(internal.examples_reduced=FALSE)
  RFopt <- RFoptOld[[2]]
  PL <- as.integer(RFopt$basic$printlevel)
  add.na <- RFopt$fit$addNAlintrend ## * !is.null(further.models) 22.10.19 ## for fitting

  ## Achtung! Reicht nicht, oben auf NULL zu setzen 9.11.20!!
  if (missing(dim)) dim <- NULL
  if (missing(grid)) grid <- NULL
  if (missing(further.models)) grid <- NULL
  orig.further.models <- further.models
  
  dist.given <- !missing(distances) && length(distances)>0
  ##  if (dist.given) {   printf("geht nicht");  bug; }

  prepRecall <- matrix.indep.of.x.assumed <- FALSE  
  rangex <- neu <- gridlist <- RFsp.info <- mindist <- C_coords <- NULL
  if (missing(data)) stop("missing data")
  missing.x <- missing(x) || is.null(x)
  missing.vdim <- missing(vdim) || length(vdim) == 0 || vdim == 0
  missing.further <- missing(further.models) || length(further.models) == 0
  
  if (isSpObj(data)) data <- sp2RF(data)
  if (!is.list(data) || is.data.frame(data)) data <- list(data)
  orig.data <- data

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

    for (i in 1:length(data)) {
      xi <- list()
      xi$grid <- isGridded(data[[i]])
      compareGridBooleans(grid, xi$grid)
      
      ## PUT data.columns here, if allow for different
      ## column labels in different sets 
          
      dimensions <-
        if (xi$grid) data[[i]]@grid@cells.dim else nrow(data[[i]]@data)
      restotal <- prod(dimensions)
      if (data[[i]]@.RFparams$vdim > 1) {
        dimensions <- c(dimensions, data[[i]]@.RFparams$vdim)      
        if (RFopt$general$vdim_close_together)
          dimensions <- dimensions[c(length(dimensions),
                                     1:(length(dimensions)-1))]
      }
      L <- nrow(data[[i]]@data) * ncol(data[[i]]@data)
      repet <- L[i] / prod(dimensions)
      if (repet != as.integer(repet))
        stop("number of calculated repetitions does not match",
             " the length of the data.")
      if (repet > 1) dimensions <- c(dimensions, repet)
 
      RFsp.info[[i]] <- list(data.params = data[[i]]@.RFparams,
                             dimensions = dimensions,
			     gridTopology=if (xi$grid) data[[i]]@grid else NULL,
			     coords = if (!xi$grid) data[[i]]@coords else NULL)
            
      tmp <- RFspDataFrame2conventional(data[[i]])
      xi$x <- tmp$x
      if (!is.null(tmp$T)) xi$T <- tmp$T
      data[[i]] <- as.matrix(tmp$data)
      dim(data[[i]]) <- c(restotal, L / restotal)
      x[[i]] <- xi
    }
 

     neu <- UnifyXT(x=x, printlevel=min(PL, PL_IMPORTANT))
    
    info <- data.columns( ## ultralangsam !!
        data=orig.data[[1]], model=model, RFopt=RFopt,
        vdim=vdim, params=params, add.na=add.na,
        x = neu,
        ...)      

    sel <- info$is.data
    if (!is.null(sel)) {	
      for (i in 1:length(data)) {
   	data[[i]] <- data[[i]][, sel, drop=FALSE]
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

        dimensions <- sapply(distances, nrow)
        spatialdim <- tsdim <- xdimOZ <- dimensions[1]
        if (length(dim) > 0 && dim != spatialdim)
          stop("unclear specification of the distances: either the distances is given as a vector or distance vectors should given, where the number of rows matches the spatial dimension")
        lcc <- sapply(distances, function(x) 0.5 * (1 + sqrt(1 + 8 * ncol(x))) )
        if (!all(diff(dimensions) == 0))
          stop("sets of distances show different dimensions")
        range_distSq <- function(M) range(apply(M, 2, function(z) sum(z^2)))
        rangex <- sqrt(range(sapply(distances, range_distSq)))
      } else {        
        xdimOZ <- 1L
        spatialdim <- tsdim <- as.integer(dim)      
        lcc <- sapply(distances, function(x) if (is.matrix(x)) -1
                                        else 0.5 * (1 + sqrt(1 + 8* length(x))))
        rangex <- range(sapply(distances, range))
      }
      mindist <- min(rangex)
      if (is.na(mindist)) mindist <- 1 ## nur 1 pkt gegeben, arbitraerer Wert
      if (mindist <= RFopt$nugget$tol) {
        if (!RFopt$eneral$allowdistanceZero)
          stop("distance with value 0 identified -- use allowdistanceZero=T?")
        mindist <- 1e-15 * (RFopt$nugget$tol == 0) + 2 * RFopt$nugget$tol

        for (i in 1:length(distances))
          if (is.vector(distances[[i]]))
            distances[[i]][distances[[i]] == 0] <- mindist
          else distances[[i]][1, apply(distances[[i]], 2,
                                       function(z) sum(z^2))] <- mindist
      }

      if (any(as.integer(lcc) != lcc))
	stop("number of distances not of form k(k-1)/2")
      coordunits <- RFopt$coords$coordunits
      has.time.comp <- FALSE      

      neu <- UnifyXT(distances=distances, dim = spatialdim) 
     info <- data.columns(data[[1]], model=model, RFopt = RFopt,
                           xdim=spatialdim,
                           force=allowFirstCols, halt=!allowFirstCols,
                           vdim=vdim, params=params, add.na=add.na,
                           x = neu,
                           ...)
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
       neu <- UnifyXT(x=x, printlevel=min(PL, PL_IMPORTANT))
      } else neu <- NULL

     ##    Print(neu)
 
     info <- data.columns(data[[1]], model=model, RFopt = RFopt,
                           xdim=if (missing.x) dim else NA,
                           force=allowFirstCols, halt=!allowFirstCols,
                           vdim=vdim, params=params, add.na=add.na,
                           x = neu,
                           ...)
    
      if (missing.x) { ## dec 2012: matrix.indep.of.x.assumed        
	x <- list()
	if (length(info$is.x) == 0) {
	  matrix.indep.of.x.assumed <- TRUE
	  if (PL > 0)
	    DetectionNote("Columns representing coordinates not found. So no spatial context is assumed.")
	  for (i in 1:sets) {
	    x[[i]] <- 1:nrow(data[[i]])
	    storage.mode(x[[i]]) <- "numeric"
	  }
	} else {
	  for (i in 1:sets) {
	    xx <- data.matrix(data[[i]][, info$is.x, drop=FALSE])#if data.frame
	    colnames(xx) <- names(info$is.x)
	    x[[i]] <- list(x=xx, grid=FALSE)
	  }
        } ## length(info$is.x) != 0
        
        for (i in 1:sets) {
          data[[i]] <- data[[i]][ , info$is.data, drop=FALSE]
 	  colnames(data[[i]]) <- names(info$is.data)
          data[[i]] <- data.matrix(data[[i]])
        }
       
        neu <- UnifyXT(x=x, printlevel=min(PL, PL_IMPORTANT))       
        C_coords <- trafo.to.C_UnifyXT(neu)

      } ## missing xgiven; KEIN ELSE, auch wenn nachfolgend z.T. gedoppelt wird
    
        ## {}
    } # ! distance
    sets <- length(data)

  } # !isRFsp 

  newmodel <- info$newmodel$model
  varnames <- names(info$is.data)
  coordnames <- names(info$is.x)
  if (!is.null(info$newmodel$C_coords)) C_coords <- info$newmodel$C_coords
  
  for (i in 1:sets) 
    if (is.data.frame(data[[i]])) data[[i]] <- data.matrix(data[[i]])


  attr(data, "Unified") <- TRUE
 
  if (!dist.given) { ##   x coordinates, not distances
 
    if (!is.list(neu[[1]])) neu <- list(neu)

    coordunits<- neu[[1]]$coordunits
    spatialdim <- as.integer(neu[[1]]$spatialdim)
    has.time.comp <- neu[[1]]$has.time.comp
    tsdim <- as.integer(spatialdim + has.time.comp)

    getrange <- function(x)
      if (x$grid) rbind(x$x[1, ], x$x[1, ] + x$x[2, ] * (x$x[3, ] - 1))
      else apply(x$x, 2, range)
    rangex <- sapply(neu, getrange)

    ## falls mehrere datasets:
    if (ncol(x[[1]]$x) > 1 || is.null(x[[1]]$dist.given) || !x[[1]]$dist.given){
      rangex <- t(rangex) 
      base::dim(rangex) <- c(length(rangex) / spatialdim, spatialdim)
    }
    rangex <- apply(rangex, 2, range)

    getmindistSq <- function(x) {
      if (x$grid) sum(x$x[2,]^2)
      else if (nrow(x$x) < 2) NA
      else if (nrow(x$x) <= mindist_pts) min(dist(x$x))
      else min(dist(x$x[sample(nrow(x$x), mindist_pts), ]))
    }

    if (has.time.comp &&
        any(sapply(neu, function(x) x$T[2]) <= RFopt$nugget$tol))
      stop("the step of the temporal grid component ",             
           if (any(sapply(neu, function(x) x$T[2]) == 0))
             "equals zero." else "is smaller than nugget tolerance 'tol'.")
    
    if (any(sapply(neu, function(x) x$grid && any(x$x[2,]<=RFopt$nugget$tol))))
          stop("the step of some spatial grid component ",             
               if (any(sapply(neu, function(x) x$grid && any(x$x[2,] == 0))))
                 "equals zero." else "is smaller than nugget tolerance 'tol'.")
  
    zaehler <- 0

    repeat {
      mindist <- sqrt(min(sapply(neu, getmindistSq)))      
      if (is.na(mindist)) mindist <- 1 ## nur 1 pkt gegeben, arbitraerer Wert
      if (mindist <= RFopt$nugget$tol) {
        if (!RFopt$general$allowdistanceZero)
          stop("Distance with value 0 identified -- use allowdistanceZero=T?")
        if ((zaehler <- zaehler + 1) > 10)
          stop("unable to scatter point pattern")
        for (i in 1:length(neu)) if (!neu[[i]]$grid)
          neu[[i]]$x <- neu[[i]]$x + rnorm(length(neu[[i]]$x), 0,
                                           10 * RFopt$nugget$tol)
      } else break;
    }
  
    xdimOZ <- ncol(neu[[1]]$x)
  }

 
  if (is.null(coordnames))
    coordnames <- SystemCoordNames(locinfo=neu[[1]], RFopt=RFopt) 
  
  restotal <- sapply(neu, function(x) x$restotal)
  ldata <- sapply(data, length)
  
  if (missing(model)) {
    if (!missing.further)
      stop("'model' is not given, but 'further.models'")
   } else {
    model.vdim <- info$newmodel$vdim
    if (model.vdim != 0) {
      if (length(vdim) == 0) vdim <- model.vdim
      else if (vdim != model.vdim)
        stop("given multivariate dimension differs from the dimensions expected by the model")
   }

    if (is.null(newmodel)) stop("new model not found.", CONTACT)

    if (!missing.further) {
      if (!missing(params) && length(params) > 0) {        
        nf <- names(further.models)
##  print("E")## danach lange

        for (m in 1:length(further.models))
          if (!is.null(further.models[[m]]) &&
                   !is.numeric(further.models[[m]])) {
            M <- further.models[[m]]
            further.models[[m]] <-
              PrepareModel2(model=model, ..., params=params,add.na=add.na,
                            further.model=M, arg = nf[m],
                            return_transform = FALSE,
                            x=neu, data=orig.data[[1]])$model
          }        
##  print("E")
      } else {
##   print("Ex")## danach lange

       for (m in 1:length(further.models)) {
          if (!is.null(further.models[[m]]) && !is.numeric(further.models[[m]]))
            further.models[[m]] <-
              PrepareModel2(further.models[[m]], ..., add.na=add.na,
                            return_transform = FALSE,
                            data=orig.data[[1]])$model
        }
##   print("Ex")
      }
   }
   }


  if (length(vdim) == 0) {
    repetitions <- sapply(data, ncol)
    if (min(repetitions) == 1) vdim <- 1 ## ggT==1 wuerde reichen
  } else repetitions <- as.integer(ldata / (restotal * vdim))

 
  if (length(vdim) == 0) {
    ## Verbesserung kann nur erreicht werden, wenn vorher
    ##            keine x-Koord fuer Modell bekannt waren. 
    if (missing(model) || !missing.x) {
      stop("Not detectable which columns are data columns.\nPlease set RFoptions(varnames=) or RFoptions(varidx=) accordingly.")
    }

    UD <- UnifyData(model=model, x=neu, data=data,
                    distances=distances, RFopt=RFopt,
                    mindist_pts=mindist_pts,
                    dim=dim, allowFirstCols=allowFirstCols,
                    vdim = vdim,
                    params=params,
                    further.models=orig.further.models, ...)
    UD$data.info$data.info <- info$data.info
    return(UD)
  } else {
    if (PL > PL_SUBIMPORTANT ||
        (PL >= PL_IMPORTANT && info$data.info != "safe")) {
      DetectionNote("Using columns ", paste(info$is.data, collapse=","),
              " as data columns",
              if (any(repetitions) > 1 && length(data) < 20)
                paste("with", paste(repetitions, collapse=",") , "repetitions"),
              ".")
    }
  }
 
  if (any(repetitions)==0) stop("no or not sufficiently many data are given")
  if (!is.na(vdim) && any(ldata != repetitions * restotal * vdim))
    stop("mismatch of data dimensions")

  vrep <- repetitions * vdim

  for (i in 1:length(data)) {
    base::dim(data[[i]]) <- c(ldata[i] / vrep[i], vrep[i])
  }


  ##  Print(vdim, vrep, repetitions); stopifnot(vdim ==2)
#  Print(repetitions)

  return(list(
    ## coords = expandiertes neu # #
    model = if (missing(model)) NULL else newmodel,
    orig.model = if (missing(model)) NULL else model,
    further.models = further.models,
    transform = info$newmodel$transform,
 ##   dimdata=dimdata, isRFsp = isRFsp,  # 3.2.0 : not given anymore
## oben auch noch loeschen!!
 ##  len = len,    # 3.2.0 : not given anymore
    data=data,
    RFsp.info = RFsp.info,
    coords = neu,    
    dist.given=dist.given,
    spatialdim=spatialdim,
    tsdim=tsdim,
    rangex = as.matrix(rangex),
    coordunits=coordunits,
    has.time.comp = has.time.comp,
    matrix.indep.of.x.assumed = matrix.indep.of.x.assumed,
    mindist = mindist,
    xdimOZ = xdimOZ,
    vdim = vdim,
    coordnames=coordnames,
    varnames=varnames,
    data.info = info,
    repetitions = repetitions,
    C_coords = C_coords,
    orig.trend = list(...)$trend,
    DataNames = info$M$data.coordnames ## info$M might be NULL already;
    ##                                       from model only!
  ))
}

