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

rfempirical <- function(x, y = NULL, z = NULL, T = NULL, data, grid,
			distances, vdim,
                        alpha=VARIOGRAM, ...
			) {
  
  ## repetition is last dimension
  
  ## bin centers will be a vector of scalar distances (in cylinder coords, e.g.)
  ## for the angles: start always with the first on negative angle, continue
  ##                 counter clockwise [0, 2pi]
  ## in 3 d the third angle is zero if vector in the (x, y) plane, positive
  ##                 angle starts with points above plane
  
  ## make sure that exactly one negative value appears, and that zero is
  ## added if bin starts with a positive value
  if ((is(data, "RFsp") || isSpObj(data)) && !missing(x))
    stop("x, y, z, T may not be given if 'data' is of class 'RFsp' or an 'sp' object")

  
  ## to do: distances
  if (!missing(distances) && length(distances)>0)
    stop("option distances not programmed yet.")

  ##   seq(from=0, to=deltaT[1], by=deltaT[2])}
  byT <- 2
  toT <- 1
 
  RFopt <- internalRFoptions(COPY=FALSE, ...)
  ## kein  if (!hasArg("COPY")) on.exit(optionsDelete(RFopt)) da nichts zu loeschen
  
  
  varunits <- RFopt$coords$varunits
  call <- match.call()

  Z <- UnifyData(x=x, y=y, z=z, T=T, distances=distances, grid=grid,
		 RFopt = RFopt,
		 data=data,  ## allowFirstCols=FALSE, 20.12.17 -- warum war dies
		 ##             gesetzt??
		 vdim=if (missing(vdim)) NULL else vdim)
  
  grid <- sapply(Z$coords, function(z) z$grid)
    
  if (Z$dist.given) stop("option distances not programmed yet.")
  
  if (missing(vdim) || length(vdim) == 0) {
    vdim <- if (!is.na(Z$vdim)) Z$vdim else 1
  } else {
        if (!is.na(Z$vdim) && vdim!=Z$vdim)
	  warning("given multivariate dimension 'vdim' does not match multivariate dimension of the data")
  }


  grid <- sapply(Z$coords, function(z) z$grid)
  data <- boxcoxIntern(Z$data, vdim=vdim, ignore.na=TRUE)  

  totpts <- sapply(Z$coords, function(z) z$totpts)
  spatialdim <- Z$spatialdim
  repetitions <- Z$repetitions
  sets <- length(Z$data)
  dist.given <- FALSE
  for (i in 1:sets) {
    ## ginge auf C ohne kopieren:
    dist.given <- dist.given || Z$coords[[i]]$dist.given
    dataX <- data[[i]]
    dim.data <- c(totpts[i], vdim, repetitions[i])    
    dim(dataX) <- dim.data    
    if (vdim > 1 && repetitions[i] > 1) {
      dataX <- aperm(dataX, c(1, 3, 2)) ## now: coords, repet, vdim
      dim(dataX) <- c(dim.data[1] * dim.data[3], dim.data[2])
      variance <- cov(dataX)  
    } else {
      dim(dataX) <- if (vdim == 1) prod(dim.data) else dim.data[1:2]    
      variance <- var(dataX)
    }
    rm(dataX)
  }
  
  pseudo <- alpha > 0
#  Print(pseudo, alpha,  PSEUDO, PSEUDOMADOGRAM)

  deltaT <- RFopt$empvario$delta
  deltaTgiven <- all(deltaT > 0)
  bins = RFopt$empvario$bins
  if(length(bins)==0) bins <- 20
  if (length(bins) == 1) {
    ## automatic bin depending on coords
    xx <- Z$coords[[1]]$x
    if(grid[1]) {
      maxi <- if (ncol(xx) == 1 && length(Z$coords[[1]]$T) == 0) {
                    bins * if (deltaTgiven) deltaT[byT] else 1
              } else max(xx[XSTEP + 1, ] * xx[XLENGTH + 1, ]) / 2
      bins <- seq(0, maxi, len = bins)
    } else {
      bins <- seq(0, sqrt(sum((apply(xx, 2, max)-apply(xx, 2, min))^2))/2,
		 len = bins)
    }
    if (RFopt$basic$printlevel >= PL_SUBIMPORTANT)
      message("Bins in RFvariogram are chosen automatically:\n", 
	      paste(signif(bins, 2), collapse=" "))
  }
  bins <- prepareBins(bins)
  if (length(bins) < 3 || !all(is.finite(bins)))
    stop("'bins' not finite or too few values")
  if (any(diff(bins)<=0)) stop("'bins' must be a strictly increasing sequence")
  
  centers <- pmax(0, (bins[-1] + bins[-length(bins)])/2)
  n.bins <- length(bins) - 1

#  Print(centers, bins)

#Print(phi0, phigiven)
  
  fft <- RFopt$empvario$fft && grid[1] && all(grid == grid[1]) &&
    alpha %in% c(VARIOGRAM, PSEUDO) && !dist.given
  
  if (alpha == COVARIANCE && dist.given && vdim > 1)
    warning("Fully symmetric covariance structure assumed as only distances are given.")

##  Print(fft, pseudo, alpha)
  
  has.time.comp <- Z$has.time.comp

  n.phi <- RFopt$empvario$nphi # 0 if automatic
  n.theta <- RFopt$empvario$ntheta # 0 if automatic
  phigiven <-  !dist.given && spatialdim > 1 && n.phi > 1
  thetagiven <- !dist.given && spatialdim > 2 && n.theta > 1
  basic <- !(has.time.comp || phigiven || thetagiven)
  
  phi <- if (!phigiven) c(0, 0) else c(RFopt$empvario$phi0, n.phi)
  n.phi <- max(1, phi[2])
  theta <- if (!thetagiven) c(0, 0) else  c(RFopt$empvario$theta0, n.theta)
  n.theta <- max(1, theta[2])
   
  if (has.time.comp) {
    T.start  <- sapply(Z$coords, function(x) x$T[XSTART + 1])
    T.step <- sapply(Z$coords, function(x) x$T[XSTEP + 1])
    T.len  <- sapply(Z$coords, function(x) x$T[XLENGTH + 1])
    if (sets > 1) {
      if (any(abs(diff(diff(T.step))) > 1e-15))
	stop("only data sets with the same time step allowed") #generalise todo
    }
    T <-  c(0, T.step[1], max(T.len))
  } else {
    T <-  c(1, 1, 1)
  }
  
  timeComponent <- T[XLENGTH + 1] > 1 && deltaTgiven ## T[XLENGTH + 1] > 1 impliziert time
 
  if (timeComponent) {
    stepT <-  deltaT[byT] / T[XSTEP + 1]
    if (abs(stepT - as.integer(stepT)) > stepT * 1e-15 || stepT <= 0)
      stop("deltaT not a multiple of the distance of the temporal grid")
    stepT <- as.integer(stepT + 0.5) ## in Gittereinheiten ; +0.5 for rounding 
    nstepT <-as.integer(min(T[XLENGTH + 1],
                            as.integer(deltaT[toT] / T[XSTEP+1])) / stepT)
    ## nstepT == n <=> set of temporal bin centers has length n + 1 elements,
    ## i.e. nstep == #( {temp bin centers} \ { 0 } )
  } else {
    stepT <- 1
    nstepT <- 0
  }
  n.delta <- 1 + nstepT
  
  n.phibins <- n.phi 
  dplt <- (!fft && !basic && !has.time.comp) ||
          ((pseudo || timeComponent) && phi[2]>0)
  if (dplt) n.phibins <- n.phibins * 2

  ##  Print(fft, basic, has.time.comp, pseudo, timeComponent, deltaT, deltaTgiven, T, phi, n.phi, n.phibins, dplt); 
  
  totalbinsOhnevdim <- as.integer(n.bins * n.phibins * n.theta * n.delta)
  totalbins <- totalbinsOhnevdim * vdim^2
  
  phibins <- thetabins <- Tbins <- NULL
  
  if (timeComponent) Tbins <- (0:nstepT) * deltaT[byT] 
  if (phi[2] > 0) phibins <- phi[1] + 0 : ((n.phibins - 1)) * pi / n.phi

  if (n.theta > 1)
    thetabins <- theta[1] + (0 : (n.theta-1) + 0.5) * pi / n.theta
  
  dims <- c(bins=n.bins, phi=n.phibins, theta=n.theta, delta=n.delta,
	    vdim=rep(vdim, 2))
  
  empirical.sd <- NULL

  if (fft) { # TO DO!!
    ## to do: das liest sich alles irgendwie komisch
    maxspatialdim <- 3
    
    if (Z$spatialdim > maxspatialdim)
      stop("fft does not work yet for spatial dimensions greater than ",
	   maxspatialdim)
    
    empirical <- N <- meanVar <- 0
    for (i in 1:sets) {
      xx <- Z$coords[[i]]$x
      if (ncol(xx)<maxspatialdim)  # not matrix(0, ...) here!
        ##                              since x is a triple
        xx <- cbind(xx, matrix(1, nrow=nrow(xx), ncol=maxspatialdim-ncol(xx)))
      TLen <- if (has.time.comp) Z$coords[[i]]$T[XLENGTH + 1] else 1
      neudim <- c(xx[XLENGTH + 1, ], if (has.time.comp) TLen)
            
      ## last: always repetitions
      ## last but: always vdim
      ## previous ones: coordinate dimensions
      dim(data[[i]]) <- c(neudim,vdim, length(data[[i]]) / vdim / prod(neudim))

      ## to achieve a reflection in x and z instead of y we transpose the
      ## array
      crossvar <- doVario(X=data[[i]], asVector=TRUE, pseudo=pseudo,
                          has.time.comp=has.time.comp)
 
                                        #
#crossvar =  List of 3
# $ : num [1:80000] -3.79426e-13 2.92282e+01 5.66963e+01 7.83819e+01 8.18708e+01 ...
# $ : num [1:80000] 100 99 98 97 96 95 94 93 92 91 ...
# $ : NULL

#      Print(pseudo)    
#      str(crossvar)
   #   print(crossvar)
   
      
#      Print(pseudo, range(crossvar[[1]], na.rm=TRUE), range(crossvar[[2]], na.rm=TRUE), if (length(crossvar[[3]]) > 0) range(crossvar[[3]], na.rm=TRUE))

      
      back <- .Call(C_fftVario3D, as.double(xx), 
		    crossvar[[1]],
                    crossvar[[2]],
                    crossvar[[3]],
		    as.double(bins), as.integer(n.bins), 
		    as.integer(TLen), 
		    as.integer(stepT), as.integer(nstepT),       
		    as.double(phi), 
		    as.double(theta), 
		    as.integer(repetitions[i]),
		    as.integer(vdim),
		    totalbinsOhnevdim,
		    as.logical(pseudo),
                    as.double(RFopt$empvario$tol))

      if (maintainers.machine() && totalbinsOhnevdim == n.bins) {
        ## 16.11.20 for debugging only
        backX <- .Call(C_fftVario3DX, as.double(xx), 
		    crossvar[[1]],
                    crossvar[[2]],
                    crossvar[[3]],
		    as.double(bins), as.integer(n.bins), 
		    as.integer(TLen), 
		    as.integer(stepT), as.integer(nstepT),       
		    as.double(phi), 
		    as.double(theta), 
		    as.integer(repetitions[i]),
		    as.integer(vdim),
		    totalbinsOhnevdim,
		    as.logical(pseudo) )

        diff <- abs(back[, 1:2] - backX)
        if (any(diff > 0)) {
          Print(pseudo) ## ok
          print(cbind(back[, 1:2],NA, backX, back[,1:2]-backX)) ## ok
          stop("RFemp error: pls contact author")
        }
      }
       
      ## the results are now reformatted into arrays
      ## the angles are given in clear text
      N <- N + back[, EV_FFT_N + 1]
      empirical <- empirical + back[, EV_FFT_EV + 1]# back contains only sums, not averages
      if (pseudo) meanVar <- meanVar + back[, EV_FFT_VAR + 1]
    }
    empirical <- empirical / N ## might cause 0/0, but OK
    if (pseudo) meanVar <- meanVar / N
    N <- as.integer(round(N))   

    if (pseudo) {
      dim(meanVar) <- dims
                                        #    Print(meanVar[1, 1, 1, 1, , ])
      meanVar <- rep(as.vector(meanVar[1, 1, 1, 1, , ]),
                     each = prod(rev(dims)[-1:-2]))
      dim(meanVar) <- dims
      empirical <- 2 * empirical + meanVar
    empirical.sd <- NULL ## kein korrekter Wert -- TODO
    }

  } else { ## ! fft
    ## #####################################################################
    ##
    ## MARTINS CODE WENN FFT == FALSE
    ##
    ## #####################################################################
    
    if (basic) {
      N <- empirical.sdSq <- empirical <- 0

      for (i in 1:sets) {

        ## ACHTUNG: nicht Z$C_coord !! (TODO, damit es schneller laeuft)
        back <- .Call(C_empirical, 
                        Z$coords[[i]]$x, ## Z definition
                        as.integer(spatialdim),
                        data[[i]],
                        as.integer(repetitions[i]), as.integer(grid[i]), 
                        as.double(bins), as.integer(n.bins),
                      as.integer(vdim),
                      as.double(alpha),
                      Z$coords[[i]]$dist.given)
##        Print("back")
        n.new <- back[, EV_N + 1]
	N <- N + n.new
	dummy <- back[, EV_EV + 1]    
        dummy[(is.na(dummy) & (centers==0)) | n.new == 0] <- 0
        empirical <- empirical + dummy * n.new
	
        dummy <- back[, EV_SDSQ + 1]
        dummy[n.new == 0] <- 0
        empirical.sdSq <- empirical.sdSq + dummy * n.new	
      }
      dummy <- N != 0
      empirical[dummy] <- empirical[dummy] / N[dummy]
      empirical.sd <- sqrt(empirical.sdSq / N)      
      rm("back")
    } else { ## anisotropic space-time
      ## always transform to full 3 dimensional space-time coordinates
      ## with all angles given. Otherwise there would be too many special
      ## cases to treat in the c program. However, there is some lost
      ## of speed in the calculations...
      if (!is.matrix(data[[1]])) stop("Bug")
    
      for (i in 1:sets) {
	coords <-  Z$coords[[i]]

        xx <- coords$x
        if (!is.matrix(xx)) stop("coordinates are not given by a matrix")
	if (ncol(xx)<4)  # not matrix(0, ...) here! since x could be a triple
          xx <- cbind(xx, matrix(1, nrow=nrow(xx), ncol=4-ncol(xx)))   
	
	## x fuer grid und nicht-grid: spalte x, y, bzw z
	N <- empirical.sdSq <- empirical <- 0

	back <-
          .Call(C_empvarioXT,
		xx,
		as.double(if (length(coords$T)>0) coords$T else rep(1,3)),
		data[[i]],
		as.integer(repetitions[i]),
		as.integer(grid[i]), 
		as.double(bins), as.integer(n.bins), 
		as.double(phi[1:2]), 
		as.double(theta[1:2]), 
		as.integer(c(stepT, nstepT)), 
		## input : deltaT[toT] max abstand, deltaT[byT]: echter gitterabst.
		##   c   : delta[toT]: index gitterabstand, deltaT[byT]:#of bins -1
		##                   (zero is the additional distance)
		as.integer(vdim),
		as.double(alpha),
                Z$coords[[i]]$dist.given
		)
	
	N <- N + back[, EV_N + 1]
	dummy <- back[, EV_EV + 1]    
        dummy[(is.na(dummy) & (centers==0)) | back[, EV_N + 1] == 0] <- 0
        empirical <- empirical + dummy * back[, EV_N + 1]
	
        dummy <- back[, EV_SDSQ + 1]
        dummy[back[, EV_N + 1] == 0] <- 0
        empirical.sdSq <- empirical.sdSq + dummy^2 * back[, EV_N + 1]
	
	rm("back")

	if (FALSE) {
	  if (!has.time.comp && vdim == 1) {
	    ## vario is symmetric in phi;
	    ## so the number of phi's can be halfened in this case
	    dim(empirical) <- dims
	    dim(N) <- dims
	    dim(empirical.sdSq) <- dims
	    
	    if (dims[2] > 1) {
	      dims[2] <- as.integer(dims[2] / 2)
	      half <- 1 : dims[2]
	      N <- N[, half,,,,, drop=FALSE] +N[, -half,,,,,drop=FALSE]
	      empirical <- empirical[, half, , , , , drop=FALSE] +
	      empirical[, -half, , , , , drop=FALSE]
	      empirical.sdSq <- empirical.sdSq[, half, , , , , drop=FALSE] +
	      empirical.sdSq[, -half, , , , , drop=FALSE]
	      phibins <- phibins[half]
	    }
	  }
	} ## end false

	
      } ## sets
	
      idx <- N > 1 & !is.nan(empirical) & empirical != 0 
      evsdSq <- empirical.sdSq[idx] / N[idx]
      
      if (any(evsdSq < -1e-14)) {
	Print(idx, N[idx] - 1, empirical.sdSq[idx], #
	      empirical.sdSq[idx] / (N[idx] - 1), empirical)
	warning(paste(evsdSq))
      }
      evsdSq[evsdSq < 0] <- 0
      empirical.sd[idx] <- sqrt(evsdSq)   
      empirical.sd[!idx] <- NaN
    }
 
      
    ## ################################################################
    ##
    ## END OF MARPINS CODE WENN FFT == FALSE
    ##
    ## ################################################################

  } # !fft

##  Print(empirical , dims)

  dim(empirical) <- dims

  
  dim(N) <- dims
  if (!is.null(empirical.sd)) dim(empirical.sd) <- dims

  name <- list()
  namedim <- names(dims)
    for (i in 1:length(dims)) {
      name[[i]] <-
      if (namedim[i] %in% c("vdim1", "vdim2")) {
	if (length(Z$varnames) == 0) NULL
	else rep(Z$varnames, length.out=dims[i])
      } else if (namedim[i] != "bins") paste(namedim[i], 1:dims[i], sep="")  
    }
  dimnames(empirical) <- name
  ##  {} else names(empirical) <- Z$varnames[1]

  dim <- has.time.comp + spatialdim
  l <- list(centers=centers,
            empirical=empirical,
            var=variance,
            sd= empirical.sd,
            n.bin=N,
            phi.centers=phibins,
            theta.centers=thetabins,
            T=Tbins,
	      vdim = vdim,
            coordunits = rep(Z$coordunits, length=dim),
            dim = dim,
            varunits = varunits,
            alpha=alpha
            )
  if (RFopt$general$spConform) l <- do.call("new", c(list(CLASS_EMPIR), l))
  else class(l) <- CLASS_EMPIR

  return(l)  
} # function rfempirical


## ############################################
## END OF MAIN FUNCTION 
## ############################################



doVario <- function(X, asVector=FALSE, pseudo=FALSE, has.time.comp=FALSE) {
  dimX <- dim(X)
  idx.repet <- length(dimX) 
  idx.vdim <- length(dimX) - 1
  
  d <- length(dimX) - 2## last two dimensions are repet & vdim
  twoD <- dimX[3] == 1
  n <- d + pseudo 
  len<- 2^(n-1)
 
  X_list <- as.list(rep(NA, len))
  X_list[[1]] <- X
  
  ##reflect the data, carefully with time reflection
  refl.order <- if(has.time.comp && !pseudo) c(1,3,4) else c(1,3,2)
  
  j <- 2
  for (i in 1:(n-1)) {
    for (k in 1:(2^(i-1))) {
      X_list[[j]] <- reflection(X_list[[k]], refl.order[i])
      j <- j + 1
    }      
  }
  
  ## to do the crossvariogram
  
  ## decide which blocks are needed
  blockidx <- rep(FALSE, 8)
  if(!has.time.comp && !pseudo){
    blockidx[1:(if (twoD) 2 else 4)] <- TRUE ## else 3 D
  } else if(has.time.comp && pseudo) {
    stop("Time component is not compatible with Pseudo variogram")
  } else { # ((has.time.comp && !pseudo) || (!has.time.comp && pseudo))
    blockidx[if (twoD) c(1:2, 5:6) else 1:8] <- TRUE
  }

#  Print("do", pseudo,asVector, len, blockidx)
    
  numbers <- cubes <- array(dim=c(dimX[1:d], len, dimX[idx.repet],
                                  rep(dimX[idx.vdim], 2)))
  meanVar <- if (pseudo) cubes else NULL
  for (i in c(1:len)){
    crossvar <- crossvario(X_list[[i]], pseudo=pseudo, dummy=!blockidx[i])
#    Print(crossvar)
    if (has.time.comp) {
      cubes[,,,,i ,,,] <- crossvar[[1]]
      numbers[,,,,i ,,,] <- crossvar[[2]]      
      if (pseudo) meanVar[,,,,i ,,,] <- crossvar[[3]]
    } else {
      cubes[,,,i ,,,] <- crossvar[[1]]
      numbers[,,,i ,,,] <- crossvar[[2]]
      if (pseudo) meanVar[,,,i ,,,] <- crossvar[[3]]
    }
  }
  
  if (asVector)
    return(list(as.vector(cubes), as.vector(numbers), as.vector(meanVar)))
 
  
  ##revert the reflection ## currently not used as asVector
  i<- n - 1
  for (i in (n-1):1) {
    parts<- len / (2^i)      
    positions <- 2^(i - 1)       
    for (j in 1:parts) {
      for (k in 1:positions) {
	idx <- 2* positions * j- positions + k
	if (has.time.comp) {
	  cubes[,,,,idx ,,,] <- reflection(cubes[,,,,idx ,,,], i)
	  numbers[,,,,idx ,,,] <- reflection(numbers[,,,,idx ,,,], i)
	  if (pseudo)
            meanVar[,,,,idx ,,,] <- reflection(meanVar[,,,,idx ,,,], i)
	} else {
	  cubes[,,,idx ,,,] <- reflection(cubes[,,,idx ,,,], i)
	  numbers[,,,idx ,,,] <- reflection(numbers[,,,idx ,,,], i)
	  if (pseudo)
            meanVar[,,,idx ,,,] <- reflection(meanVar[,,,idx ,,,], i)
	}
      }
    }
  }
  return(list(cubes, numbers, meanVar))
} 

crossvario <- function(f, pseudo = FALSE, dummy = FALSE) {
  d <- dim(f)
  idx.repet <- length(d) 
  idx.vdim <- length(d) - 1
  repetvdim <- c(idx.vdim, idx.repet)
  vdim <- d[idx.vdim]
  repet <- d[idx.repet]
  CVd <- c(d[-repetvdim], repet, vdim, vdim)
  if (dummy)
    return(list(array(NA, dim=CVd), array(NA, dim=CVd), array(NA, dim=CVd)))
  
  idx <- rep(TRUE, length(d) - 2)
  idx.data <- paste("[", paste(1, ":", d, collapse=", "), "]")
  idx.vario <- paste("[", paste(rep(",", length(d)-2), collapse=""), "r,i,j]")
  idx.w <- paste("[", paste(1, ":", d[-repetvdim], collapse=", "), "]")
  
  dim.coord <- 2 * d[-repetvdim]-1
  F <- If <- array(0, dim=c(dim.coord, d[repetvdim]))
  eval(parse(text=paste("If", idx.data, "<- !is.na(f)")))
  f[is.na(f)] <- 0
  eval(parse(text=paste("F", idx.data,  "<- f")))
  LIf <- list(If)
  LF <- list(F)
  
  nbvals <- Crossvario <- array(0, CVd)
  meanVar <- if (pseudo) array(0, CVd) else NULL
  
  for (i in 1:vdim) {
    for (j in 1:vdim) {
      for (r in 1:repet) {
	If <- do.call("[", c(LIf, idx, i, r))
	dim(If) <- dim.coord
	Ig <- do.call("[", c(LIf, idx, j, r))        
	dim(Ig) <- dim.coord
	F <- do.call("[", c(LF, idx, i, r))
	dim(F) <- dim.coord
	G <- do.call("[", c(LF, idx, j, r))
	dim(G) <- dim.coord

	if (!pseudo) {    
	  fftIfIg <- fft(If * Ig)
	  fftFG <- fft(F * G)
	  fftIfG <- fft(G * If)
	  fftIgF <- fft(F * Ig)   
	  z <- fft(Conj(fftFG) * fftIfIg
		   + Conj(fftIfIg) * fftFG
		   - Conj(fftIgF) * fftIfG
		   - Conj(fftIfG ) * fftIgF, inverse=TRUE)
	  N <- fft( Conj(fftIfIg) * fftIfIg, inverse=TRUE )
	} else {
	  F2 <- F^2
	  G2 <- G^2
	  fftIf <- fft(If)
	  fftIg <- fft(Ig)
          mV <- fft(Conj(fft(F2))* fftIg + Conj(fftIf) * fft(G2), inverse=TRUE)
	  z <- fft(- 2* Conj(fft(F)) * fft(G), inverse=TRUE)
	  ## N <- 2* fft(Conj(fftIf)*fftIg, inverse=TRUE)
	  N <- fft(Conj(fftIf)*fftIg, inverse=TRUE)
          w <- Re(mV) / (2 * prod(dim(N)))
          eval(parse(text=paste("meanVar", idx.vario, "<- w", idx.w)))
 	}
	
	w <- Re(z) / (2 * prod(dim(N))) 
	eval(parse(text=paste("Crossvario", idx.vario, "<- w", idx.w)))
  	eval(parse(text=paste("nbvals", idx.vario,
			      "<- Re(N", idx.w, ") / prod(dim(N))")))
      }
    }
  }  
  return(list(C=Crossvario, n=as.array(round(nbvals)), mV=meanVar))
}


prepareBins <- function(bins) {
  if(missing(bins)) return(NULL)
  if (bins[1] > 0) {
    if (getRFoptions(getoptions_="basic")$printlevel>1)
      message("empirical variogram: left bins border 0 added\n")
    bins <- c(0, bins)
  }
  if (bins[1]==0) bins <- c(-1, bins)
  if (bins[1] < 0) bins <- c(bins[1], bins[bins>=0])
  
  bins
}



