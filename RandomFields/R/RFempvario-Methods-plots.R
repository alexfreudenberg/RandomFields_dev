## Authors 
## Martin Schlather, schlather@math.uni-mannheim.de
##
##
## Copyright (C) 2012 -- 2014 Alexander Malinowski & Martin Schlather
##               2015 -- 2017 Martin Schlather
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



## Methods for classes 'RFempVario' and CLASS_FITLIST  #######################

as.matrix.RFempVariog <- function(x, ...) {
  z <- as.array.RFempVariog(x, ...) 
  if (length(dim(z)) > 2)
    stop("the data set cannot be turned into a matrix. Try 'as.array'")
  z
}

as.array.RFempVariog <-
  as.array.RFspatialGridDataFrame <- function(x, ...) {
    z <- as.matrix(if (is.list(x)) x$empirical else x@empirical)    
    dim <- dim(z)
    dim <- dim[dim > 1]
    if (length(dim) == 1) return(as.vector(z))
    dim(z) <- dim
  }

as.vector.RFempVariog <-
  as.vector.RFspatialGridDataFrame <- function(x, ...) {
    as.vector(as.array.RFempVariog(x))
  }


summary.RFempVariog <- function(object, ...) {
  OP <- if (isS4(object)) "@" else "$"
  empirical <- do.call(OP, list(object, "empirical"))
  if (is.array(empirical)) {
    dims <- dim(empirical)
    dims <- dims[dims > 1]
    dim(empirical) <- dims
  }
  l <- list(centers=do.call(OP, list(object, "centers")),
            empirical=empirical)
  obj <- do.call(OP, list(object, "sd"))
  if (length(obj) > 0) {
    l$sd <- obj
    if (is.array(empirical)) dim(l$sd) <- dims
  }

  obj <- do.call(OP, list(object, "phi.centers"))
  if (length(obj) > 1) l$phi <- obj
  obj <- do.call(OP, list(object, "theta.centers"))
  if (length(obj) > 1) l$theta <- obj
  obj <- do.call(OP, list(object, "vdim"))
  if (obj > 1) l$vdim <- obj
  
  obj <- do.call(OP, list(object, "T"))
  if (length(obj)>0) l$T <- obj
  
  class(l) <- "summary.RFempVariog"
  l
}


print.summary.RFempVariog <- function(x, ...) {
  cat("Object of class '", CLASS_EMPIR, "'\n", sep="")
  str(x, no.list=TRUE) #
  invisible(x)
}

print.RFempVariog <- function(x, ...) {
  print.summary.RFempVariog(summary.RFempVariog(x, ...))
}


setMethod(f="show", signature=CLASS_EMPIR,
          definition=function(object) print.RFempVariog(object)) #

## coercion methods
setAs(CLASS_EMPIR, "list",
      def=function(from) {
        li <- list()
        for (name in slotNames(from)) {
          li[[name]] <- eval(parse(text=paste("from@", name, sep="")))
        }
        return(li)
      })

setAs(CLASS_FITLIST, CLASS_EMPIR, def=function(from) from@ev)

## plot method


plotRFempVariogUnbinned <- function(x, coordunits, varunits, varnames,
                                    legend, ..., plotmethod="plot") {
  stop("currently not used ") # !!
  graphics <- RFoptions()$graphics
  dots = mergeWithGlobal(list(...))
  dotnames <- names(dots)
  coords <- GridTopology2gridVectors(cbind(x@centers$x, x@centers$T))  
  if (length(coords)>1) {    
    coords[[1]] <- sort(unique((coords[[1]] - min(coords[[1]])) *
			       rep(c(-1, 1), each=length(coords[[1]]))))
    lab.names <- dimnames(x@centers$x)[[2]]
    
    if (!(x@centers$spatialdim==1 && x@centers$has.time.comp) &&
	!(x@centers$x[3,2]==dim(x@empirical)[2]))
      coords[[2]] <- sort(unique((coords[[2]] - min(coords[[2]])) *
				 rep(c(-1, 1), each=length(coords[[2]]))))
    else
      lab.names[2] <- "T"
  } else {
    x@centers <- coords[[1]]-min(coords[[1]])
    do.call(graphics::plot, args=c(dots, list(x=x, plot.nbin=FALSE)))
    return()
  }
  
  if (!("main" %in% dotnames)) {
    main <- "Variogram image plot"
    if (length(varnames)>0) main <- paste(main, "for", varnames)
    dots$main <- main
  }
  lab.names <- paste(lab.names, "-distance", sep="")
  
  idx <- lab.names != "T-distance"
  if (any(idx) && all(coordunits[idx] != ""))
    lab.names[idx] <-
  paste(lab.names[idx], " [", coordunits[idx], "]", sep="")
  if (!all(idx) && all(coordunits[!idx] != ""))
    lab.names[!idx] <-
  paste(lab.names[!idx], " [", coordunits[!idx], "]", sep="")
  
  if (!("xlab" %in% dotnames)) dots$xlab <- lab.names[1]
  if (!("ylab" %in% dotnames)) dots$ylab <- lab.names[2]
  if (!("xlim" %in% dotnames)) dots$xlim <- range(coords[[1]])
  if (!("ylim" %in% dotnames)) dots$ylim <- range(coords[[2]])
  
  idx1 <- coords[[1]] >= dots$xlim[1] & coords[[1]] <= dots$xlim[2]
  idx2 <- coords[[2]] >= dots$ylim[1] & coords[[2]] <= dots$ylim[2]
  coords[[1]] <- coords[[1]][idx1]
  coords[[2]] <- coords[[2]][idx2]
  
  dims <- dim(x@empirical)
  ev.plot <- matrix(x@empirical, nrow=dims[1], ncol=dims[2])
  ev.plot <- ev.plot[idx1, idx2]
  
  range.ev <- range(ev.plot)  #x@empirical
  col <- (if ("col" %in% dotnames) dots$col else
          default.image.par(NULL, NULL)$default.col)

  if (graphics$split.screen && legend) {
    scr <- split.screen(erase=FALSE, rbind(c(0,.85,0,1), c(.85,1,0,1)))
    screen(scr[1], new=FALSE)
  } else scr <- NULL  
  par(mar=c(2,2,3,1))

  dots$type <- NULL
  do.call(plotmethod,
	  args=c(dots, list(x=coords[[1]], y=coords[[2]], z=ev.plot)))
					#      xlab=lab.names[1],
					#      ylab=lab.names[2],
					#      main="Variogram image plot"
  if (legend) {
    stop("das ist doch quark?!")
    
    if (graphics$split.screen) {
      screen(scr[2], new=FALSE)
      par(mar=c(2,0,3,3))
      z.legend <- seq(range.ev[1], range.ev[2], length=length(col))
      image(y=z.legend, x=c(0,1), z=rbind(z.legend, z.legend), axes=F, col=col,
	    xlab="")
      axis(4)
      box()
    } else {
      my.legend(min(coords[[1]]), max(coords[[2]]), range(ev.plot),
		col=col, bg="white")
    }
  }
  
  if (graphics$close_screen) {
    close.screen(scr)
    scr <- NULL
  }
  return(invisible(scr))
}

RFplotEmpVariogramX <- function (x, y=y,...) stop("hier")

RFplotEmpVariogram <- function(x, model = NULL, nmax.phi = NA, nmax.theta = NA,
			       nmax.T = NA,
			       plot.nbin = TRUE, plot.sd=FALSE,
                               method,
			       boundaries = TRUE,
 #                              cloud = FALSE,
			       ..., all=FALSE) {

  n.points <- 150 ## number of locations, the theoretical curve is evaluated
  mar.phi <- c(1,0.5,1,0.5)
  
  dots = list(...)
#  str(x, max=2)
#  Print(dots)
  ##  print(model)

  dotnames <- names(dots)
  graphics <- RFoptions()$graphics
  OP <- c("$", "@")[1 + isS4(x)]
  newx <- list()
  fitmethod.names <- character(0)
  if (is(x, CLASS_FITLIST)) {
    if(length(do.call(OP, list(x, "ev")))==0)
      stop("The fit does not contain an empirical variogram.")

    for (i in c("autostart", "self", "plain", "sqrt.nr", "sd.inv",
                "internal1", "internal2", "internal3", "ml"))
      newx[[i]] <- do.call(OP, list(x, i))

    fitmethod.names <-  getRFfitmethod(newx, method, all=all)
    newx <- newx[fitmethod.names]
    x <- do.call(OP, list(x, "ev"))
  } else {
    if (!is(x, CLASS_EMPIR))
      stop("method only for objects of class '", CLASS_EMPIR, "' or '",
           CLASS_FITLIST, "'")
  }
  
  for (i in c("empirical", "sd", "centers", "phi.centers", "theta.centers",
              "n.bin", "coordunits", "vdim", "dim",
              "alpha")) ## emprical variogram/covariance/etc
    assign(i, do.call(OP, list(x, i)))
  T <- do.call(OP, list(x, "T")) ## 21.11.21 -> bug report da in assign
  ## irregulaeres Verhalten
  ## Print(T)
  
  varnames <- if (is.matrix(empirical)) dimnames(empirical)[[2]][1]
              else names(empirical)[1]
#  if(alpha < 0 && vdim > 1)  warning("theoretical madogram model is not available at the moment")

  params <- if ("params" %in% dotnames) dots$params else NULL
  if (!("type" %in% dotnames)) dots$type <- "b" 
  cex <- if ("cex" %in% dotnames) dots$cex else .8
  if (!("pch" %in% dotnames)) dots$pch <- 19
  if(!("xlim" %in% dotnames)) dots$xlim <- range(centers)
  ylim.not.in.dotnames <- !("ylim" %in% dotnames)
  xlab <- if ("xlab" %in% dotnames) dots$xlab else "distance"
  ylab.ev <- if ("ylab" %in% dotnames) dots$ylab
             else if (alpha == as.integer(alpha)) FCTN_TYPE_NAMES[alpha + 3]
             else paste0("alpha-pseudomadogram with alpha=", abs(alpha))
  main0 <- if ("main" %in% dotnames) dots$main
           else if (alpha == as.integer(alpha))
             paste(FCTN_TYPE_NAMES[alpha + 3], "plot")
           else paste0("alpha-pseudomadogram plot with alpha=", abs(alpha))
  if ("oma" %in% dotnames) oma <- dots$oma
  else {
    if (!is.null(main0)) oma.top <- 2 else oma.top <- 0
    oma.left <- 6
    oma <- c(4,oma.left,oma.top,0)
  }
 
  dots$params <- dots$cex <-  dots$main <- dots$xlab <- dots$ylab <-
    dots$oma <- NULL

  

  lab <- xylabs("bin centers", NULL, units=coordunits)

  
  has.sd <- !is.null(sd)
  
  n.phi <- min(nmax.phi, l.phi <- max(1,length(phi.centers)), na.rm=TRUE)
  n.theta <- min(nmax.theta, l.theta <- max(1,length(theta.centers)),
		 na.rm=TRUE)
  n.T <- min(nmax.T, l.T <- max(1,length(T)), na.rm=TRUE)

  if(!is.null(model)) {
    if (!is.list(model)) model <- list(model)
    if (!all(sapply(model, FUN=function(x) !isS4(x) || ## model as list
                             is(x, CLASS_CLIST) || ## model as RMmodel
                             is(x, "formula") ## model as formula
                    )))
      stop("model must be (a list of elements) of class 'CLASS_CLIST'")
    modelnames <-
      if(length(names(model)) > 0) names(model)
      else paste("model", 1:length(model))
    fitmethod.names <- c(fitmethod.names, modelnames)
    names(model) <- modelnames
    newx <- c(newx, model)
  }
  n.methods <- length(fitmethod.names)
  
  
  if (n.phi > 6 || n.theta > 3 || n.T > 3)
    message("'If you feel the picture is overloaded, set the parameters 'nmax.phi', 'nmax.theta' and 'nmax.T'")
  
  halfphi.sector <- pi/(2*l.phi)
  halftheta.sector <- pi/(2*l.theta)
  phi.angles <- c(halfphi.sector, 0, -halfphi.sector) 
  theta.angles <- seq(-halftheta.sector, halftheta.sector, len=5) # len=odd!
  if (n.phi>1 && boundaries) {
    phi.angles <- phi.angles * 0.96
    theta.angles <- theta.angles * 0.96
  } 
  
  TandV <- (n.T > 1 && vdim > 1) && graphics$split_screen
  
  if (vdim>1 && length(varnames)==0)
    varnames <- paste("v", 1:vdim, sep="")
  
  range.nbin <- range(c(0, n.bin), na.rm=TRUE)
  ylim.nbin <- range.nbin * c(1,1.2)

  n.col <- 1 + max(1, n.methods)
  col.model <- if ("col.model" %in% dotnames) rep(dots$col.model, len=n.col)
               else 2:n.col
  dots$col.model <- NULL
  n.col <- max(1, n.phi)
  col.v <- col <-
    if ("col" %in% dotnames) rep(dots$col, len=n.col) else 1:max(n.col)
  dots$col <- NULL
  
  if (n.methods > 0){
    dotsRFfit <- dots
    dotsRFfit$type <- "l"
    dotsRFfit$lwd <- 2
    ltyRFfit.v <- 1:n.methods
    dotsRFfit$lty <- NULL
    x.base <- seq(from = max(dotsRFfit$xlim[1],1e-3),
                             to = dotsRFfit$xlim[2], len = n.points)
    spatialdim <- dim - (length(T) > 0) ## ja nicht n.T > 1, da
    ## stille Zeitkomponent moeglich, die aber das Modell braucht

    x.spt <- x.space <-
      do.call(cbind, c(list(x.base), as.list(rep(0, spatialdim-1))))
    total <- n.points *  vdim^2

    ##  [range/mean], x, vdim, vdim, n.methods, n.T, n.theta, nph, 
    VALUES <- array(dim=c(if (n.phi>1 && boundaries) 2 else 1,
                          n.points, vdim, vdim, n.methods,
                          n.T, n.theta, n.phi))
  }

  par(cex=cex, xaxs="i")
  Screens <- if (TandV) c(n.T, n.theta) else c(n.T * vdim * vdim, n.theta)
  ArrangeDevice(graphics, Screens)

  if (graphics$split_screen) {
    all.scr <- scr <-
      split.screen(if (TandV) c(vdim, vdim) else Screens, erase=FALSE)
  } else  all.scr <- scr <- NULL

  for (v1 in 1:vdim) {
    for (v2 in 1:vdim) {
      if (TandV) {
	scr <- split.screen(erase=FALSE, Screens, all.scr[v1 + (v2 - 1) * vdim])
	all.scr <- c(all.scr, scr)          
      }
      main <- if (vdim == 1) {
                if (is.null(main0) || length(varnames)==0) main0
                else paste(main0, "for", varnames)
              } else {
                if (!TandV) main0
                else paste(main0, "for", varnames[v1], "vs.",  varnames[v2])
              }
      if (ylim.not.in.dotnames)
        dots$ylim <- range(empirical[,,,, v1, v2], na.rm=TRUE)
      for (iT in 1:n.T) {
        for (ith in 1:n.theta) {
          ## plot n.bin
          if (plot.nbin) {
            screen(scr[1], new=FALSE)
            par(oma=oma)
            scr2 <- split.screen(erase=FALSE, rbind(c(0,1,.2,1), c(0,1,0,.2)),
                                 screen=scr[1])
            all.scr <- c(all.scr, scr2)
            screen(scr2[2], new=FALSE)
            par(mar=c(0,.5,0,.5))
            for (iph in 1:n.phi) {
              if (n.phi > 1) col <- col.v[iph]
              if (iph==1) {
                plot(centers,
                     n.bin[ ,iph, ith, iT, v1, v2],
                     xlim=dots$xlim, ylim=ylim.nbin,
                     type=if (n.phi>1) "p" else "h",
                     col =if (n.phi>1) col else "darkgray", lwd=4,
                     pch=16, axes=FALSE, ylab="",
                     xlab = lab$x)
                box()
                at <- seq(range.nbin[1], range.nbin[2], len=3)
                if (ith==1)
                  axis(2, las=1, at=at, labels=format(at, scient=-1, dig=2),
                       outer=TRUE)
                else  axis(2, las=1, at=at, labels=FALSE)
                
                if (iT==n.T && (n.T > 1 || (v1==vdim && v2==vdim))) axis(1)
                if (ith==1) title(ylab="n.bin", line=5, outer=TRUE, adj=0)
                
              } else {
                points(centers, n.bin[ ,iph, ith, iT, v1, v2],
                       type="p", col=col, pch=16)
              }
            }
            screen(scr2[1], new=FALSE)
          } else {
            screen(scr[1], new=FALSE)
            par(oma=oma)
          }
          
          ## plot empirical
          ##if (ith==1) par(mar=c(0,6,1,1)) else par(mar=c(0,.5,1,.5))
          
          plotted.meth <- NULL  # needed for legend
          for (iph in 1:n.phi) {
            if (n.phi>1) col <- col.v[iph]
            par(mar=mar.phi)
            if (iph==1) {
              do.call(graphics::plot,
                      c(dots, list(x=centers,
                                   y=empirical[,iph,ith,iT,v1,v2],
                                        #ylim=ylim.ev, type=type, pch=19,
                                   col=col, axes=FALSE, ylab="", xlab=lab$x)))
              box()
              axis(2, las=1, labels=(ith==1), outer=(ith==1))
              if (!plot.nbin) axis(1)          
              if (l.theta > 1 || l.T > 1 || vdim > 1) {
                L <- c(if (!TandV && vdim > 1)
                         paste(varnames[c(v1,v2)], collapse=":"),
                       if (l.T>1) paste("T=",signif(T[iT], 3), " ", sep=""),
                       if (l.theta>1)
                         paste(sep="", "theta=", signif(theta.centers[ith],3)))
                legend("topleft", legend=paste(L, collapse=","))
              }
              if (ith == 1) title(ylab=ylab.ev, line=5, outer=TRUE)
              if (has.sd && plot.sd)
                legend("topright", bty="n", #col=NA, lwd=NA,
                       legend="arrows represent +-1 sd intervals")
              
            } else {  ## iph > 1
              do.call(graphics::points,
                      c(dots, list(x=centers,
                                   y=empirical[,iph,ith,iT,v1,v2],
                                   col=col))
                      )   #type=type,  pch=19
	    } ## for iph 

	    if (n.methods > 0) {
             
              if (v1 == 1 && v2 == 1) {              
                if (!is.null(phi.centers)) {
                  ## sehr genau Abschaetzung, indem mehrere (3)
                  ## Winkel angeschaut werden und dann der mittlere
                  ## Wert der Variogramme angeschaut wird, da ja auch das
                  ## emp. Variogramm ueber ein Winkelintervall gemittelt wird
                  if (is.null(theta.centers)) { # 2D
                    angles <- phi.centers[iph] + phi.angles
                    x.space <- cbind(rep(cos(angles), each=n.points),
                                     rep(sin(angles), each=n.points))
                  } else { # 3D
                    angles <- expand.grid(phi.centers[iph] + phi.angles,
                                          theta.centers[iph] + theta.angles)
                    phis <- angles[ ,1]
                    thetas <- angels[ ,2]
                    ct <- cos(thetas)
                    x.space <- cbind(rep(ct * cos(phis), each=n.points),
                                     rep(ct * sin(phis), each=n.points),
                                     rep(sin(thetas), each=n.points))
                  }
                  x.spt <- x.space <- x.base * x.space
                }
                if (length(T) > 0) x.spt <- cbind(x.space, T[iT])
              }
            
              for(i in 1:n.methods) {
                if (v1 == 1 && v2 == 1) {
                  values <- covETC(model = newx[[fitmethod.names[i]]],
                                   x = x.spt, grid = FALSE, params=params,
                                   internal.examples_reduced=FALSE,
                                   alpha=alpha)
                  if(length(values) == 0 || all(is.na(values)))
                    stop("internal calcualation error in RFplotEmpVario")

                  if (length(values) > total) { ## the mean over
                    ## a boundle of angles around the precise one
                    dim(values) <- c(n.points, length(values) / total,vdim,vdim)
                    values <-
                      if (n.phi>1 && boundaries) apply(values, c(1, 3, 4),range)
                      else apply(values, c(1, 3, 4), mean)
                  }
                  VALUES[, , , , i, iT, ith, iph] <- values
                }

                values <- VALUES[ , , v1, v2, i, iT, ith, iph]		
                if (n.phi>1 && boundaries) {
                  do.call(graphics::matplot,
                          args=c(dotsRFfit,
                                 list(x=x.base,
                                      y=if (is.matrix(values)) t(values)
                                        else values,
                                      xlab=lab$x, ylab=lab$y,
                                      add = TRUE, col=col.model[i], lty = 3)))
                } else {
                  do.call(graphics::points,
                          args=c(dotsRFfit, list(x=x.base, y = values, 
                                                 col=col.model[i],
                                                 lty = ltyRFfit.v[i])))
                }    
                if(iph == 1) plotted.meth <- c(plotted.meth, fitmethod.names[i])
	      } ## for i in n.methods
	    } ## if n.methods > 0	  
            
	    if (has.sd && plot.sd) {
	      sdnot0 <-  sd[ ,iph, ith, iT] != 0
	      arrows(centers[sdnot0],
		     empirical[sdnot0 ,iph, ith,iT] - sd[sdnot0,iph,ith,iT],
		     centers[sdnot0],
		     empirical[sdnot0 ,iph, ith,iT] + sd[sdnot0,iph,ith,iT],
		     code=2, angle=90, length=0.05, col=col)
	    }
	    
	  } # iph in nphi	  
          
	  pifactor <- signif((phi.centers[1:n.phi]) / pi, 2)
	  
	  len.mnames <- length(plotted.meth)
	  string.emp <- "empirical"
	  if(len.mnames > 0) {
	    labels <-
              if (n.phi>1) paste("\"phi=\",", rep(pifactor, each = len.mnames+1),
                                 ", pi, ", "\", ") else "\""
	    labels <- paste("c(",
			    paste("expression(paste(",labels,
				  rep(c(string.emp, plotted.meth), l.phi),
				  "\"", "))", collapse=","), ")")
	    labels <- eval(parse(text=labels))
            colors <- c(col[1], col.model[1:length(plotted.meth)])
	  } else {
	    labels <- if (n.phi > 1) 
			eval(parse(text=paste(
                                       "c(",
                                       paste("expression(paste(\"phi=\",", pifactor,", pi))",
                                             collapse=","),
                                       ")")))
            colors <- col.v
	  }
	  
	  if (l.phi > 1 || len.mnames > 0) {# && iT==1 && ith==1 # auskommentiert auf Sebs wunsch
                                        #
            legend("bottomright",
                   col=colors,
                   lwd=1,
                   pch=rep(c(19,rep(NA, len.mnames)), l.phi),
                   bty="n",
                   legend=labels,
                   lty = rep(c(1, if(len.mnames==0) NULL
                                  else ltyRFfit.v[1:len.mnames]),
                             l.phi))
          }
	  scr <- scr[-1]
	} # n.theta
      } # T
    } # vdim 1
  } # vdim 2
  dots$type <- NULL
  
  if (!is.null(main))
    do.call(graphics::title, args=c(dots, main=main, outer=TRUE))
  if (!is.null(xlab))
    do.call(graphics::title, args=c(dots, xlab=xlab, outer=TRUE))
  
  if (graphics$close_screen) {
    close.screen(all.scr)
    all.scr <- NULL
  }
  
  return(invisible(all.scr))
}

setMethod(f="plot", signature(x=CLASS_EMPIR, y="missing"),
          function(x, y, ...) RFplotEmpVariogram(x, ...))
