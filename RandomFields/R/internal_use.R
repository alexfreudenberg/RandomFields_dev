
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



########################################################################
## Hier sind Funktionen, die vermutlich kaum ein Anwender interessiert
## und wohl eher zum internen Gebrauch, wenn auch im NAMESPACE stehend;
##
## Ausnahme: rfGenerate* stehen in generatemodels.R
########################################################################

PrintModelList <-function (operators=FALSE, internal=FALSE,
                           newstyle=TRUE) {
   stopifnot(internal <=2, internal >=0, operators<=1, operators>=0)
    .C(C_PrintModelList, as.integer(internal), as.integer(operators),
       as.integer(newstyle))
    invisible()
}



showManpages <- function(path="/home/schlather/svn/RandomFields/RandomFields/man") {
  files <- dir(path)
  result <- character(length(files))
  for (i in 1:length(files)) {
    cat("\n\n\n\n\n")
    system(paste("grep \"\\examples\" -A 10000 ", path, "/", files[i], sep=""))
    result[i] <- readline(files[i])
  }
  result <- cbind(files, result)
  result <- result[result[, 2]!= "", ]
  result[result[, 2]=="t", 2] <- "To Do"
  result[result[, 2]=="s", 2] <- "dontshow"
  return(result)
}
# showManpages()

compareVersions <- function(path = "~/svn", package="RandomFields",
			    subpaths=c("cran", package)) {  
  path <- paste(path, package, subpaths[1], sep="/")
  l <- list(list(subpath=paste(path, package, subpaths[2], "R", sep="/"),
		 pattern="*R"),
            list(subpath=paste0(path, package, subpaths[2], "/src", sep="/"),
		 pattern="*\\.h"),
            list(subpath = paste0(path, package, subpaths[2], "/src", sep="/"),
		 pattern = "*\\.cc"))
  
  for (i in 1:length(l)) {
    subpath <- l[[i]]$subpath
    pattern <- l[[i]]$pattern
    files <- dir(path=subpath, pattern=pattern)
    
    for (f in files) {
      cmd <- paste("diff ",  path, "/", subpath, "/", f, " ",
		   subpath,"/", f, sep="")
      Print(cmd)  #
      system(cmd)
    }
  }
}

jpegORpdf <- function(graphics=NULL) {
  GD <- getOption("device")
  
  ispdf <- is.logical(all.equal(GD, "pdf"))
  isjpg <- is.logical(all.equal(GD, "jpeg"))
  

  if (!ispdf && !isjpg) {
    ispdf <- is.function(GD) && is.logical(all.equal(args(pdf), args(GD)))
    isjpg <- is.function(GD) && is.logical(all.equal(args(jpeg), args(GD)))
  } 
  
  if (is.null(graphics)) graphics <- getRFoptions(getoptions_="graphics")        
  file <- graphics$file
  if (file != "") {
    f <- strsplit(file, "\\.")[[1]]
    ext <- tolower(f[length(f)])
    is_pdf <- ext == "pdf"
    is_jpg <- ext == "jpeg" || ext == "jpg"
    if (is_pdf || is_jpg) graphics$file <- paste(f[-length(f)], sep=".")
   
    if (ispdf || isjpg) {
      if ((is_pdf && isjpg) || (is_jpg && ispdf)) stop("'extension' does not fit given device")
    } else {
      ispdf <- is_pdf
      isjpg <- is_jpg
    } 
  }

  if (ispdf) GD <- pdf
  else if (isjpg) GD <- jpeg

  graphics$ispdf <- ispdf
  graphics$isjpg <- isjpg
  graphics$ext <- c(".pdf", ".jpg")[isjpg + 1]
  graphics$GD <- GD

#  Print(graphics)

  return(graphics)
  
}

ScreenDevice <- function(height, width) {
  graphics <- jpegORpdf()    
  if (graphics$ispdf || graphics$isjpg) {
    if (graphics$filenumber == 99999)
      stop("number of pictures exceeds max. number of pictures")        

    file <- graphics$file
    if (file != "") {
      if (!graphics$onefile) {
        file <- paste(file,
                      formatC(format="d", flag="0", graphics$filenumber,
                              width=5),
                      graphics$ext, sep="")
        setRFoptions(graphics.filenumber = graphics$filenumber + 1)
      } else file <- paste(file, graphics$ext, sep="")
      ##            Print(file, dev.list())
      if (graphics$ispdf)
        graphics$GD(height=height, width=width, onefile=FALSE, file=file)
      else if (graphics$isjpg) graphics$GD(height=height, width=width,
                                           file=file, units="in", res=graphics$resolution)
      else stop("device not recognized by 'RandomFields'")
      return()
    } else {
      if (graphics$isjpg) stop("'file' in 'RFoptions' not set.")
    }
  }
  
  args <- names(as.list(args(graphics$GD)))    
  if (all(c("width", "height") %in% args) &&
      ( !any(c("file", "filename") %in% args)) || graphics$ispdf) {
    graphics$GD(height=height, width=width)        
    ##    Print("OK", height, width, RFoptions()$graphics, par()$cex)
    return()
  }
  Note("aspect_ratio")
}


#restore_par <- function(oldpar) {
#  do.call(graphics::par, oldpar)
#  graphics::par(cex = oldpar$cex) ## muss extra aufgerufen werden. noch mehr?
#}

ArrangeDevice <- function(graphics, figs, dh=2.8, h.outer=1.2,
                          dw = 2.5, w.outer=0.7) {
#  plot(1,1); return;  Print("ARRANGE", graphics)
  
  H <- curH <- graphics$height
  curW <- graphics$width
  if (is.finite(H) && H>0) {
    if (is.na(graphics$width) || graphics$width<=0.0) {
      H <- H * pmin(1, graphics$increase_upto[1] / figs[1],
                    graphics$increase_upto[2] / figs[2])
      DH <- H * dh / (dh + h.outer)
      HO <- H - DH
      curH <- figs[1] * DH + HO
      W <- H * (dw + w.outer) / (dh + h.outer)
      DW <- W * dw / (dw + w.outer)
      WO <- W - DW
      curW <-  figs[2] * DW + WO
    }

    graphics <- jpegORpdf(graphics)
    jpg_pdf <- graphics$isjpg || graphics$ispdf
    
   if (is.na(graphics$always_open_device)) {
      open <- interactive() || jpg_pdf
      setRFoptions(graphics.always_open_device = open)
    } else open <- graphics$always_open_device

    if (length(dev.list()) > 0 && open) {      
      close <- graphics$always_close_device
      if (is.na(close)) close <- jpg_pdf
      if (close) dev.off()
    }
    
    if (open) { # && raphics$grDefault) {      
      ScreenDevice(height=curH, width=curW)
    }
    
  } else warning("device could not be opened as 'height' (=", H,
                 ") is not positive.")
 
  if (graphics$split_screen) {
    if (length(dev.list()) > 0 && any(par()$mfcol != 1)) par(mfcol=c(1,1))
  } else {
    if (any(figs != 1) && length(dev.list()) > 0) par(mfcol=figs)
  }
  return(c(curH, curW))
}


#     p <- c("DSpat", "GeoGenetix", "LS2Wstat", "ProbForecastGOP", "scpm", "UncerIn2", "	AHMbook", "CompRandFld", "constrainedKriging", "geoR", "georob", "geozoning", "lgcp", "NLMR", "raptr", "spagmix", "spatsurv", "windfarmGA", "fractaldim", "geostatsp", "IsoriX", "prioritizr", "rpanel", "SAMCpack", "SpaDES", "SpaDES.core", "SpaDES.tools", "spatstat")

# for (i in p) if (!do.call("require", list(i))) install.packages(i)

#  find -maxdepth 1 -type d -exec cat {}/00install.out \; > xxxx



StartExample <- function(reduced = TRUE, save.seed=TRUE) {
  if (save.seed) {
    RFopt <- internalRFoptions()    
#    L <- list(se ed=RFopt$ba sic$s eed)
 #   if (length(RFopt$coords) != 0) L$coord_system= RFopt$coords$coord_system
 #   assign("RandomFields_options", envir=.RandomFields.env, L)
    assign("RandomFields_options", envir=.RandomFields.env, RFopt)
  }
# Print(ls( envir=.RandomFields.env), .RandomFields.env)
  if (!interactive() || hasArg("reduced")) {
    ## do not touch next lines
    ## REDUCED <- reduced, reduced
    REDUCED <- reduced
    setRFoptions(examples_reduced = REDUCED, RETURN=FALSE)
  }
}

FinalizeExample <- function() {
  ## print("finalize")
  if (!interactive()) {
    close.screen(all.screens = TRUE)  
    n <- length(dev.list()) - 2 ## otherwise R CMD check complains
    ##                             for missing graphic device
    if (n > 0) {
      for (i in 1:n) {      
        dev.off() ## OK
      }
    }
  }
#  Print(ls( envir=.RandomFields.env), .RandomFields.env)
  options <- get("RandomFields_options", envir=.RandomFields.env)
                                        #
#  Print(options)
  options$examples_reduced <- FALSE 
  RFoptions(list_ = options) ## internalRFoptions geht nicht, da
  ##                            printlevel, seed auch gesetzt werden muss
  #do.call("RFoptions", options)
  #print("done finalize")
}


maintainers.machine <- function() .Call(C_maintainers_machine)
maintainer.info <- function(...)
  if (maintainers.machine()) paste0(" [Maintainer info: ", paste0(...), "] ")
