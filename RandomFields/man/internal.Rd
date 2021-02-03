\name{Internal functions}
\alias{rfGenerateModels}
\alias{rfGenerateConstants}
\alias{rfGenerateTest}
\alias{rfGenerateMaths}
\alias{checkExamples}
\alias{StartExample}
\alias{FinalizeExample}
\alias{Dependencies}\alias{ScreenDevice}
\alias{plotWithCircles}
\alias{showManpages}
\alias{RM_DEFAULT}
\alias{maintainers.machine}
\title{Internal functions}
\description{
 These functions are internal and should not be used.
}
\usage{
rfGenerateModels(package = "RandomFields", assigning,
                 RFpath ="~/svn/RandomFields/RandomFields",
                 RMmodels.file=paste(RFpath, "R/RMmodels.R", sep="/"),
                 PL = RFoptions()$basic$printlevel)

rfGenerateConstants(package="RandomFields", aux.package = "RandomFieldsUtils",
	            RFpath = paste0("~/svn/",package, "/", package),
                    RCauto.file = paste(RFpath, "R/aaa_auto.R", sep="/"),
	            header.source = c(if (length(aux.package) > 0)
                                  paste0("../../", aux.package,"/",
	                          aux.package, "/src/Auto", Aux.package, ".h"),
	                          paste0("src/Auto", Package,".h")),
	            c.source = paste0("src/Auto", Package, ".cc"))

rfGenerateTest(package = "RandomFields", files = NULL,
               RFpath = paste0("~/svn/", package, "/", package))

rfGenerateMaths(package = "RandomFields",
			    files = "/usr/include/tgmath.h",
			    do.cfile = FALSE,
                            ## copy also in ../private/lit
                            Cfile = "QMath",
                            Rfile = "RQmodels",
                            RFpath = paste0("~/svn/",package,"/", package))

ScreenDevice(height, width)

FinalizeExample()
StartExample(reduced = TRUE, save.seed = TRUE)
%Dependencies(pkgs = all.pkgs, dir = "Dependencies",
%             install = FALSE, check=TRUE, reverse=FALSE, package =
%                 "RandomFields")
showManpages(path="/home/schlather/svn/RandomFields/RandomFields/man")

maintainers.machine()
}
\arguments{
  \item{package,assigning, RFpath, RMmodels.file, PL}{internal}
  \item{aux.package,RCauto.file,header.source,c.source}{internal}
  \item{files}{internal}
  \item{Cfile, Rfile, do.cfile }{internal}
%  \item{pkgs, dir,install, check, reverse}{internal}
  \item{height,width}{window sizes}
  \item{reduced, save.seed, path}{internal}
 % \item{REGISTER, COVREGISTER, RELAX}{internal}
}

\me
\keyword{spatial}

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

## internal functions: no examples given

\dontshow{FinalizeExample()}}
