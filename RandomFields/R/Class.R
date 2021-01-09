
## Authors 
## Martin Schlather, schlather@math.uni-mannheim.de
##
##
## Copyright (C) 2012 -- 2015 Alexander Malinowski & Martin Schlather
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



##########################################################################
## classes for 2- and higher-dimensional data objects, based on ##########
## Spatial classes from 'sp'-package                            ##########


setClass("RFspatialGridDataFrame", contains ="SpatialGridDataFrame",
         representation(.RFparams="list") )
setValidity("RFspatialGridDataFrame", 
            function(object) {
              return(check.validity.n.vdim(object))
            })


setClass("RFspatialPointsDataFrame", contains ="SpatialPointsDataFrame",
         representation(.RFparams="list") )
setValidity("RFspatialPointsDataFrame", 
            function(object) {
              return(check.validity.n.vdim(object))
            })


## classes for 1-dimensional data objects ################################

setClass("RFgridDataFrame", 
         representation(data="data.frame", grid="GridTopology",
                        .RFparams="list"))
setValidity("RFgridDataFrame", 
            function(object) {
              if (length(object@grid@cells.dim) > 1)
                return("this class is only for 1-dimensional coordinates")
              if (nrow(object@data) != object@grid@cells.dim)
                return("data must have the same length as the grid length")
              return(check.validity.n.vdim(object))
            })

setClass("RFpointsDataFrame", 
         representation(data="data.frame", coords="matrix",
                        .RFparams="list"),
         prototype(data=data.frame(NULL), coords=NULL, .RFparams=list()))
setValidity("RFpointsDataFrame", 
            function(object) {
              if (nrow(object@data) != length(object@coords))
                return("data and coords must have the same length")
              if (ncol(object@coords)!=1)
                return("coords must have exactly 1 column, otherwise use class 'RFspatialPointsDataFrame'")
              return(check.validity.n.vdim(object))
            })




setClassUnion("RFspatialDataFrame",
              c("RFspatialGridDataFrame", "RFspatialPointsDataFrame"))
setClassUnion("RFdataFrame", c("RFgridDataFrame", "RFpointsDataFrame"))
setClassUnion("RFsp", c("RFspatialGridDataFrame", "RFspatialPointsDataFrame",
                       "RFgridDataFrame", "RFpointsDataFrame"))

check.validity.n.vdim <- function(object) {
  if (!all(c("n", "vdim") %in% names(object@.RFparams)))
    return("slot '.RFparams' must contain 'n' and 'vdim'")
  var.given <- !is.null(object@.RFparams$has.variance) &&
    object@.RFparams$has.variance
  nc <- (object@.RFparams$n + var.given) * object@.RFparams$vdim

  if (nc != ncol(object@data)) {
    stop("number of data at each location (=", ncol(object@data),
         ") does not match the expected ones (=", nc,
         ").\n  The latter is based on the information or assumption that there are\n  ",
         object@.RFparams$n, " repetition(s) of ",
         object@.RFparams$vdim, " variable(s)",
         if (var.given) " and the variances",
         " given.",
         if (object@.RFparams$n * object@.RFparams$vdim == 1) "\nEither wrong parameters have been given for 'RFspatialGridDataFrame' or 'RFspatialPointsDataFrame' or an explicite transformation of the data from 'sp' objects to 'RFsp' objects has not been performed with 'sp2RF'.")
  }
  return(TRUE)
}




## definition of class CLASS_CLIST
setClass(CLASS_CLIST, 
         representation(
                        # call='RMexp(var=1, sclae=1, Aniso=id, proj=id)
                        call = "language",
                        # name='RMexp'
                        name = "character",
                        # submodels=NULL, submodels=list(RMmodel1, RMmodel2)
                        submodels = "list",
                        # model specific parameter 
                        par.model = "list",
                        # var=1, scale=1, Aniso=id, proj=id 
                        par.general = "list"
                        )
         )

## rules for validity checking of CLASS_CLIST objects
isModel <- function(model) return(is.list(model) && is.character(model[[1]]) )
isRMmodel <- function(x) isS4(x) && is(x, CLASS_CLIST)
setValidity(CLASS_CLIST, 
            function(object){
              if (length(object@submodels) > 0 &&
                  !all(sapply(object@submodels, function(x) isRMmodel(x)
                                        # || is(x, "formula")
                                     )))
                return(paste("submodels must be of class '", CLASS_CLIST,
                             "' or a formula", sep=""))
                              
              return(TRUE)
            })


## definition of class CLASS_RM ################################
setClass(CLASS_RM, contains ="function",
         representation(
                        type = "character",
                        domain = "character",
                        isotropy = "character",
                        operator = "logical",
                        monotone = "character",
                        finiterange = "logical",
                        maxdim = "numeric",      # max possible dimension
                        simpleArguments = "logical",
                        vdim = "numeric"         # ??
                        )
         )


## definition of class 'RMmodelFit'
setClass(CLASS_SINGLEFIT, # contains=CLASS_CLIST,
         representation(
             model = CLASS_CLIST,
             modelAsList = "list",
             formel = "ANY",
             likelihood = "numeric",
             variab = "matrix",
             coordsystem = "character",
             param = "matrix",
             globalvariance  = "ANY",
             covariat = "ANY",
             hessian = "ANY",
             AIC = "numeric",
             AICc = "numeric",
             BIC = "numeric",
             residuals = "ANY"
             )
         )


setClass(CLASS_EMPIR, 
         representation(centers = "ANY",
                        empirical = "ANY",
                        var = "ANY",
                        sd = "ANY",
                        n.bin = "ANY",
                        phi.centers = "ANY",
                        theta.centers = "ANY",
                        T = "ANY",
                        vdim = "numeric",
                        coordunits = "character",
                        dim = "numeric",
                        varunits = "character",
                        call = "ANY",
                        alpha = "numeric"
                        )
         )


setValidity(CLASS_EMPIR, 
            function(object){
              if(!(is.null(object@call)) && !(is(object@call, "language")))
                return("slot 'call' must be NULL or of class 'language'")
              return(TRUE)
            })




setClass(CLASS_FITLIST, 
         representation(Z = "list",
                        ev=CLASS_EMPIR,
                        table = "data.frame",
                        n.variab = "integer",
                        n.param = "integer",
                        n.covariates = "integer",
                        lowerbounds = CLASS_CLIST,
                        upperbounds = CLASS_CLIST,
                        transform = "list",
                        #vario = "character",
                        coordunits = "character",
                        varunits = "character",
                        number.of.data = "integer",
                        modelinfo = "data.frame",
                        number.of.parameters = "integer",
                        p.proj = "integer",
                        v.proj = "integer",
                        x.proj = "ANY",  ## integer or logical
                        fixed = "ANY",
                        true.tsdim = "integer",
                        true.vdim = "integer",
                        report = "character",
                        submodels = "ANY",
                        autostart = CLASS_SINGLEFIT,
                        users.guess = CLASS_SINGLEFIT, # Martin: 2.4.: eingefuegt
                        self = CLASS_SINGLEFIT,
                        plain = CLASS_SINGLEFIT,
                        sqrt.nr = CLASS_SINGLEFIT,
                        sd.inv = CLASS_SINGLEFIT,
                        internal1 = CLASS_SINGLEFIT,
                        internal2 = CLASS_SINGLEFIT,
                        internal3 = CLASS_SINGLEFIT,
                        ml = CLASS_SINGLEFIT
                        #ml.residuals = "ANY" # matrix or RFsp
                        )
         )

## generic S4 method for 'plot'
setGeneric("plot")
setGeneric("atan2")
