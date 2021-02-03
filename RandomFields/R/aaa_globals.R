## Authors 
## Martin Schlather, schlather@math.uni-mannheim.de
## Sebastian Gross
##
##
## Copyright (C) 2017 - 2018 Martin Schlather
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

SYMBOLS <- c(SYMBOL_PLUS, SYMBOL_MULT , RM_S[2])
iR_P <- paste0("i", R_P)
RC_MAXSTABLE_NAMES <- c(BR_NAME, OPITZ_NAME, EG_NAME, SMITH_NAME)


CLASS_CLIST <- 'RMmodel'
CLASS_RM <- 'RMmodelgenerator'
CLASS_SINGLEFIT <- 'RMmodelFit'
CLASS_FITLIST <- 'RFfit'
CLASS_EMPIR <- "RFempVariog"
CLASS_PLOT <- "RFplot"
 
VOID <- "..void.."

## if names are changed fitgauss.R, weights, must be changed
LSQMETHODS <- c("plain", "self", "sqrt.nr", "sd.inv",
                "internal1", "internal2", "internal3") 
MLMETHODS <- c("ml") # "reml", "rml1"),
PRIMMETHODS <- c("users.guess", "autostart")
CROSSMETHODS <- NULL
METHOD_PREFLIST <- c(MLMETHODS, LSQMETHODS, PRIMMETHODS)

par.storage <- ".RandomFields.par"
.RandomFields.env <- new.env()

## ACHTUNG! Die in C definierten PL_* haben andere Bedeutung
PL_IMPORTANT 	<- as.integer(1)
PL_SUBIMPORTANT 	<- as.integer(2)
PL_DETAILSUSER <- as.integer(3)## currently unused
PL_RECURSIVE 	<- as.integer(4)
PL_STRUCTURE 	<- as.integer(5)
PL_ERRORS 	<- as.integer(6)

PL_FCTN_DETAILS 	<- as.integer(7)
PL_FCTN_SUBDETAILS 	<- as.integer(8)

PL_DETAILS 	<- as.integer(9)
PL_SUBDETAILS 	<- as.integer(10)



isPosDef <- function(type) {
  if (is.character(type)) type <- pmatch(type, TYPE_NAMES, duplicates.ok=TRUE)-1
  ##  .C(C_isPosDef, as.integer(type))$type
  type==TcfType | type == PosDefType | type == ManifoldType
}
isVariogram <- function(type) { 
  if (is.character(type)) type <- pmatch(type, TYPE_NAMES, duplicates.ok=TRUE)-1
  ##  .C(C_isNefDef, as.integer(type))$type
  isPosDef(type) | type == VariogramType
}
