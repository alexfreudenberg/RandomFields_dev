
## Authors 
## Martin Schlather, schlather@math.uni-mannheim.de
##
## Copyright (C) 2018 -- 2019 Martin Schlather
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


##.C <- function(...) { cat("\ndebugging version: "); x <- ...elt(1); print(if (is.list(x)) x[1] else x); base::.Call(C_DebugCall); z <- base::.C(...);  cat("\nend call: ");  base::.Call(C_DebugCall); z}

##.Call <- function(...) { cat("\ndebugging version: "); x <- ...elt(1); print(if (is.list(x)) x[1] else x); base::.Call(C_DebugCall); z <- base::.Call(...); cat("\nend call: "); base::.Call(C_DebugCall); z}



.onLoad <- function(lib, pkg) {
  base::.C("loadoptions")
}

GetMessage <-function() .Call("attachoptions")

.onAttach <- function (lib, pkg) {
  packageStartupMessage(GetMessage());
}

.onDetach <- function(lib) {
## do not use the following commmands in .onDetach!
}

.onUnload <- function(lib, pkg){
#   RFoptions(storing=FALSE) ## delete everything
  base::.C("detachoptions")
}
