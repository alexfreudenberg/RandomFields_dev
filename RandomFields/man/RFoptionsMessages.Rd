\name{RFoptionsMessages}
\alias{RFoptionsMessages}
\alias{notes}
\alias{hints}
\alias{warnings}
\alias{Notes}
\alias{Hints}
\alias{Warnings}
\alias{messages}
\title{Setting control arguments for hints, notes and warnings.}
\description{
 \command{\link{RFoptions}} sets and returns control arguments for the analysis
 and the simulation of random fields. It expands the functionality of
 \link[RandomFieldsUtils]{RFoptions}. Usage:\cr
 \code{RFoptions(...)}
}


\section{Details}{

Most of the options below change their value in a session
without the user's notice.

  \bold{Options for hints}
  \describe{
    \item{\code{help_addNA}}{Logical.
      If \code{TRUE} help is shown for
      \code{\link[fit]{RFoptions(addNAlintrend=)}}. 
      
	Default: \code{TRUE}.
      }

      \item{\code{help_colour_palette}}{Logical.
    If none of the packages \pkg{RColorBrewer} and \pkg{colorspace}
    are available and graphics are displayed, a message is displayed.

    Default: \code{TRUE}.
  }

  \item{\code{help_help}}{Logical.
    If \code{TRUE} help is shown for managing the appearance of messages
    from \pkg{RandomFields}.

    Default: \code{TRUE}.
  }

  \item{\code{help_mle}}{Logical.
    \command{\link{RFfit}} returns a list of fitted models according to
    set of fitting methods used. Usually the user is interested in the
    result obtained by MLE. So this result is extracted.
    If \code{note.mle=TRUE} the user is given information how to extract
    and to use other results.
  
    Default: \code{TRUE}.
  }

  
  \item{\code{help_newstyle}}{Logical.
    If \code{TRUE} a message is displayed the by the argument
    \code{spConform=FALSE} oldstyle return values are obtained instead
    of S4 objects.
    
    Default: \code{TRUE}.
  }

  \item{\code{help_normal_mode}}{Logical.
    If \code{TRUE} then the function \command{\link{RFfit}}
    displays the message that other values for the option
    \code{modus_operandi} are available.
    
     Default: \code{TRUE}.
   }

   \item{\code{help_onlyvar}}{Logical.
    If \code{TRUE} a warning is given if \command{\link{RFfit}} is used
    to fit the variance only, what should be done by
    \command{\link{RFlikelihood}}.
    
    Default: \code{TRUE}.
  }
}

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  \bold{Options for notes}
\describe{
   
  \item{\code{note_ambiguous}}{Logical.
    Usually, the argument \code{grid} in \command{\link{RFsimulate}},
    for instance, can or should be given. If not given,
    the system takes a default definition.
    Minor problem when trying to detect whether it is a grid
    are reported if \code{note_ambiguous=TRUE}.
    
    Default: \code{FALSE}.
  }

 \item{\code{note_aspect_ratio}}{Logical.
    if \code{TRUE} then a warning is given not a standard graphical
    device is used and the package plots try to keep a certain aspect
    ratio. 
    
    Default: \code{TRUE}
  }
 
  \item{\code{note_coordinates}}{Logical.
    If \code{TRUE} then a transformation from earth coordinates to
    cartesian coordinates is reported.
    
    Default: \code{TRUE}.
  }

  \item{\code{note detection}}{Logical. 
    If \code{TRUE} comments on the detection of coordinates and
    variables are given.
    
    Default: \code{TRUE}.
  }
  
  \item{\code{note_no_fit}}{Logical.
    If \code{TRUE} then the function \command{\link{RFfit}}
    displays the message if, maybe except linear function parts, no
    parameter has to estimated. Note that purely linear function parts can be
    estimated by \command{\link{RFlikelihood}}.
    
    Default: \code{TRUE}.
   }
   
 }
 \bold{Options for warnings}
 \describe{
  \item{\code{warn_ambiguous}}{Logical.
    Usually, the argument \code{grid} in \command{\link{RFsimulate}},
    for instance, can or should be given. If not given,
    the system takes a default definition.
    More severe problem when trying to detect whether it is a grid
    are reported if \code{warn_ambiguous=TRUE}.
    
    Default: \code{TRUE}.
  }
        
  \item{\code{warn_constant}}{Logical.
    The definition of \command{\link{RMconstant}} has changed.
    A warning is displayed if the command is used. \code{warn_constant}
    will become obsolete in future versions.
    
    Default: \code{TRUE}.
  }
  
   \item{\code{warn_coord_change}}{Logical.
    If \code{TRUE} then it is reported if the covariance model relies on
    to different coordinate systems.
    
    Default: \code{TRUE}.
  }
  
 % \item{\code{warn_mathdef}}{Logical. PROPABLY obsolete}
   
  \item{\code{warn_missing_zenit}}{Logical. % to do
    Only for Earth systems: a missing zenit is frequently a cause
    for errors that are difficult to understand. Therefore, in such
    cases an additional warning message is displayed.
    
    Default: \code{TRUE}
  }

  % warn_mode internal
  % warn_modus_operandi  internal

  
   \item{\code{warn_negvar}}{Logical. In very special constructions of
    the argument \code{var} the non-negativiy of the variance cannot be
    detected, e.g. when the variance is given spatially varying function.
    If \code{warn_negvar=TRUE} it is reported when the determination of
    non-negativity is not possible.
     
    Default: \code{TRUE}.
  }

  \item{\code{warn_newAniso}}{ obsolete.\cr
    Logical.
    If \code{newAniso=TRUE} and the argument \code{Aniso} is used in the model
    definition, then a message is displayed that the matrix \code{Aniso}
    is multiplied from the right by \eqn{x}, where up to Version 2.0
    the argument \code{aniso} was available which was multiplied from
    the left by \eqn{x}.
     
    Default: \code{TRUE}.
  }

   
%  \item{\code{warn_oldstyle}}{Logical.
%    If \code{TRUE} a warning is given if an obsolete function
%    from Version 2 is used.\cr   Default: \code{TRUE}.  }

  \item{\code{warn_on_grid}}{Logical.
    If a (one-dimensional) grid is given, but the argument
    \code{grid=FALSE}, e.g. in \code{RFsimulate}, this contraction is
    reported if \code{warn_on_grid=TRUE}
    
    
    Default: \code{TRUE}.
  }
  
  \item{\code{warn_parallel}}{Logical.
    The package \pkg{RandomFields} should not be called from within
    functions of the packages \pkg{parallel} and alike.
    If \code{TRUE} and \pkg{RandomFields} parallel calls a warning is
    given.     
    
    Default: \code{TRUE}.
  }

  
  \item{\code{warn_scale}}{Logical.
    If \code{warn_scale=TRUE} then a scale less than 10 [km] is reported
    if earth coordinates are transformed to cartesian coordinates.
    
    Default: \code{TRUE}.
  }

  \item{\code{warn_seed}}{integer.
    If \code{warn_scale}\eqn{> 0} a note is given when the seed is fixed.    
    Default: \code{2}.
  }

  \item{\code{warn_singlevariab}}{logical.
    If a covariance model is not stationary and pairs of
    coordinates are not given, the second part of the pair is set to
    zero. If \code{warn_singlevariab=TRUE}, additionally a warning is
    given. 
   
    Default: \code{TRUE}.
  }

}
  %%%%%%%%%%%
  \bold{Warning options taken from the package \pkg{RandomFieldsUtils}}
  \describe{
    \item{\code{warn_unknown_option}}{integer. This options is
      best set by a separate call, i.e. by
      \code{RFoptions(warn_unkown_option= )}. Values
      \itemize{
	\item{0}{no warning at all}
	\item{+1,-1}{unkown option with lower case initial leads to an
	  error}
	\item{+2,-2}{unknown option that does not consist of a single
	  capital letter leads to an error}
	\item{+3}{any unknown options leads to an error}
      }
      Positive values lead to a message about the ignored unknown option.

      Further values are used for internal purposes and should not be used.   
      
      Default: \code{3}.
    }
  }
}

\seealso{
  \command{\link{RFoptions}}
  }

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

RFoptions()$messages
\dontshow{FinalizeExample()}}


\keyword{spatial}

