\name{RFoptions}
\alias{RFoptions}
\alias{RFOPTIONS}
\alias{RC_DUPLICATEDLOC_NAMES}
\title{Setting control arguments}
\description{
 \command{\link{RFoptions}} sets and returns control arguments for the analysis
 and the simulation of random fields. It expands the functionality of
 \link[RandomFieldsUtils]{RFoptions}. Usage:\cr
 \code{RFoptions(...)}
}
%\section{Usage}{
% RFoptions(...)
%
%section{Arguments}{
%  \describe{  
%    \item{...}{arguments in \code{tag = value} form, or a list of tagged
%      values.}
%  }
 %\item{no.readonly}{If \command{\link{RFoptions}} is called without
%   argument then all arguments are returned in a list. If
%   \code{no.readonly=TRUE} then only rewritable arguments are returned.
% }
%}
\section{Details}{
  The subsections below comment on\cr
  \bold{0. \code{basic}: \RFU}\cr
  \bold{1. \code{general}: General options}\cr
  \bold{2. \code{br}: Options for Brown-Resnick
    Fields}\cr
  \bold{3. \code{circulant}: Options for circulant embedding methods
    \command{\link{RPcirculant}}}\cr
  \bold{4. \code{coords}: Options for coordinates and units, see
    \link{coordinate systems}}\cr
  \bold{5. \code{direct}: Options for simulating by simple matrix decomposition}\cr
  \bold{6. \code{distr}: Options for distributions, in particular \command{\link{RRrectangular}}}\cr
  \bold{7. \code{empvario}: Options for calculating the empirical variogram}\cr
  \bold{8. \code{\link{fit}}: Options for \command{\link{RFfit}},
    \command{\link{RFratiotest}}, and \command{\link{RFcrossvalidate}}}\cr
  \bold{9. \code{gauss}: Options for simulating Gaussian random fields}\cr
  \bold{10. \code{graphics}: Options for graphical output}\cr
  \bold{11. \code{gui}: Options for \command{\link{RFgui}}}\cr
  \bold{12. \code{hyper}: Options for simulating hyperplane tessellations}\cr
  \bold{13. \code{krige}: Options for Kriging}\cr
  \bold{14. \code{maxstable}: Options for simulating max-stable random fields}\cr
  \bold{15. \code{mpp}: Options for the random coins (shot noise) methods}\cr
  \bold{16. \code{nugget}: Options for the nugget effect}\cr
  \bold{17. \code{registers}: Register numbers}\cr
  \bold{18. \code{sequ}: Options for the sequential method}\cr
  \bold{19. \code{solve}: Options for solving linear systems}\cr
  \bold{20. \code{special}: Options for some special methods}\cr
  \bold{21. \code{spectral}: Options for the spectral (turning bands) method}\cr
  \bold{22. \code{tbm}: Options for the turning bands method}\cr
  \bold{23. \code{\link{message}}: Options for hints, notes and
    warnings}\cr
 \bold{24. \code{internal}: Internal options, which should not be
   changed by the user}\cr
  \cr\cr
  
  

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% \bold{16. Options for RFloglikelihood}\cr
% 
% "auto", "full", "composite", "selection"




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bold{1. General options}
\describe{
  % only works for RFsimulate, not for RFsimulate
  
  \item{\code{duplicated_locations}}{Allowed values are
    given by \code{RC_DUPLICATEDLOC_NAMES}:
%    \code{"error"},   \code{"risk error"}, \code{"regard as repeated"} and \code{"scatter"}:
    \describe{
      \item{\code{"error"}}{
	In case of kriging, locations that appear in both the given data
	and the locations to be predicted, are searched and eliminated.
	This can be very time consuming.

	In other situations (nearly)
	duplicated locations will mostly end up in an
	error, e.g. as covariance matrices cannot be inverted, or by
	explicit detection (in some cases).
      }

      \item{\code{"risk error"}}{
	No strategy is applied to avoid errors.
	But also no early error checks are performed.
	Duplicated locations will mostly end up in an
	error.
      }

      \item{\code{"regard as repeated"}}{ duplicated locations are
	considered 
	as repeated measurement. The (standard) nugget effect becomes a
	non-stationary 
	model in an abstract space that cannot be extended outside the given
	locations. See also \link{RPnugget} for the distinction between
	measurement error (leading always to a non-stationary model!)
	and spatial nugget.
	When kriging is performed, an error
	model must be given and kriging is not an exact interpolator
	anymore.
      }

      \item{\code{"scatter"}}{
	Coordinates that are not on a grid are slightly scattered, so that
	the points have some tiny distances.

	The exact behaviour behaves heavily on
	\code{RFoption(nugget.tol=)}.
	If the latter is zero, scattering will lead to different locations
	that have different (spatial) nugget value.
	If \code{nugget.tol}\eqn{>0}, then scattering will leave
	identical locations within a neighbourhood less than
	\code{nugget.tol}, so that they have still the same nugget
	value.
      }
    }
    
    Default: \code{"risk error"}.
  }

  \item{\code{cPrintlevel}}{
    \code{cPrintlevel} is automatically set to \code{printlevel}
    when \code{printlevel} is changed.
    Standard users will never use a value higher than 3.
    
    0 : no messages\cr
    1 : messages and warnings when the user's input looks odd\cr
    2 : messages (and internal errors) documenting the choice of the
    simulation method\cr
    3 : further user relevant informations\cr
    4 : information on recursive function calls\cr
    5 : function flow information of central functions \cr
    6 : errors that are internally treated\cr
    7 : details on building up the covariance structure\cr
    8 : details on taking the square root of the covariance matrix\cr
    9 : details on intermediate calculations\cr
    10 : further details on intermediate calculations\cr
    
    Note that \code{printlevel} works
    on the R level whereas \code{cPrintlevel} works on the C level.
    
    Default: 1 \cr % [also  do].\cr
  }
  
  \item{\code{detailed_output}}{logical.
    if \code{TRUE} some function, e.g. \code{\link{RFcrossvalidate}}
    will return additional information.
  }
  
  \item{\code{every}}{integer.
    if greater than zero, then every \code{every}th iteration is
    printed if simulated by TBM or random coin method. The value zero
    means that nothing is printed.
    
    Default: \code{0} % [do].
  }
  
  
  \item{\code{exactness}}{logical or NA. Currently only used when simulating
    Gaussian random fields.
    \itemize{
      \item \code{TRUE}: \command{\link{RPcoins}},
      \command{\link{RPhyperplane}}, \command{\link{RPsequential}},
      \command{\link{RPspectral}} and \command{\link{RPtbm}} and
      \emph{approximative} circulant embedding  are excluded.
      If the circulant embedding method is considered as badly
      behaved, then the matrix decomposition methods are preferred.
      
      \item \code{FALSE}: all the methods are allowed.
      If the circulant embedding method is
      considered as badly behaved or the number of points to be
      simulated is large, the turning bands methods are
      rather preferred.
      
      \item \code{NA}: Similar to \code{FALSE}, but
      some inexact algorithms get less preference.
    }
    
    Default: \code{NA} .
  }
  
  \item{\code{expected_number_simu}}{positive integer which is usally set
    internally as the value of the argument \code{n} in
    \command{\link{RFsimulate}}. The argument \code{expected_number_simu}
    should be set only by an advanced users and only if 
    \command{\link{RFsimulate}} will be called with argument \code{n} alone.
  }
  
  
  \item{\code{gridtolerance}}{
    used in \command{\link{RFsimulate}} to see if the coordinates build a
    grid for x, y, z, T-values. This argument is also used
    in case of conditional
    simulation where the data locations might ly on
    a grid.
    
    Default: \code{1e-6}
  }

  \item{asList}{logical. Lists of arguments are treated slightly
    different from non-lists. If \code{asList=FALSE} they are treated the
    same way as non-lists. This options being set to \code{FALSE} after
    calling \command{RFoptions} it should be set as first element of a list.
    
    Default: \code{TRUE}
  }
  
  \item{\code{modus_operandi}}{character. One of the values
    \code{"careless"}, \code{"sloppy"}, \code{"easygoing"},
    \code{"normal"}, \code{"precise"}, \code{"pedantic"},
    \code{"neurotic"} .
    \bold{This argument is in an experimental stage and its definition
      and effects will change very likely in near future.}%to do
    This argument sets a lot of argument at once related to estimation
    and simulation. \code{"careless"} prefers rather fast algorithms,
    but the results
    might be very rough approximations. By way of contrast,
    \code{"neurotic"} will try very
    hard to return exact result at the cost of hugh computing times.
    
    Default: \code{"normal"}
  }
  
  
  \item{\code{na_rm_lines}}{
    logical. If \code{TRUE} then a line of the data that contains a
    \code{NA} value is deleted. Otherwise it is tried to deal with the
    \code{NA} value at higher costs of computing time. (Only used for
    kriging -- estimation can fully deal with \code{NA}s.)
    
    Default: \code{FALSE}.
  }
  
  \item{\code{output}}{character.
    one of the values \code{"sp"} (if and only if
    \code{spConform=TRUE}), 
    \code{"RandomFields"} (if and only if \code{spConform=FALSE}),
    \code{"geoR"}.

    The output mode \code{geoR} currently adds some attributes such as
    the call of the function.

    NOTE: \code{output} is in an experimental stage, whose effects might
    change in future. Currently, \code{output} changes the values of
    \code{reportcoord}, \code{returncall} and \code{spConform}.
  }
  
  \item{\code{pch}}{character.
    \command{\link{RFfit}}: shown before evaluating any method;
    if \code{pch!=""} then one or two
    additional steps in the MLE methods are
    marked by \dQuote{+} and \dQuote{#}.
    
    Simulation: 
    
    The character is printed after each
    performed simulation if more than one simulation is performed at
    once. If \code{pch='!'} then an absolute
    counter is shown instead of the character.
    If \code{pch='\%'} then a
    counter of percentages is shown instead of the character.
    Note that also \sQuote{\eqn{\mbox{\textasciicircum}}{^}H}s are printed in
    the last two cases, 
    which may have undesirable interactions with some few other R
    functions, e.g. \command{\link[utils]{Sweave}}.
    
    Default: \code{'*'}. % [do]. 
  }
 
 
 

  \item{\code{practicalrange}}{logical or integer.
    If not \code{FALSE} the range of primitive
    covariance functions is
    adjusted so that cov(1) is zero for models with finite range.
    (Operators are too complex to be adjusted; for anisotropic
    covariance the practical range is not well defined.)
    
    The value of cov(1) is about 0.05 (for \code{scale=1})
    for models without range. See \command{\link{RMmodel}} or type
    \cr
    \code{\link{RFgetModelNames}(type="positive definite",
      domain="single variable", isotropy="isotropic", operator=FALSE, vdim=1)}
    \cr
    for the list of primitive models.
    \itemize{
      \item \code{FALSE} : the practical range ajustment is not used.
      \item \code{TRUE} : \code{practicalrange} is applicable only if
      the value is known exactly, or, at least, can be approximated by
      a closed formula.
      \item \code{2} : if the practical range is not known exactly it
      is approximated numerically.	
    }
    
    Default: \code{FALSE} .
  }
  
  \item{\code{printlevel}}{If \code{printlevel}\eqn{\le0}{<=0}
    there is not any output on the screen. The
    higher the number the more tracing information is given. 
    Standard users will never use a value higher than 3.

    0 : no messages\cr
    1 : important (error) messages and warnings\cr
    2 : less important messages\cr
    3 : details, but still for the user\cr
    4 : recursive call tracing (only used within \command{\link{RFfit}})\cr
    5 : function flow information of large functions\cr
    6 : errors that are internally treated\cr
    7 : details on intermediate calculations\cr
    8 : further details on intermediate calculations\cr
    
    Default: 1 %[also do].\cr
  }

  \item{reportcoord}{character.
    Current values are \code{"always"}, \code{"important"},
    \code{"warn"}, \code{"never"},

    Both \code{"warn"} and \code{"important"} have any effect only
    if the coordinate system is changed internally. In this case
    \code{"warn"} yields a displayed warning message whereas
    \code{"important"} adds an attribute to the result as in the
    case \code{"always"}.
    
    If \code{"always"} or \code{"important"}
    the reports are added as attribute to the results.
    Note that in this case the class of the result may change
    (e.g. from \code{"numeric"} to \code{"atomic"}).

    Default: \code{"warn"}
  }

  \item{returncall}{logical. If \code{TRUE}
    then the call is returned as an attribute

    Default: \code{TRUE}
  }
   
  \item{seed}{integer. If \code{NULL} or \code{NA}
    \command{\link{set.seed}} is \bold{not} called.
    Otherwise, \code{\link[base]{set.seed}(seed)} is set
    before simulations are performed, e.g. by
    \command{\link{RFsimulate}} or \command{\link{RFdistr}}.
    
    If the argument is set locally, i.e., within a function,
    it has the usual local effect. If it is set globally, i.e. by
    \command{RFoptions} the \code{seed} is fixed
    for \bold{all subsequent} calls.
    
    If the number of simulations \code{n} is greater than one
    and if \code{RFoptions(seed=seed)} is set, the \eqn{i}th
    simulation is started with the seed \sQuote{\code{seed}\eqn{+i-1}}.
    %   The function \code{set.seed} should not be used in case \code{n}
    %    is greater than 1.
    %
    %Vgle!   
    %set.seed(5)
    %RFsimulate(RPschlather(RMmatern(nu=2), xi=1, mu=1, s=1), x, grid=F, n=5)@data
    %set.seed(5)
    %RFsimulate(RPschlather(RMmatern(nu=2.01), xi=1, mu=1, s=1), x,grid=F,n=5)@data
    %RFoptions(cPr=3, seed=5)
    %RFsimulate(RPschlather(RMmatern(nu=2), xi=1, mu=1, s=1), x, grid=F, n=5)@data
    %RFsimulate(RPschlather(RMmatern(nu=2.01), xi=1, mu=1, s=1), x, grid=F,n=5)@data

    
    Note also that \command{\link{RFratiotest}} has its own argument
    \code{seed} with a slightly different meaning.
  }

   
  \item{seed_incr, seed_sub_incr}{ (does not work yet)
    This argument is important iff \command{RFsimulate} is used within
    a function from package \pkg{parallel}.
    The value of \code{seed_incr} should be set only locally, i.e. not
    by \code{RFoptions()}.
 
    If \code{seed_incr != 0} (or the number of simulations \code{n} is
    greater than 1) and \code{!is.na(seed)}
    then the seed for each simulation is calculated as
    
    \code{seed} \eqn{+ (k-1) *} \code{seed_sub_incr}  \eqn{+} \code{seed_incr} \eqn{* n}
    
    where \eqn{k} runs from 1 to \code{n}.

    
    Default: 0
  }

  \item{\code{set}}{integer.
    Certain models (e.g. \command{\link{RMfixcov}} and
    \command{\link{RMcovariate}}) 
      allow for lists as arguments.
    \code{set} selects a certain list element.
    If necessary the list is recycled.
  }

  \item{\code{spConform}}{logical.
    \code{spConform=TRUE} might be used by
    a standard user as this allows the comfortable use of \command{plot},
    for instance, while \code{spConform=FALSE} is \bold{much} faster and
    and consumes \bold{much less memory}, hence might
    be used by programmers or advanced users.
    
    Details: if \code{spConform=TRUE} then \command{\link{RFsimulate}} and
    many other functions 
    return an \code{sp}-object (which is an S4 object). Otherwise, matrices
    or lists are 
    returned as defined in RandomFields 2.0, see the manuals for the
    specific functions. Frequently, the latter have now a class attribute
    to make the output nicer.
    
    Note: for large data sets (to be generated),
    \code{spConform=TRUE} should \bold{not} be used.

    See also \code{output}.
    
    Default: \code{TRUE} %[do].
  }
  

  \item{\code{skipchecks}}{logical.
    If \code{TRUE}, several checks whether the given parameter values
    and the dimension are within the allowed range is skipped.
    Do not change the value of this variable except you really
    know what you do.
    
    Default: \code{FALSE} $ %[also do]. 
  }
  
  \item{\code{storing}}{Logical.
    If \code{FALSE} then the intermediate results are
    destroyed after the call of a function \code{\link[=RF]{RFfoo}}.
    
    If \code{storing=TRUE}, in particular
    additional simulations can be performed by calling
    \command{\link{RFsimulate}} with at most the argument \code{n}.
    This call can then be much faster, but a rather large
    amount of memory will be kept.
    
    When \code{storing} is turned from \code{TRUE} to \code{FALSE} by
    global call then all is deleted.
    
    Default: \code{FALSE} 
  }
  
  \item{\code{Ttriple}}{Logical or \code{NA}.
    If \code{TRUE}, then triple for the time argument \code{T} is
    expected, containing start, step (by), length.
    If \code{FALSE} a sequence on a grid is expected.
    If \code{NA} then the decision is automatic, but will lead to an
    error if ambiguous.
  }
  
  \item{\code{vdim_close_together}}{logical. Used especially in functions that
    create covariance matrices. If the model is multivariate, then two
    ways of ordering the matrix exist. To consider first all variables at
    a certain location (\code{vdim_close_together=TRUE}) or to consider first
    all locations keeping the variable fixed
    (\code{vdim_close_together=FALSE}).
    Note that several simulation methods rely on the value \code{FALSE},
    so that these methods will not work anymore if
    \code{vdim_close_together=FALSE}. 
    
    Default: \code{FALSE}.
  }
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bold{2. Options for Brown-Resnick Fields}
\describe{

  \item{\code{deltaAM}}{ 
    integer; only used for simulation of BR processes via \code{RPbrmixed}
    with \code{optim_mixed=2}. In this case, \code{deltaAM} is the number 
    of additionally simulated Gaussian processes used for an update of
     \code{\link[=RPbrmixed]{areamat}} in the optimization pricedure.
     
     Default: \code{300}
  }
  
  \item{\code{maxtrendmem}}{
    integer; the maximal number of real valued variables used for
    intermediate storage:
    \itemize{
      \item
      RPbrshifted: trends for shifted locations that may
      be stored at the same time when simulating BR processes.
      \item
      RPbrnormed: Let \eqn{n} be the number of locations.
      Then a \eqn{n\times n} (covariance) matrix has to be evaluated
      at random columns. 
    }
    if \code{maxtrendmem} is large (and \eqn{n} small, \eqn{n\le 10^4}),
    multiple evaluations can be avoided. 
    
    Default: \code{1e7} .
  }
  \item{\code{meshsize}}{
    positive; width of the grid on which the shape functions in the M3
    representation of BR processes are simulated; only used for
    simulation of BR processes via \code{RPbrmixed}.
    
    Default:  \code{0.1} .
  }

  \item{\code{optim_mixed}}{\code{0, 1, 2}; only used for simulation of BR
    processes via \code{RPbrmixed}.\cr
    If \code{optim_mixed=0}, the arguments
    \code{\link[=RPbrmixed]{lambda}} and 
    \code{\link[=RPbrmixed]{areamat}} of \command{\link{RPbrmixed}}
    are used for the simulation.\cr
    If \code{optim_mixed=1}, \code{\link[=RPbrmixed]{lambda}} is estimated for
    \code{\link[=RPbrmixed]{areamat=1}}.\cr
    If \code{optim_mixed=2}, \code{\link[=RPbrmixed]{areamat}} is optimized and
    \code{\link[=RPbrmixed]{lambda}} is estimated.
    
    Default: \code{1} .
  }
  \item{\code{optim_mixed_tol}}{
    value in \eqn{[0,1]}; only used for simulation of BR processes via
    \code{RPbrmixed} with \code{optim_mixed=2}. In this case,
    \code{\link[=RPbrmixed]{areamat}} is optimized under the constraint that the
    probability of drawing the shape function incorrectly is bounded by
    \code{optim_mixed_tol} (cf. Oesting et al., 2012).
    
    Default: \code{0.01} .
  }
  
  \item{\code{variobound}}{
    positive; the shape functions in the mixed moving maxima
    representation are cut off where the variogram belonging
    to \code{phi} exceeds \code{variobound}.
    
    Default: \code{8.0} .
  }
  
  \item{\code{vertnumber}}{
    positive integer; for an efficient simulation of the shape functions
    in the M3 representation of BR processes, the component \eqn{E} from
    of the domain \eqn{[x_0, \infty] \times E}{[x_0, Inf] x E} of the
    underlying Poisson point process is sub-dividedinto cubes
    (cf. Oesting et al., 2012); \code{vertical}  is the number of
    vertical breaks of \eqn{E}; only used for simulation of BR processes
    via \code{RPbrmixed} with \code{optim_mixed=2}.
    
    Default: \code{7} .
  }
   
}


 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 \bold{3. \code{circulant}: Options for circulant embedding methods, cf. \command{\link{RPcirculant}}}\cr
   These options influence the standard circulant embedding
   method, cutoff circulant embedding  intrinsic circulant embedding.
   It can also influence \command{\link{RPtbm}} if the line is simulated
   with any circulant embedding method.
   \cr
 \describe{
 \item{\code{approx_maxgrid}}{See \command{\link{RPcirculant}}}
 \item{\code{approx_step}}{See \command{\link{RPcirculant}}}
 \item{\code{dependent}}{See \command{\link{RPcirculant}}}
 \item{\code{force}}{See \command{\link{RPcirculant}}}
 \item{\code{maxGB}}{See \command{\link{RPcirculant}}}
 \item{\code{maxmem}}{See \command{\link{RPcirculant}}}
 \item{\code{mmin}}{See \command{\link{RPcirculant}}}
 \item{\code{strategy}}{See \command{\link{RPcirculant}}}
 \item{\code{tolIm}}{See \command{\link{RPcirculant}}}
 \item{\code{tolRe}}{See \command{\link{RPcirculant}}}
 \item{\code{trials}}{See \command{\link{RPcirculant}}}
 \item{\code{useprimes}}{See \command{\link{RPcirculant}}}
 }
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bold{4. \code{coords}: Options for coordinates and units}
\describe{
  \item{\code{coord_system}}{character. See \link{coordinate systems}
  }
  
  \item{\code{coordunits}}{ See \link{coordinate systems}
  }
 
  \item{\code{coordnames}}{See \link{coordinate systems}
  }

  
  \item{\code{new_coord_system}}{See \link{coordinate systems}
  }

  \item{\code{new_coordunits}}{See \link{coordinate systems}
  }

  \item{\code{polar_coord} }{See \link{coordinate systems}
  }

  \item{\code{varnames}}{See \link{coordinate systems}
  }  

  \item{\code{varunits}}{See \link{coordinate systems}
  }
    
  \item{\code{xyz_notation}}{See \link{coordinate systems}
  }
 
  \item{\code{zenit}}{See \link{coordinate systems}
  }
  
}





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bold{5. \code{direct}: Options for simulating by simple matrix decomposition}
\describe{
  \item{\code{max_variab}}{Maximal size of the covariance matrix.

  Default: 12000
  }
}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bold{6. \code{distr}: Options for distributions, in particular \command{\link{RRrectangular}}}
\describe{
  \item{\code{innermin}}{
    Default value to simulate from the
    \command{\link{RRrectangular}} distribution.
    The minimal length of the interval where the Taylor expansion shall
    be valid.
   
    Default: \code{1e-20} .
  }

  \item{\code{maxit}}{
    Default value to simulate from the
    \command{\link{RRrectangular}} distribution.

    The number of iterative steps where the
    the constant of the Taylor development is increased,
    to find an upper bound for the given function.
    
    Default: \code{20} .
  }

  \item{\code{maxsteps}}{
    Default value to simulate from the
    \command{\link{RRrectangular}} distribution.

    \code{maxsteps} is usually the number of steps in the middle part of
    the approximation. From this value and the length between
    the determined endpoints for  the approximation at the origin and in
    the tail, the step length is calculated. If the step length is less
    than \code{minsteplen} the number of steps is reduced.
    
    Default: \code{1000} .
  }

  
  \item{\code{mcmc_n}}{
    In case of the use of MCMC it leaves out \eqn{n-1}
    member of the Markov chain bevor the \eqn{n} member
    is returned. See also maxsteps.
    
    Default: \code{15} .
  }

  
  \item{\code{minsteplen}}{
    Default value to simulate from the
    \command{\link{RRrectangular}} distribution.
    The minimal step length
    for the middle part of approximation, which is a step function,
    
    Default: \code{0} (i.e. not used as a criterion.)
  }
  
  \item{\code{outermax}}{
    Default value to simulate from the
    \command{\link{RRrectangular}} distribution.
    The largest possible endpoint for the middle part that
    approximates the function by a step function. See also \code{innermax}.

    Default: 20.
  }
 
  \item{\code{parts}}{
    Default value to simulate from the
    \command{\link{RRrectangular}} distribution.

    \code{parts} determines the number of tests that are performed to
    check whether a proposed power function is an upper bound for
    the given function, at the origin and the tail.
    
    Default: \code{8} .
  }

  \item{\code{repetitions}}{
    Minimal number of realisations to determine a quantity of the
    distribution by MCMC. E.g. to determine the integral value \eqn{c}
    in the paper of Oesting, Schlather, Zhou.
    
    Default: 1000.
  }
  
  \item{\code{safety}}{ 
     Default value to simulate from the
     \command{\link{RRrectangular}} distribution.

     First, at the origin, the first power function of the Taylor
     expansion is taken as potential upper function.
     The constant of the power function are increased by factor
     \eqn{1 + }\code{safety} and the exponent of the function
     similarly decreased. A number of test evaluations
     is performed to check whether this modified function is indeed
     a upper bound. If not, the considered interval at the origin
     is reduced iteratively, the constants of the power function
     further increased and the exponent decreased.
     If \code{maxit} iteration have been performed without success,
     the search for an upper bound fails.
     The search at the origin also fails if the interval around
     the origin has become less than \code{innermin}.
     

     Similar procedure is performed for the tail.
     
   
    Default: \code{0.08} .
  }
}


 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bold{7. \code{empvario}: Options for calculating the empirical variogram}
\describe{
   \item{\code{fft}}{
    Logical. Determines whether FFT should be used for data on a grid

    Default: \code{TRUE}.
  }
  \item{\code{phi0}}{
    numeric. In case of anisotropic fields directional cones are
    considered. The argument \code{phi0} determines the starting angle.
    
    Default: \code{0}.
  }
%  \item{\code{pseudovariogram}}{
%    logical. Only in the multivariate case. Whether the
%    pseudovariogram or the crossvariogram should be calculated.
%    
%    Default: \code{FALSE}.
%  }
  \item{\code{theta0}}{
    numeric. In case of anisotropic fields directional cones are
    considered. The argument \code{theta0} determines one of the
    boundaries, hence all boundaries for a given fixed number of cones.
    The argument \code{theta0} determines the starting value of the
    second anglue in polar coordinate representation in 3 dimensions.
     
    Default: \code{0}.
  }
  \item{\code{tol0}}{
    numeric. Estimated values of the empirical variogram
    below \code{tol0} times the grid step in the third dimension
    are considered to be zero. Hence the respective values are set
    to zero.
    
    Default: \code{1e-13}.
  }
}
 
 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 \bold{8. \code{\link{fit}}: Options for \command{\link{RFfit}},
   \command{\link{RFratiotest}}, and \command{\link{RFcrossvalidate}}}

 see \code{\link{fit}}.

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bold{9. \code{gauss}: Options for simulating Gaussian random fields}
\describe{
  \item{\code{approx_zero}}{
    Value below which a correlation is considered to be essentially zero.
    This argument is used to determine the practical range of covariance
    function with non-compact support.
    
    Default: \code{0.05}
  }
 
 \item{\code{boxcox}}{
   real vector of one or two components.
   If the first component is \code{Inf} then no transformation
   is performed. Otherwise the BoxCox transformation is performed.
   Note that Box Cox only works in a Gaussian framework.
   Note further that either \code{boxcox} or \code{loggauss}
   may be given.
  
   Default \code{c(Inf, 0)}
  }


  \item{\code{direct_bestvar}}{integer.
    When searching for an appropriate simuation method
    the matrix decomposition method (\code{method="direct"})
    is preferred if the number of variables is less than or equal to 
    \code{direct_bestvariables}.
    
    Default is \code{1200}.
  }

  \item{\code{loggauss}}{
    logical. Whether a log-Gauss random fields should be returned.
    See also \code{boxcox} for a generalisation.
  }

  \item{\code{paired}}{
    (\dQuote{Antithetic pairs}.)
    Logical. If \code{TRUE} then the second half of the
    simulations is logical. If \code{TRUE} then the second half of the
    simulations is obtained by
    only changing the signs of all the standard Gaussian random variables,
    on which the first half of the 
    simulations is based. Default is \code{FALSE}.
  }
  
  \item{\code{stationary_only}}{
    See \command{\link{RPgauss}}
  }
 
}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bold{10. \code{graphics}: Options for graphical output}
These parameters should be set always globally by \command{RFoptions}.
Options for \command{\link[graphics]{par}} are passed through
\command{\link{RFpar}}. 
\describe{
  \item{\code{always_close_device}}{logical.
    If \code{FALSE} the current device is kept as it is;
    otherwise the current device is closed before the next
    device is opened. If \code{NA} it closes the preceding device
    if the opened device is pdf or jpeg.
      
    Default: \code{NA}.
  }
  \item{\code{always_open_device}}{logical.
    If \code{TRUE} a new graphical window is opened for every
    \command{\link[graphics]{plot}} if a standard graphical output is
    used, trying to respect the aspect ratios for the plots.
    The devices pdf and jpeg are always opened.

    If \code{NA} then the value is set
    to \code{\link[base]{interactive}()}. 
    
    Default: \code{TRUE}.
  }

  \item{\code{close_screen}}{logical; only relevant if
    \code{split_screen = TRUE} and \code{always_close_screen = FALSE}.
    If \code{FALSE} the windows opened by
    \command{\link[graphics]{split.screen}} are left open.
    
    Default: \code{TRUE}.
  }

  \code{file}{character; only relevant if
    \code{split_screen = TRUE}.
    argument \code{file} in \command{\link{pdf}}
    If \code{""} then no internal naming is performed.

    Default: \code{""}.
  }

  \item{\code{filenumber}}{integer; only relevant if
    \code{split_screen = TRUE}. Starting number of the file if
    \code{onefile=FALSE}. It is set to 0 whenever \code{file} is
    changed and \code{onefile=FALSE}. 

    Default 0.
  }

  \item{\code{grDefault}}{
    logical. If \code{FALSE} the graphic style up to Version 3.2
    is used. Otherwise, the changes of th graphical style are reduced to
    a minimum.
    
    Default: \code{FALSE} 
  }

  \item{\code{grPrintlevel}}{
    integer values 0, 1, 2; only relevant when simulations are
    plotted. The higher the more text is shown in the 
    plot.
    
    Default: \code{0}.
  }
  
  \item{\code{height}}{real number; only relevant if
    a new device is opened, see \code{alwyas_open_screen}.
    \itemize{
      \item \code{height=NA} or \code{height} is not positive: no device
      is opened.
      \item \code{width = NA}
      If \code{height} is greater than zero then it gives the height
      of a single figure in a plot created by \pkg{RandomFields};
      See also \code{close_screen}.
      
      If plots with multiple figures are
      shown, the height and width of the plot
      will be increased by a factor up the
      ones given by \code{increase_upto}.

      The width is calculated so that the aspect ratio is correct.
     \item \code{width} not \code{NA}
      \code{height} and \code{width} give the size of the
      whole window.
    }
    Default: \code{6}.
  }
  
   
  \item{\code{increase_upto}}{
    See \code{height}.
    
    Default: \code{c(3,4)}.
  }

   \item{\code{maxchar}}{positive integer; (rough) maximum length of
   plot title for model plots.
       
    Default: \code{20}.
  }

  \item{\code{n.points}}{vector of one or two positive integers;
    resolution of plot in each dimension.
    If two values are given the second one is used for one-dimensional
    plots and contour plots (which are usually much denser than image
    plots)
    
    Default: \code{c(100,200)}.
  }

  \item{\code{split_screen}}{logical.
     If \code{TRUE} \command{\link[graphics]{split.screen}}
     is used to split the screen.
     Otherwise \code{par(mfcol)}.
     When using \code{split_screen} then the figures tend to be fancier.
     
    Default: \code{TRUE}.
  }

  \item{\code{onefile}}{logical; only relevant if
    \code{split_screen = TRUE}.
    About the behaviour of argument \code{onefile} in
    \command{\link{pdf}}

    Default: \code{FALSE}.
  }
  
 \item{\code{width}}{real number or NA; only relevant if
    \code{always_open_screen=TRUE}. See \code{height} for details.
       
    Default: \code{NA}.
  }

  \item{\code{xlim}}{real; default right limit point of model plots
       
    Default: \code{2.0}.
  }
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bold{11. \code{gui}: Options for c\command{\link{RFgui}}}
\describe{
  \item{\code{alwaysSimulate}}{
    logical. If \code{TRUE} then a new random field is simulated
    whenever a parameter is changed. Otherwise only the covariance
    function or the variogram is re-plotted; simulations are performed
    only when the correponding button is pressed.
    
    Default: \code{TRUE}.
  }
  \item{\code{simu_method}}{
    \code{"RPcirculant"},
    \code{"RPcutoff"},
    \code{"RPintrinsic"},
    \code{"RPtbm"}, 
    \code{"RPspectral"},
    \code{"RPdirect"},
    \code{"RPsequential"},
    \code{"RPaverage"},
    \code{"RPnugget"},
    \code{"RPcoins"},
    \code{"RPhyperplane"},
    \code{"RPspecific"},
    \code{"any method"}.
    
				      
    Default: \code{"RPcirculant"}.
  }
  \item{\code{size}}{vector of 2 components.
    Grid size of the simulated stochastic processes.
    The two components of the vector correspond to one-dimensional and
    two-dimensional processes, respectively.
    
    Default: \code{c(1024, 64)}.
  }
}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bold{12. \code{hyper}: Options for simulating hyperplane tessellations}
\describe{
  \item{\code{mar_distr}}{integer.
    This argument should not be changed yet.
    
    It  codes the marginal distribution used in the
    simulation:
    
    0 : uniform distribution\cr
    1 : Frechet distribution with form argument \code{mar_param}\cr
    2 : Bernoulli distribution (Binomial with \eqn{n=1}) with
    argument \code{mar_param}
    
    
    Default: \code{0} .
  }

  \item{\code{mar_param}}{Argument used for the marginal
    distribution. The argument should not be changed yet.
    
    Default: \code{NA} .
  }

  \item{\code{maxlines}}{integer.
    Maximum number of allowed lines.
    
    Default: \code{1000} .
  }
  
  \item{\code{superpos}}{integer.
    number of superposed hyperplane tessellations.
    
    Default: \code{300} .
  }
}



 
\bold{13. \code{krige}: Options for Kriging}
\describe{
  
 % \item{\code{apart_linear_fit}}{
%    logical value in case the linear part
%    (and/or the global variance) of the model are not determined yet,
%    i.e. contain \code{NA}s.
 %   
%    If \code{TRUE} repeated measurements are always split up
%%2    into different data sets, to which different linear parts of
%    the model are fitted. See also \command{\link{RFinterpolate}}.
%    
%    Default: \code{FALSE}.
%  }
 
  \item{\code{cholesky_R}}{ obsolete }

  \item{\code{fillall}}{
    logical value for imputing.
    If true all the components are estimated whether they are
    \code{NA} or not.
    
    Default: \code{TRUE}.
  }
   
  \item{\code{locmaxn}}{
    Kriging is conditions on maximal \code{locmaxn} points.
    If the data contain more points, neighbourhood kriging is performed.
    
    Default: \code{8000}. 
  }
  
  \item{\code{locsplitfactor}}{
    In case of neighbourhood kriging, the area is split into small
    boxes. The complete neighbourhood contains (2 *
    \code{locsplitfactor} +1) boxes in each direction.
    
    
    Default: \code{2}. 
  }
  
%  \item{\code{locsplitn}}{vector of 3 components.
%    A box should contain no more than \code{locsplitn[1]}
%    points, but never less than \code{locsplitn[2]}. If
%    a box had originally less than \code{locsplitn[2]} points,
%    then the box is increased until at least \code{locsplitn[3]}
%    points are in the box.    }
   \item{\code{locsplitn}}{vector of 3 components.
    A box should contain no more than \code{locsplitn[3]}
    points, but never less than \code{locsplitn[1]}. If
    a box had originally less than \code{locsplitn[1]} points,
    then the box is increased until at least \code{locsplitn[2]}
    points are in the box.
   
    Default: \code{c(200, 1000, 5000)}. 
  }
  
  \item{\code{method}}{ obsolete }
  
  \item{\code{return_variance}}{logical.
    If \code{FALSE} the kriged field is
    returned. If \code{TRUE} a list of two elements, \code{estim} and
    \code{var}, i.e. the kriged field and the kriging variances,
    is returned.
    
    Default: \code{FALSE}. 
  }
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bold{14. \code{maxstable}: Options for simulating max-stable random fields}
\describe{ % extreme
  \item{\code{check_every}}{
    integer. In order to get a precise simulation result, by definition,
    the maximum must be taken, for each shape function, over alle
    locations of interest. Clearly, small values will not play a role.
    To this end, the global minimum has to be determined.
    The calculation of the global minimum is expensive and therefor
    should not be done too frequently. On the other hand,
    rare updates increases the computing times for taking the maximum
    over a single shape functions. Here, after every \code{check_every}
    considered shape function, the global minimum is calculated.
    It is expected that a good choice for  \code{check_every} is in
    in the interval \eqn{[10, 100]}.
     
    (For ease and for concerns of efficiency, the more adequate, local
    minimum is not considered.)
    
    Default: 30 .
  }
    
  \item{\code{density_ratio}}{
    value in \eqn{[0,1]}. This argument is considered only
    if \code{flat=-1} and the simulation is performed on a grid.
    Then, the ratio between the highest and the lowest value is
    calculated within the convex hull of the grid. If the
    value is less than \code{density_ratio} then the grid points
    are considered separately. Else the density is considered to be
    constant in the convex hull of the grid.
    
    Default: 0.0.
  }
  
  \item{\code{eps_zhou}}{positive real number, which 
    gives the aimed relative precision when the constant \eqn{c}
    in the paper of Oesting, Schlather, Zhou (2018) has to be estimated.
    E.g. if \code{eps_zhou=0.01} then the first 2 digits should be
    correct.
   
    Default: 0.01
  }


  \item{\code{flathull}}{\code{NA, FALSE, TRUE}. Only used in M3
    modelling in the algorithm by Oesting, Schlather, Zhou (2018).
    The argument is considered only if the simulation is performed on a
    grid. 
    If \code{flat=TRUE} , then the density
    is considered to be \code{flat} in the convex hull of the grid,
    i.e. the simulation method of Schlather (2002) is used.
    If \code{flat=NA} the choice is done automatically.
    
    Default: FALSE .
  }
  

  \item{\code{max_gauss}}{
    The simulation of the max-stable process by the old-fashioned
    method of Schlather (2002) and by older methods for Brown-Resnick
    processes 
    uses
    a stopping rule that necessarily needs a finite upper endpoint
    of the marginal distribution of the random field.
    In the case of
    \command{\link[=RPbrownresnick]{Brown-Resnick processes}},
    \command{\link[=RPschlather]{extremal Gaussian fields}},
    and
    \command{\link[=RPopitz]{extremal t fields}},
    the upper endpoint is approximated by \code{standardmax}. 
    
    Default: \code{3.0} .
  }
  
  \item{\code{max_n_zhou}}{positive integer.
    The overall constant \eqn{c} in the paper of
    Oesting, Schlather, Zhou (2018) has to be determined
    by MCMC, if the shape functions are random.
    
    The two arguments, \code{min_n_zhou} and \code{max_n_zhou},
    give the minimal and the maximal
    number of simulations that are performed. To economize
    computer time the values of \eqn{c} is partially estimated 
    when the shape functions are simulated. If the number
    of shape functions is larger than the number of simulations
    given by \code{eps_zhou} then
    no further simulation is performed to determine \eqn{c}.
    So, it is advantageous to simulate all fields at once by 
    \code{RFsimulate(..., n = )}. 
   
    Default:  1000 and 10000000, respectively.
  }

  \item{\code{maxpoints}}{
    positive integer; the maximal number of Poisson points to be simulated
    for one realization of the max-stable random field. This option will
    not be considered for most of the users. This option allows 
    the simulation to interrupt after \code{maxpoints} shape function
    have been placed.
     
    Default: \code{2e9} (never).
  }

  \item{\code{mcmc_zhou}}{positive integer.
    In case of random shape functions, an MCMC step is required.
    \code{mcmc_zhou}-1 equals the number of members of the MCMC chain
    that are left out before the next value of the chain is returned.
    
    Default: 20
  }
  
  \item{\code{min_n_zhou}}{see \code{max_n_zhou}}

 
  \item{\code{mcmc_zhou}}{positive integer.
    In case of random shape functions, an MCMC step is required.
    \code{mcmc_zhou}-1 equals the number of members of the MCMC chain
    that are left out before the next value of the chain is returned.
    
    Default: 20
  }
  
  \item{\code{min_n_zhou}}{see \code{max_n_zhou}}

  \item{\code{min_shape_gumbel}}{To increase speed,
    the minimum field value is assumed to be \code{min_shape_gumbel}
    for calculation of threshold values for simulation short cuts.
    During a simulation, its value becomes void as soon as the real
    (current) minimum of the field being simulated exceeds
    \code{min_shape_gumbel}
%    values of shape
%    functions smaller than \code{min_shape_gumbel} are considered as
%    being zero. Necessarily, the shape function must be essentially
%    isotropic and decreasing with known inverse function.

    Default: \code{-1e15}.
  }

  \item{\code{scatter_method}}{logical. 
    If 
    
    Default: NA;
  }
  
  \item{\code{xi}}{
    Extreme value index. Default: \code{2e9} .
    While \eqn{\xi} can be set globally, the shift \eqn{\mu} and the
    scale \eqn{s} can be given only locally within the process
    definitions, e.g., \command{\link{RPsmith}}.

   Default: \code{1.0}.
  }

}
 
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bold{15. \code{mpp}: Options for the random coins (shot noise) methods}
\describe{
\item{\code{about_zero}}{
    In certain cases (\link{Coins},\link{RMtruncsupport}),
    functions are assumed to zero if the value is less than \code{about_zero}.
    
    Default: \code{0.001} .
  }
  \item{\code{n_estim_E}}{integer. Number of draws from the
    distribution of the scale to estimate the mean of the distribution.
    This is used only if the mean of the scale distribution
    is not explicitely given.
    Default: \code{50000} .
  }

  \item{\code{scatter_method}}{
  }
  
  \item{\code{scatter_size}, \code{scatter_max}}{
    Real valued and integer valued, respectively, or \code{NA}.

    Used in the internal function \code{RMscatter} that calculates
    \eqn{\sum_{i=1}^n f(x + h_i)} for some function \eqn{f} and
    for some distances \eqn{h_i}.
      
    %    see \code{size} and \code{max} in \code{\link{RMscatter}}.
    Let \eqn{\varepsilon=}\code{about_zero}, \eqn{s=}\code{scatter_size} and \eqn{m=}\code{scatter_max}.
    We distinguish 4 cases:
    \itemize{
      \item \code{scatter_size > 0} and \code{scatter_max >= 0}\cr
      Here, \eqn{n} equals \eqn{(2m)^d}.
      and \eqn{h_i \in M = \{ (k s, \ldots, k s),\ldots, (m s, \ldots, m
	s)\}}
      with \eqn{k=-m}.
      
      \item \code{scatter_size > 0} and \code{scatter_max < 0}\cr
       same as the previous case, but \eqn{m} is chosen such that
      \eqn{f(k_i e_i s_i) \approx \varepsilon}, \eqn{-k_i\in N},
      \eqn{i=1,\ldots,d} and 
      \eqn{f(m_i e_i s_i) \approx \varepsilon}, \eqn{m \in N}.
      
      \item \code{scatter_size <= 0} and \code{scatter_max >= 0}\cr
      This option is possible only for grids.
      Here \eqn{h_i} runs on the given grid \eqn{i=1,\ldots,d},
      but at most \code{scatter_max} steps.
      
      \item \code{scatter_size <= 0} and \code{scatter_max < 0}\cr
      this option is possible only for grids.
      Here, \eqn{h_i} runs over the whole grid.
    }
  }
  \item{\code{shape_power}}{
    Shape functions are powered by \code{shape_power} before used as
    intensity function for the point process.
    
    Default: \code{2.0}.
  }
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bold{16. \code{nugget}: Options for the nugget effect}\cr
Simulating a nugget effect is per se trivial.
However, it gets complicated
and best methods (including \code{direct} and \code{circulant
  embedding}!) fail if zonal anisotropies are considered,
where sets of points have to be identified that belong to the
same subspace of eigenvalue 0 of the anisotropy matrix.
\describe{
  \item{\code{tol}}{
    The nugget tolerance influences two different kind of models
    \itemize{
      \item \command{\link{RPnugget}}
      \item \command{\link{R.is}}
    }
    See there for more information.
  }
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bold{17. \code{registers}: Register numbers}\cr
Model for different purposes are or can be stored at
different places. They are called registers and have non-negative
numbers up to 21 (currently).
The user can use the registers 0..9.
\describe{
  \item{\code{register}}{number in 0:9; place where intermediate calculation
    for random field simulation are stored; 
    the number refers to 10 internal registers 0..9.
    
    Changing the register number only makes sense, when
    two different random fields, say, are to be simulated
    alternatingly, several times in a row. Then the
    simlulation speed can be increased if several registers
    are used, \code{storing=TRUE} and \command{\link{RFsimulate}}
    is used with the only argument \code{n}.
     
    Default: \code{0} %[also do].
  } 
}




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bold{18. \code{sequ}: Options for the sequential method}
\describe{
  \item{\code{back_steps}}{See \command{\link{RPsequential}}}
  \item{\code{initial}}{See \command{\link{RPsequential}}}
  \item{\code{max_variables}}{See \command{\link{RPsequential}}}
}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bold{19. \code{solve}: Options for solving linear systems}

\describe{
  \item{\code{det_as_log}}{\RFU }
  \item{\code{eigen2zero}}{\RFU }
  \item{\code{max_chol}}{integer. Maximum number of rows of a matrix in
    a Cholesky decomposition
    
    Default: \eqn{8192}
  }
  \item{\code{max_svd}}{integer. Maximum number of rows of a matrix in
    a svd decomposition
    
    Default: \eqn{6555}
  }
    
  \item{\code{pivot}}{Type of pivoting for the Cholesky
      decomposition. Possible values are 
      \describe{
	\item{PIVOT_NONE}{No pivoting.}
	\item{PIVOT_AUTO}{If the  matrix has a size greater than
	  3x3 and Choleskey fails without pivoting, privoting
	  is done. For matrices of size less than 4x4, no pivoting and
	  no checks are performed.}
 	\item{PIVOT_DO}{Do alwaoys pivoting.
	  NOTE: privoted Cholesky decomposition yields only very approximately
	  an upper triangular matrix L, but still L^t L = M holds true.}
 	\item{PIVOT_IDX}{uses the same pivoting as in the previous
	  pivoted decomposition. This option becomes relevant only when
	  simulations with different parameters or different models shall be
	  performed with the same seed so that also the pivoting must be
	  coupled.
	}
% 	\item{PIVOT_IDXBACK}{ same as \code{PIVOT_IDX}, but
%	the sequence of indices of the pivoting is returned via
%	\code{RFoptions()$solve$pivot_idx}.}
      }

      Default: \code{PIVOT_auto}}  
  \item{\code{pivot_actual_size}}{\RFU }     
  \item{\code{pivot_check}}{logical. Only used in pivoted Cholesky
    decomposition. 
    If \code{TRUE} and a numerically zero diagonal element is detected,
    it is checked whether the offdiagonal elements are numerically zero 
    as well.
    (See also \code{pivot_max_deviation} and
    \code{pivot_max_reldeviation}.)
    if \code{NA} then, in  \link{RPdirect}, the value is equivent to
    \describe{
      \item{\code{FALSE}}{if the model is positive (semi-)definite.}
      \item{\code{TRUE}}{if the model is genuinely negative definite.}
    }
    
    Default: \code{NA}
  }

  \item{\code{pivot_idx}}{\RFU}
  \item{\code{pivot_relerror}}{\RFU  }    
  \item{\code{pivot_max_deviation}}{\RFU }    
  \item{\code{pivot_max_reldeviation}}{\RFU  }   
  \item{\code{solve_method}}{\RFU }
  \item{\code{spam_factor}}{\RFU  }    
  \item{\code{spam_min_n}}{\RFU  }
  \item{\code{spam_min_p}}{\RFU  }
  \item{\code{spam_pivot}}{\RFU}
  \item{\code{spam_sample_n}}{\RFU}  
  \item{\code{spam_tol}}{\RFU}
  \item{\code{svdtol}}{\RFU}
  \item{\code{use_spam}}{\RFU}
}
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bold{20. \code{special}: Options for specific methods}
\describe{
  \item{\code{multicopies}}{Only used by \command{\link{RMmult}}.
    The covariance functions are multiplied if the corresponding
    independent random fields are multiplied. To get
    an approximative Gaussian random fields with a multiplicative
    covariance functions the average over \code{multicopies}
    products of random fields is calculated.
  }
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bold{21. \code{spectral}: Options for the spectral (turning bands) method}
\describe{
  \item{\code{ergodic}}{
    In case of an additive model and \code{ergodic=FALSE},
    the additive component are chosen proportional to their
    variance. In total \code{lines} are simulated. If
    \code{ergodic=TRUE}, the components are simulated
    separately and then added.
    
    Default: \code{FALSE}.
  }

 \item{\code{prop_factor}}{see \command{\link{RPspectral}}
    
  }

    
  \item{\code{sigma}}{see \command{\link{RPspectral}}
  }

  \item{\code{sp_grid}}{   see \command{\link{RPspectral}}
  }
 
  \item{\code{sp_lines}}{
    see \command{\link{RPspectral}}
  }
 
}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\bold{22. \code{tbm}: Options for the turning bands method}
\describe{
  \item{\code{center}}{Scalar or vector.
    If not \code{NA}, the \code{center} is used as the center of
    the turning bands for \code{TBM2} and \code{TBM3}.
    Otherwise the center is determined
    automatically such that the line length is minimal.
    See also \code{points} and the examples below.
    
    Default: \code{NA} .
  }

  \item{\code{fulldim}}{
    positiv integer. The dimension of the space into which the
    simulated field is embedded. So, the value \code{fulldim}
    must be at least the dimension of the field.

    Default: 3.
  }
  
  \item{\code{grid}}{Logical. 
    The angle of the lines is random if
    \code{grid=FALSE}, 
    and \eqn{k\pi/}{k*pi/}\code{lines}
    for \eqn{k}{k} in \code{1:lines},
    otherwise.
    
    This option is used by both \command{\link{RPspectral}}
    and \command{\link{RPtbm}}, the latter only when the dimension is 2.
    
    Default: \code{TRUE} .
  }
  
  \item{\code{layers}}{
    Logical or integer. If \code{TRUE} then the turning layers are used whenever
    a time component is given.
    If \code{NA} the turning layers are used only when the
    traditional TBM is not applicable.
    If \code{FALSE} then turning layers may never be used.
    
    Default: \code{TRUE} .
  }

  \item{\code{lines}}{
    Number of lines used.
    
    Default: \code{60} .
  }
  
  \item{\code{linesimustep}}{
    If \code{linesimustep} is positive the grid on the line has lag
    \code{linesimustep}. 
    See also \code{linesimufactor}.
    
    Default: \code{0.0} .
  }

 \item{\code{linesimufactor}}{ \code{linesimufactor} or
    \code{linesimustep} must be non-negative; if
    \code{linesimustep}
    is positive then \code{linesimufactor} is ignored.
    If both
    arguments are naught then \code{points} is used (and must be
    positive).
    The grid on the line is \code{linesimufactor}-times
    finer than the smallest distance. 
    See also \code{linesimustep}.
    
    Default: \code{2.0} .
  }

  \item{\code{points}}{integer. If greater than 0,
    \code{points} gives the number of points simulated on the TBM
    line, hence 
    must be greater than the minimal number of points given by
    the size of the simulated field and the two paramters
    \code{TBMx.linesimufactor} and \code{TBMx.linesimustep}.
    If \code{points} is not positive the number of points is
    determined automatically.
    The use of \code{center} and \code{points} is highlighted
    in an example below.
    
    Default: \code{0}.
  }

  \item{\code{reduceddim}}{
    if positiv integer, then the value itself. If negativ, then
    the value is substracted from fulldim.

    Default: \code{-2}.
  }
  
}
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 \bold{23. \code{\link{messages}}: Options for specific messages}

 See \code{\link{messages}}.





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 \bold{24.A Internal options with prefix}
  Internal options should not be changed by the user.
  \describe{    
    \item{\code{do_tests}}{
      Internal variable. Do not use it.
      Default: \code{FALSE}.
    }

    \item{\code{examples_reduced}}{non-negative integer.
      If positve, then the design of any simulation in \pkg{RandomFields}
      is internally reduced in size (roughly downto the given value in each
      direction). Warnings report this behaviour.
      This option is necessary to run the examples of \pkg{RandomFields}
      under the time constraint of CRAN.
    }

    \item{\code{warn_unknown_option}}{integer. See
      \code{?RandomFieldsUtils::RFoptions} for details.
      
      Default for \pkg{RandomFields}: \code{1}
    }
  }
    
\bold{24.B Internal options without prefix}
  \describe{    
    \item{\code{getoptions_}}{(from \pkg{RandomFieldsUtils})
      string vector of prefixes that indicate
      classes of options. In this package it
      can be \code{"general"} and \code{"fit"}, for instance.
      See \code{str(RFoptions(), max=1, give.attr=FALSE)} to get the
      whole list of options.
      The given classes of options are then
      returned by \code{RFoptions()}. Note that the values are the
      previous values.
      
      \code{getoptions_} must always be the first argument.
    }
    \item{\code{local_}}{logical.
      In parallel computing, each instance of R (identified by its pid),
      gets its own copy of the set of of options.
      If \code{TRUE}, \code{RFoptions} will refer to the local set.
      Note that \code{local_} must be set as the very first argument,
      even before \code{getoptions_}.
      
      Default: \code{FALSE}.
    }
    \item{\code{parent}}{environment, that should be used to solve the
    model definition}
  }
%\bold{General comments}\cr
%Most of the above arguments determine the basic settings of a
%simulation, e.g. \code{matrix_method} (which chooses the method to
%calculate a square
%root of a positive definite matrix). The values of
%such arguments are read by
%\command{\link{RFsimulate}} and stored in an internal register
%and used in preliminary, deterministic calculations (e.g. calculation
%of the square of a covariance matrix)
%Some few other arguments like \code{lines} (which determines the number of
%i.i.d. processes to be simulated on the line)
%are only relevant when generating
%random numbers.
%These arguments are read only when actually simulations are performed, and
%are marked by \dQuote{[do]}.

%The difference of these two classes of arguments might become
%important only if\code{general.storing = TRUE} and
%\command{\link{RFsimulate}} is called with the only argument
%\code{n}.

}
 
 
\section{Value}{
 \code{NULL} if any argument is given, and the full list of
 arguments, otherwise.
}

\references{
  \itemize{Basic}
  \itemize{
    \item General
    \itemize{
      \item
      Schlather, M. (1999) \emph{An introduction to positive definite
	functions and to unconditional simulation of random fields.}
      Technical report ST 99-10, Dept. of Maths and Statistics,
      Lancaster University.

      \item Schlather, M. (2011) Construction of covariance functions and
      unconditional simulation of random fields. In Porcu, E., Montero, J.M.
      and Schlather, M., \emph{Space-Time Processes and Challenges Related
	to Environmental Problems.} New York: Springer.
      % \item Schlather, M. (2002) Models for stationary max-stable
    }

    \item rectangular distribution; \code{eps_zhou}
    \itemize{
      \item
      Oesting, M., Schlather, M. and Zhou, C. (2013) On the Normalized
      Spectral Representation of Max-Stable Processes on a compact set. 
      \emph{arXiv},  \bold{1310.1813}
    }

    \item \code{shape_power}
    \itemize{
      \item
      Ballani, F. and Schlather, M. (2015) In preparation.
    }
  }
}

\me

\seealso{
  \link{fit},
  \link{messages},
  \command{\link{RFsimulate}},
  \link{RFoptionsAdvanced},
 \code{\link[RandomFieldsUtils]{RFoptions}},
 and \command{\link{RFgetMethodNames}}.}

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

RFoptions()

\dontshow{\dontrun{
############################################################
##                                                        ## 
## use of exactness                                       ##
##                                                        ##
############################################################
x <- seq(0, 1, 1/30)
model <- RMgauss()

for (exactness in c(NA, FALSE, TRUE)) { 
  readline(paste("\n\nexactness: `", exactness, "'; press return"))
  z <- RFsimulate(model, x, x, exactness=exactness,
                  stationary_only=NA, storing=TRUE)
  print(RFgetModelInfo(RFsimulate,which="internal")) # $internal$name)
}
}}

\dontshow{FinalizeExample()}}

\keyword{spatial}

