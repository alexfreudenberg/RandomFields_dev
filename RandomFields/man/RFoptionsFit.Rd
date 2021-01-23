\name{RFoptionsFit}
\alias{RFoptionsFit}
\alias{fit}
\title{Setting control arguments for \command{\link{RFfit}},
    \command{\link{RFratiotest}}, and \command{\link{RFcrossvalidate}}}
\description{
 \command{\link{RFoptions}} sets and returns control arguments for the analysis
 and the simulation of random fields. It expands the functionality of
 \link[RandomFieldsUtils]{RFoptions}. Usage:\cr
 \code{RFoptions(...)}
}


\section{Details}{Options for \command{\link{RFfit}},
   \command{\link{RFratiotest}}, and \command{\link{RFcrossvalidate}}


\describe{
  \item{addNAlintrend}{logical or non-negative integer.
    If \code{FALSE} then all linear coefficients to be estimated must be
    given explicitely. If \code{TRUE} then only the constant \code{1} is
    replaced by \code{NA}. If the value is greater than eqn{1}
    then, for any linear term without
    leading constant, a linear coefficient is estimated.

    Default: \code{2}
    %\code{TRUE} in \pkg{RandomFields} and \code{2} for 
%    functions in \pkg{RandomFieldsLight}
  }

  \item{\code{algorithm}}{
    See \link{RFfitOptimiser}.

    Default: \code{NULL}
  }
  
  \item{\code{approximate_functioncalls}}{
    In case the parameter vector is too close to the given
    bounds, the ML target function is evaluated on a grid
    to get a new initial value for the ML estimation. 
    The number of points of the grid is approximately
    \code{approximate_functioncalls}.
 
    Default: \code{50}
  }

  \item{\code{boxcox_lb}}{
    lower bound for the Box-Cox transformation
    
    Default: \code{-10}.
  }
  \item{\code{boxcox_ub}}{
    upper bound for the Box-Cox transformation
    
    Default: \code{10}.
  }

  \item{\code{bin_dist_factor}}{
    numeric. The empirical variogram is calculated up the distance
    \code{bin_dist_factor} times (maximum distance among any pair of locations)
    
    Default: \code{0.5}.
  }
  
  \item{\code{bins}}{vector of explicit boundaries for the bins or the
    number of bins for the empirical variogram (used in the
    LSQ target function, which is described at the beginning
    of the Details).
    Note that for anisotropic models, the value of \code{bins} might
    be enlarged.
    
    Default: \code{20}.
  }

 \item{\code{critical}}{logical or signed integer.
    
    If \code{critical=FALSE} and if the result of
    any maximum likelihood method 
    is on a borderline, then the optimisation is redone
    in a modified way (which takes about double extra time)
    
    If \code{critical=TRUE} and if the result of
    any maximum likelihood method
    is on a borderline, then a kind of profile likelihood
    optimization is done (which takes about 10 times extra time)
    
    If \code{critical>=2} then a kind of profile likelihood
    optimization is always done (which takes about \code{n_crit}
    times extra time) for an automatically chosen selection
    of the model parameters.
    
    If \code{critical>=3} then a kind of profile likelihood
    optimization is always done (which takes about \code{n_crit}
    times extra time) for all the parameters.
    
    If \code{critical<0} then none of the refined methods
    are performed.
    
    Default: \code{TRUE}.
  }

  \item{\code{cross_refit}}{logical. 
    For each of the subset of the cross-validation method
    the parameters have to be fitted to the given model.
    If \code{cross_refit} is \code{TRUE}, this is done, but takes a huge
    amount of time. If \code{FALSE}, the model is fitted only once to
    the data and the value at
    each point is predicted with the same model given
    the values of the other points.
    
    Default: \code{FALSE}.
  }

 \item{\code{emp_alpha}}{integer within \code{RC_VARIOGRAM},
   \code{RC_PSEUDO},  \code{RC_PSEUDOMADOGRAM},
   \code{RC_MADOGRAM},
   \code{RC_COVARIANCE} or a real value within \eqn{(0,2]}
   for (cross-)pseudovariograms.
   This controls the kind of performed least squares fit to the
   empirical estimates.

    Default: \code{RC_VARIOGRAM}
  }
  \item{estimate_variance}{
    see \command{\link{RFlikelihood}}.
  }

  \item{factr, factr_recall}{
    See the argument \code{control} in \link[stats]{optim}.
    \code{factr_recall} is used for intermediate calculations.
  }
  
  \item{likelihood}{character -- not programmed yet.
    types of likelihood are \code{"auto"}, \code{"full"},
    \code{"composite"}, \code{"tesselation"}; 
    
    Default: \code{"auto"}    
  }

  \item{\code{lowerbound_scale_factor}}{
    The lower bound for the scale is determined as
    
    (minimum distance between different pairs of points) /\cr
     \code{lowerbound_scale_factor}.
    
    Default: \code{3}.
  }

  \item{\code{lowerbound_scale_ls_factor}}{ For the LSQ target
    function a different lower bound 
    for the scale is used. It is determined as
    
    (minimum distance between different pairs of points) / \cr
    \code{lowerbound_scale_ls_factor}.
    
    Default: \code{5}.
  }

 % \item{\code{lowerbound_sill}}{absolute lower bound for variance
%    and nugget. See \code{lowerbound_var_factor}.
%    
%    Default: \code{1E-10}.
%  }
  
  \item{\code{lowerbound_var_factor}}{
    The lower bound for the nugget and the variance is determined
    as var(\code{data}) / \code{lowerbound_var_factor}.
    If a standard model definition is given and
    either the nugget or the variance is fixed,
    the parameter to be estimated
    must also be greater than \code{lowerbound_sill}.
    
    Default: \code{10000}.
  }

  \item{\code{maxmixedvar}}{OBSOLETE.
    upper bound for variance in a mixed model;
    so, the covariance model for mixed model part might
    be calibrated appropriately
  }

  \item{\code{max_neighbours}}{integer.
    Maximum number of locations (with depending values)
    that are allowed.
    
    Default: \code{5000}.
  }
  
  \item{\code{minbounddistance}}{
    If any value of the parameter vector
    returned from the ML estimation
    is closer than \code{minbounddistance}
    to any of the bounds or if any value
    has a relative distance smaller than
    \code{minboundreldist}, then it is assumed that
    the MLE algorithm has dropped into a local minimum,
    and it will be continued with evaluating the
    ML target function on a grid, cf. the beginning paragraphs
    of the Details.
    
    Default: \code{0.001}.
  }
  
  \item{\code{minboundreldist}}{relative distance to the bounds
    below which a part of the algorithm is considered as
    having failed. See \code{minbounddistance}.
    
    Default: \code{0.02}.
  }

  \item{\code{min_diag}}{
    Minimal value of any estimated diagonal matrix element.    
    
    Default: \code{1e-7}.
  }

  \item{\code{n_crit}}{integer.
    The approximate profiles that are considered.
    
    Default: \code{10}.
  }

    \item{\code{nphi}}{scalar or vector of 2 components.
    If it is a vector then the first component gives the first angle
    of the xy plane
    and the second one gives the number of directions on the half circle.
    If scalar then the first angle is assumed to be zero.
    Note that a good estimation of the variogramm by LSQ with a
    anisotropic model a large value for \code{ntheta} might be needed
    (about 20).
    
    Default: \code{1}. 
  }
  \item{\code{ntheta}}{scalar or vector of 2 components.
    If it is a vector then the first component gives the first angle
    in the third direction
    and the second one gives the number of directions on the half circle.
    If scalar then the first angle is assumed to be zero.
    
    Note that a good estimation of the variogramm by LSQ with a
    anisotropic model a large value for \code{ntheta} might be needed
    (about 20). 
    
    Default: \code{1}. 
  }
  \item{\code{ntime}}{scalar or vector of 2 components.
    if \code{ntimes} is a vector, then the first component are the
    maximum time distance (in units of the grid length \code{T[3]}) and the
    second component gives the step size (in units of the grid length
    \code{T[3]}). If scalar then the step size is assumed to 1 (in units
    of the grid length \code{T[3]}).
    
    Default: \code{20}. 
  }

  \item{\code{only_users}}{boolean.
    If true then only \code{users_guess} is used as a
    starting point for the fitting algorithms
    
    Default: \code{FALSE}. 
  }

  \item{\code{optimiser}}{
    See \link{RFfitOptimiser}.
    
    Default: \code{"optim"}.
  }
  
% \item{\code{optim_var_elimination}}{This argument takes the values
%    \code{'never'}, \code{'respect bound'}, \code{'try'},
%    \code{'yes'}, and should only be 
%    set by the advanced user. Background of this option is that
%    a global variance can optimized analytically.
%    
%    The meaning of the values is as follows.
%    \itemize{
%      \item{\code{'never'}}{
%	A global variance is never tried to be eliminated
%      }
%      \item{\code{'respect bound'}}{
%	A global variance is eliminated if such a variance
%	is detected and the user did not indicate bounds
%	for the parameters.
%      }
%      \item{\code{'try'}}{
%	A global variance is eliminated if such a variance
%	is detected.
%      }
%      \item{\code{'yes'}}{
%	A global variance is tried to be eliminated altough
%	the algorithm did not find an indication. Here, the full
%	responsibility is left to the user. (This option might
%	make sense if \code{transform} is given.)
%	This option is only overwritten when it does not make sense,
%	e.g. no variance is estimated.
%      }
%    }
%    Default: \code{'respect bound'}.
%  }

  \item{pgtol, pgtol_recall}{
    See the argument \code{control} in \link[stats]{optim}.
    \code{pgtol_recall} is used for intermediate calculations.
  }
  
  \item{\code{refine_onborder}}{logical.
    If \code{TRUE} and an estimated parameter of the model
    is close to the boundary, a second search for the optimum
    is started.

    Default: \code{TRUE}
    
  }
 
  \item{\code{minmixedvar}}{
    lower bound for variance in a mixed model;
    so, the covariance model for mixed model part might
    be calibrated appropriately

    Default:  1/1000
  }
  
%  \item{\code{solvesigma}}{Logical. -- experimental stage!
%    If a mixed effect part is present where the variance
%    has to be estimated, then this variance parameter is solved
%    iteratively within the profile likelihood function, if
%    \code{solvesigma=TRUE}.This makes sense
%    if the number of independent variables is very small.
%    If \code{solvesigma=FALSE} then the variance parameter is
%    treated as any other parameter to be estimated.
%    Default: \code{FALSE}.
%  }
  
  \item{\code{ratiotest_approx}}{logical.
    if \code{TRUE} the approximative formula that twice the
    difference of the likelihoods follow about a \eqn{\chi^2}
    distribution is used. The parameter of freedom equals
    the number of parameters to be estimated for the covariance
    function, including those for the covariates.

    Default: \code{TRUE}
  }
  \item{\code{reoptimise}}{logical.
    If \code{TRUE && !only_users} then at a very last step,
    the optimisation is redone with currently best parameters
    and likelihood as scale parameter for \command{\link{optim}}.
    
    Default: \code{TRUE}. 
  }

    \item{\code{return_hessian}}{logical.
      If \code{TRUE} then the hessian matrix is calculated
      for the MLE and returned.
    
      Default: \code{TRUE}. 
  }

  

  \item{\code{scale_max_relative_factor}}{ If the initial scale
    value for the ML estimation 
    obtained by the LSQ target function is
    less than
    \eqn{(minimum distance
      between different pairs of points) / }
    \code{scale_max_relative_factor}
 
    a warning is given that probably a nugget effect
    is present. 
    Note: if \code{scale_max_relative_factor} is greater
    than \code{lowerbound_scale_ls_factor} then
    no warning is given as
    the scale has the lower bound \eqn{(minimum distance
      between different pairs of points) / }
    \code{lowerbound_scale_ls_factor}.
    
    Default: \code{1000}
  }

  \item{\code{scale_ratio}}{
    \command{\link{RFfit}} uses \code{parscale} and \code{fnscale}
    in the calls of \command{\link{optim}}. As these arguments should
    have the magnitude of the estimated values, \command{\link{RFfit}}
    checks this by calculating the absolute log ratios.
    If they are larger than \code{scale_ratio},
    \code{parscale} and \code{fnscale} are reset and the optimisation
    is redone. 
    
    Default: \code{0.1}. 
  }
   
  \item{\code{shortnamelength}}{
    The names of the variables in the returned table are
    abbreviated by taking the first \code{shortnamelength}
    letters.
    
    Default: \code{4}. 
  }
 
%  \item{\code{sill}}{ currently not maintained anymore.    
%    Additionally to estimating \code{nugget} and \code{variance}
%    separately, they may also be estimated together under the
%    condition that \code{nugget} + \code{variance} = \code{sill}.
%    For the latter a finite value for \code{sill} has to be supplied,
%    and \code{nugget} and \code{variance} are set to \code{NA}.
%    
%    \code{sill} is only used for the standard model. 
%    
%    Default: \code{NA}. 
%  }

  \item{\code{smalldataset}}{
    If the number of locations is considered as small, then some more data
    are kept in the storage to accelerate the estimation algorithm.
    
    Default: \code{2000}.
  }
  
  \item{\code{split}}{integer.
    If the number of parameters to be numerically optimised is larger
    than or equal to \code{split} then \command{\link{RFfit}} checks whether a
    space-time covariance model or a multivariate covariance model
    can be split into components, so that certain parameters
    can be estimated separately.
    
    Default: \code{4}. 
  }
   
  \item{\code{cliquesize}}{integer.
    \command{\link{RFfit}} tries to split the data set
    into parts of size splitn_neighbours[2] or less, but never more than 
    \code{splitn_neighbours[3]} and never less than
    splitn_neighbours[1].     
    
    Default: \code{c(200, 1000, 3000)}.
  }
  
  \item{\code{splitfactor_neighbours}}{
     The total number of neighbouring boxes in each direction
     \eqn{1 + 2\code{splitfactor}}, including the current box itself.
    
    Default: \code{2}.
  }
  \item{\code{split_refined}}{logical.
    If \code{TRUE} then also submodels are fitted if splitted.
    This takes more time, but \command{\link[=anova.RF_fit]{anova}} and
    \command{\link{RFratiotest}}, for instance,
    will give additional information.

    Default: \code{TRUE}.
  }

   \item{\code{suggesting_bounds}}{logical.
    If \code{FALSE} then upper and lower bounds for the unknown
    parameters that are given by the user 
    are adopted as such. If \code{TRUE} these bounds are adopted only if
    they narrow the interval calculated internally. If
    the interval given by the user is outside the internally calculated
    interval, the user's interval is always adopted.
    
    Default: \code{FALSE}.
  }
  
  \item{\code{trace}}{integer.
    Tracing with \command{\link{RFfit}}.
    If negative, numbers are shown.
    If positive and at least one LSQ method is used, pictures are shown.
    The larger the absolute number, the more more tracing is shown.
    If the modulus is one, only the tracing is updated only if a new
    optimum is found.
    If the modulus equals 2, all function calls are traced.
   
    \code{upperbound_scale_factor} * (maximum distance
      between all pairs of points).
    
    Default: \code{3}.
  }

  \item{\code{upperbound_scale_factor}}{
    The upper bound for the scale is determined
    as

    \code{upperbound_scale_factor} * (maximum distance
      between all pairs of points).
    
    Default: \code{3}.
  }
  \item{\code{upperbound_var_factor}}{ The upper bound for the
    variance and the nugget is determined 
    as \code{upperbound_var_factor} * var(\code{data})
    
    Default: \code{10}.
  }

  \item{\code{use_naturalscaling}}{
    logical. Only used if model is given in standard (simple) way.
    If \code{TRUE} then \emph{internally}, rescaled
    covariance functions will be used for which
    cov(1)\eqn{\approx}{~=}0.05.
    \code{use_naturalscaling} has the advantage that \code{scale}
    and the form parameters of the model get \sQuote{orthogonal},
    but \code{use_naturalscaling} does not work for all models.
    
    Note that this argument does not influence
    the output of \command{\link{RFfit}}: the parameter vector
    returned by \command{\link{RFfit}} refers
    \emph{always} to the standard covariance model as given in
    \command{\link{RMmodel}}. (In contrast to \code{practicalrange}
    in \command{\link{RFoptions}}.)\cr
    Advantages if \code{use_naturalscaling=TRUE}:
    
    \itemize{
      \item \code{scale} and the shape parameter of a parameterised
      covariance model can be estimated better if they are estimated
      simultaneously.
      \item The estimated bounds calculated by means of
      \code{upperbound_scale_factor} and \code{lowerbound_scale_factor},
      etc. might be more realistic.
      \item in case of anisotropic models, the inverse of the elements
      of the anisotropy matrix should be in the above bounds.
    }
    Disadvantages if \code{use_naturalscaling=TRUE}:
    \itemize{
      \item For some covariance models with additional parameters, the
      rescaling factor has to be determined numerically.
      Then, more time is needed to perform \command{\link{RFfit}}.
      \item note the \code{use_naturalscaling} only affects simple
      models, no operators. Also functions that define a parameter of
      the model are not changed.      
    }
    
    Default: \code{FALSE}.
  }
  }
}


\seealso{
  \command{\link{RFoptions}},
   \command{\link{RFfit}},
   \command{\link{RFratiotest}},
   \command{\link{RFcrossvalidate}}
}

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

RFoptions()$fit
\dontshow{FinalizeExample()}}


\keyword{spatial}

