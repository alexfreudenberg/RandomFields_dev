
## ja nicht RFopt uebergeben!! Da RFoptions() neu gesetzt wird,
## aber nicht RFopt, so dass es zu Mehrfachaufrufen kommen kann

Help <- function(which, ...) {
  w <- paste0("help_", which)
  RFopt <- getRFoptions(getoptions_=c("basic", "messages"))
  if (RFopt$basic$printlevel == 0) return()
  opt <- RFopt$messages[[w]]
#  Print(RFopt, opt, w)
  if (RFopt$basic$helpinfo && opt) {
    value <- 0
    txt <- switch(which,
      "addNA" = "Currently, 'NA's are added to any linear part of the model definition. This behaviour can be changed by 'RFoptions(addNA=FALSE)'. See the man pages for further options.",
      "colour_palette" = "Better install one of the packages 'colorspace' or 'RColorBrewer'.",
      "help" = "Hints can be generally suppressed by setting 'RFoptions(helpinfo=FALSE)'. Any messages, including hints and notes, can be suppressed by 'RFoptions(printlevel=0)'.",
      "mle" = "The MLE is extracted from a list of results. Use argument 'method' in any function of the package 'RandomFieldsLight' to extract different results, see ?RFfit for possible values of 'method'. Set RFoptions(helpinfo=FALSE) to avoid this message.",
      "newstyle"="Currently, an S4 object of class 'RFsp' is returned that allows an easy handling of the result by the user, e.g. for plotting;\nProgrammers may prefer a bare, but faster array format and set 'RFoptions(spConform=FALSE)'.",
      "normal_mode"=paste0("The modus_operandi='", MODE_NAMES[normal + 1],
                           "' is save, but slow. If you like the MLE running\n",
                           "faster (at the price of being less precise!) ",
                           "choose modus='", MODE_NAMES[easygoing + 1],
                           "' or\neven modus='", MODE_NAMES[sloppy + 1], "'."),
      "onlyvar"="Only the variance has to be estimated (except for some parameter in the linear model part). Note that 'RFlikelihood' does already this job and is much simpler.",
      stop("Bug", CONTACT)
      )
    m <- list()
    m[[paste0("messages.", w)]] <- value
    do.call(setRFoptions, m)
    message("Hint: ", txt, " [This hint appears only once per session and can be fully suppressed by 'RFoptions(", w, "=FALSE)'.]\n")
    Help("help")
  }
}

Note <- function(which, ...) {
  w <- paste0("note_", which)
  RFopt <- getRFoptions(getoptions_=c("basic", "messages"))
  if (RFopt$basic$printlevel == 0) return()
  opt <- RFopt$messages[[w]]
#  Print(RFopt, opt, w)
  if (opt) {
    value <- 2 * as.integer(opt / 2)
    txt <- switch(which,
    "ambiguous" = "Interpretation as degenerated grid. Better give 'grid' explicitely.",
    "coordinates" = {
      f <- function(cur, units) {
         if (opt == 3)
           paste0("Earth coordinates are detected.\n",
                  "If this is not what you wish, change option 'coord_system'",
                  "\nto RFoptions(coord_system = \"cartesian\").\n",
                  "Note further that angles in R.cos, R.sin, R.tan, RMangle",
                  " are now expected\nin DEGREE and ",
                  "R.acos, R.asin, R.atan, R.atan2 return results in degree.\n",
                  "\nCurrent units acore ",
                  if (cur=="") "not given." else paste0("'", cur, "'."),
                  "\nIn rare cases earth coordinates will be transformed ",
                  "within submodels.",
                  "\nThen it will be transformed into units of '", units,
                  "'. In particular, the\nvalues of all scale parameters of ",
                  "any of these submodels, defined in R^3,\nare ",
                  "understood in units of '", units,
                  "'. Change option 'units' if necessary.\n")
         else paste0("Earth coordinates detected. In case it is necessary ",
                     "they will be transformed into units of ",  units, ".")}
      f(...)},
    "detection" = {paste0(...)},
    "aspect_ratio" = "The graphical device does not seem to be a standard screen device. Hence the\naspect ratio might not be correct.",
    "no_fit" = if (...) "No genuine variable has to be estimated (except possibly some parameter in the linear model part). Note that 'RFlikelihood' does already the job." else "No genuine variable has to be estimated. You should use 'RFlikelihood' instead.", ## ...=n.covariate
    stop("Bug", CONTACT)
    )
    m <- list()
    m[[paste0("messages.", w)]] <- value
    do.call(setRFoptions, m)

    ## opt = 0: keine Meldung
    ## opt = 1: "only once meldung"
    ## opt = 2: Meldung, aber ohne "[This message ...]"
    ## opt = 3: Meldung mit "[This message]"
    
    message("Note: ", txt,
            if (opt == 1) " [This note appears usually only once per session and can be fully "
            else if (opt == 3) " [This message can be",
            if (opt %% 2 == 1)
              paste0(" suppressed by 'RFoptions(", w, "=FALSE)'.]"), "\n")
  }
}

Warning <- function(which, ...) {
  w <- paste0("warn_", which)
  RFopt <- getRFoptions(getoptions_=c("messages"))
  opt <- RFopt[[w]]
 # Print(RFopt, opt, w)
  if (opt) {
    value <- 0
    txt <- switch(which,
       "ambiguous" = "Ambiguous interpretation of the coordinates. Better give the logical parameter 'grid=TRUE' explicitely.",
       "on_grid" = paste("coordinates", if (...) " do not", # ... = grid
                         " seem to be on a grid, but grid = ", ...),
       "oldstyle" = {
         value <- 1
         opt <- 2
          paste0("Your are using an obsolete function definition. Use '", ...,
                "' instead.")
       },
       stop("Bug", CONTACT)
       )
    m <- list()
    m[[paste0("messages.", w)]] <- value
    do.call(setRFoptions, m)
    warning(txt,
            if (opt == 1) paste0(" [This warning appears only once per session and can be fully suppressed by 'RFoptions(", w, "=FALSE)'.]"))
  }
}



warn.seed.not.na <- function(RFoptOld, oldstyle=FALSE) {
  RFopt <- RFoptOld[[2]]
  basic <- RFopt$basic
  if (!is.na(basic$seed)){  
    o.seed <- RFoptOld[[1]]$basic$seed
    allequal <- all.equal(o.seed, basic$seed)
    allequal <- is.logical(allequal) && allequal
    if (basic$printlevel >= PL_IMPORTANT &&
        (is.null(o.seed) || (!is.na(o.seed) && allequal)
         )
        ) {
      warn_seed <- RFopt$messages$warn_seed
      if (warn_seed > 0) {
        if (warn_seed > 1) {
          setRFoptions(messages.warn_seed = warn_seed - 1)
          txt <- "\nSet 'RFoptions(seed=NA)' to make the seed arbitrary."
        } else txt <- ""
        message("NOTE: simulation is performed with fixed random seed ",
                basic$seed, ".", txt)
      }
    }
    if (oldstyle) {
      warning("Fixed seeds in the old style result in a different behaviour of R itself! While in the old style, the state of .Random.seed is influenced for fixed seed, it is not in the new style. The user is urged to switch to the new style.")
    }
  }
}
