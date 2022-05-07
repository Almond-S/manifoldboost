
#' Factorize tensor product model with object response
#' 
#' Factorize an mfboost tensorporduct model into the response and covariate parts for 
#' effect visualization.
#'
#' @param x a model object of class mfboost.
#' @param blwise logical, should the factorization be carried out base-learner-wise (TRUE, default)
#' or for the whole model simultaneously.
#' @param ... other arguments passed to methods (no arguments so far).
#'
#' @return a list of two mboost models of class mboost_fac containing basisfunctions
#' for response and covariates, respectively, as base-learners.
#' @export
#' 
#' @name factorize
#' @aliases factorize.FDboost
#' 
#' @example tests/mfboost_cells.R 
#'
factorize.mfboost <- function(x, blwise = TRUE, newdata = NULL, newformula = NULL, newobj.formula = NULL, ...) {
  
  wghts <- NULL
  if(!is.null(newdata)) {
    # get new weights from new geometry 
    newformula <- if(is.null(newformula)) x$formula
    # newobj.formula <- if(is.null(newobj.formula)) x$obj.formula
    
    newdata <- as_FD(newdata, 
          formula = newformula, 
          obj.formula = newobj.formula)
    
    mf <- x$family@mf$clone(deep = TRUE)
    if(is.null(mf$formula) | !is.null(newobj.formula)) 
      mf$initialize(data = newdata, formula = newobj.formula) else 
        mf$initialize(data = newdata)
    
    wghts <- mf$unstructure_weights(mf$weights_)
    # predict new pole
    mf$pole_ <- mf$structure(x$family@pole$predict(newdata = newdata))
    # update family
    x$family <- RiemannL2(mf)
    
    # attach new tangent space constraint vectors
    newdata <- c(newdata, 
                      attr(x$family@update_formula(~0, pole_ = mf$pole_),
                           "update_formula_vars"))
  }
  
  ret <- factorize.FDboost(x = x, 
                           newdata = newdata, 
                           newweights = wghts, 
                           blwise = blwise, ...)
  ydim <- ret$resp$ydim
  if(length(ydim)==2) { # (i.e. the regular case)
    # adjust geometry in family
    ret$resp$family@mf <- ret$resp$family@mf$clone(deep = TRUE)$slice(seq_len(ydim[1]))
    old <- environment(ret$resp$family@response)
    e <- rlang::env_clone(old)
    e$mf <- ret$resp$family@mf
    # change all environments from old to e
    for(i in ls(envir = e)) {
      this <- environment(e[[i]])
      if(!is.null(this) & identical(this, old))
        environment(e[[i]]) <- e
    }
    environment(ret$resp$family@response) <- e
  }
  class(ret$resp)[1] <- "mfboost_fac"
  # missing information for prediction
  ret$resp$formula <- x$formula
  ret$resp$obj.formula <- x$obj.formula
  ret
}

#' @importFrom methods setOldClass
#' @exportClass mfboost_fac

setOldClass("mfboost_fac")

#' `mfboost_fac` S3 class for factorized mfboost model component
#'
#' @description Model factorization with `factorize()` decomposes an 
#' `mfboost` model into an object of class `mfboost_fac` for the 
#' response and and object of class `FDboost_fac` for the covariate predictor. 
#' The first is essentially an `mfboost` object and the second an `mboost` object, 
#' however, 
#' in a 'read-only' mode and slightly adjusted methods (method defaults).
#'
#' @name mfboost_fac-class
#' @seealso [factorize(), factorize.mfboost(), factorize.FDboost()]
NULL

#' @export
plot.mfboost_fac <- function(x, multiplier = 1, which = NULL, ids = 1, main = NULL, ...) {
  # determine which bls to plot
  w <- x$which(which, usedonly = TRUE)
  if(any(is.na(w)))
    stop(paste("Don't know which base-learner is meant by:", 
               which[which.min(is.na(w))]))
  # determine plot title
  if(is.null(main)) {
    nm <- names(x$baselearner)[w]
    main <- paste("Response ->", regmatches(nm, regexpr("\\[.*\\]", nm)))
    nm <- names(x$family@mf$y_)[ids]
    main <- lapply(main, paste, nm, sep = "\n")
  }
  
  for(i in seq_along(w)) {
    plot.mfboost(x, which = w[i], main = main[[i]], ids = ids, 
                 multiplier = multiplier, ...)
  }
}
