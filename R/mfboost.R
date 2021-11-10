#################################################################################
#' Model-based Boosting for Manifold Valued Object Data
#' 
#' An interface for model-based gradient boosting when response observations
#' are naturally represented as vectors, matrices, 
#' or smooth potentially multidimensional functions. 
#' The function is a wrapper for the function \code{\link[mboost]{mboost}} 
#' in the model-based boosting package \code{mboost} via its functional pendant 
#' \code{\link[FDboost]{FDboost}} from the \code{FDboost} package.
#' For manifold valued responses, like shapes, appropriate loss functions 
#' can be fitted using an \code{\link[mfFamily]{mfFamily}}.
#' 
#' @param formula a symbolic description of the model formula \code{y ~ ...} 
#' on covariate level, i.e. specified as if the response was scalar, 
#' where \code{y} refers to an object (in \code{data}) containing all response data.
#' @param obj.formula intrinsic model formula for the internal representation 
#' of the response. E.g., for a functional response \eqn{y_i(t)} typically use 
#' \code{value ~ bbs(t)} to obtain smooth effects over \eqn{t}, where \code{value}
#' refers to the response evaluations (in \code{y}). See details for further 
#' information.
#' @param data a list containing the variables in the model. 
#' The response should be either provided as \code{tbl_cube} (regular case) or
#' as list of \code{tbl_cube}s or \code{data.frame}s (irregular case). See details.
#' @param family an \code{\link[mfFamily]{mfFamily}} object.
#' @param ... additional arguments passed to \code{\link[FDboost]{FDboost}}/\code{\link[mboost]{mboost}}.
#' 
#' @details 
#' While the response observations \eqn{y_i} and the corresponding predictions 
#' \eqn{\mu_i} might live on a non-linear manifold \eqn{M}, they are modeled with
#' an additive predictor \eqn{\eta_i} living in a linear space by
#' \deqn{\mu_i = g_p(\eta_i) = g_p(\sum_j h_j(x_i))}
#' where \eqn{g_p} is a response function, which might depend on a pole \eqn{p \in M}
#' and is typically chosen as \eqn{g_p = Exp_p}, the manifold exponential function
#' at \eqn{p}. The additive predictor \eqn{\eta_i} is composed of partial effects
#' \eqn{h_j(x_i)} as provided by the R packages \code{mboost} and \code{FDboost}, 
#' but potentially constraining them to a respective linear subspace, e.g. 
#' corresponding to the tangent space at \eqn{p}.
#' 
#' For further details on available covariate effects see 
#' \code{\link[FDboost]{FDboost}} and the \code{\link[baselearners]{baselearners}}
#' help of \code{mboost}.
#' 
#' Computationally, it might make a huge difference whether response observations 
#' are measured on a common regular grid or on irregular individual grids. In the 
#' regular case, the linear array model can be utilized 
#' (Brockhaus et al. 2015, Currie et al. 2006) for the design
#' matrix and the regular structure might also be used to speed up pole and 
#' gradient computation.
#' This distinction is reflected - and controlled by - the data format and 
#' in particular the format the response is provided in. For a data set with 
#' \eqn{N} observations, \code{data} should be 
#' provided as list containing scalar covariate vectors (of length \eqn{N}). 
#' The response should be contained as follows:
#' \itemize{
#'  \item in the regular case, the response should be 
#'  a \code{\link[tbl_cube]{tbl_cube}} with response 
#'  values as measures and remaining variables contained in
#'  \code{obj.formula} as dimensions, listed according to the covariates.
#'  \item in the irregular case, the response should be a \code{data.frame} 
#'  (or list) with the response variables appearing in the \code{obj.formula}. 
#'  \item for backward compatibility, data can also be in the format of 
#'  \code{FDboost}, i.e., \code{data} is again a list
#'  and the response is provided in separate list elements:
#'  \itemize{
#'   \item in the regular case, the response measurements are provided as matrix 
#'   with \eqn{N} rows corresponding to the observations and the columns 
#'   containing all measurements in long format and the \code{dim} and \code{obj.formula}
#'   variables as vectors along the columns of the matrix.
#'   \item in the irregular case, the columns of the last irregular option above
#'   are just separately added to the list.  
#'  } 
#' }
#' 
#' @seealso \code{\link{mfFamily}}
#' 
#' @return An object of class \code{mfboost} inheriting from \code{FDboost} and
#' \code{mboost}.
#'
#' Brockhaus, S., Scheipl, F., Hothorn, T. and Greven, S. (2015): 
#' The functional linear array model. Statistical Modelling, 15(3), 279-300. 
#' 
#' Currie, I.D., Durban, M. and Eilers P.H.C. (2006):  
#' Generalized linear array models with applications to multidimensional smoothing. 
#' Journal of the Royal Statistical Society, Series B-Statistical Methodology, 68(2), 259-280.
#'
#'
#' @export
#' @import FDboost mboost
#' @importFrom MASS Null
#' @importFrom dplyr bind_rows
#' @import formula.tools 
#' @example tests/mfboost_cells.R 

mfboost <- function(formula,               # response ~ xvars 
                    obj.formula = NULL,    # response.value^dim ~ blearner(arg) | id
                    data = NULL,           # list with response and covariates
                    family = Gaussian(),
                    ...) {                 # further arguments passed to mboost
  
# Interpret response data according to obj.formula ----------------------

  data_FDboost <- as_FD(data, formula = formula, obj.formula = obj.formula, ...)
  
  # Interpret formulae for FDboost --------------------------------------
  
  formulae_FDboost <- mfFormulae_for_FDboost(obj.formula, formula = formula, ...)
  
  # set no weights if nothing else is given
  ngradient_weights <- "equal"
  
# equip manifold geometry in family ------------------------
  
  if(inherits(family, "mfboost_family")){
    # already extract pole here, because initialization might remove it
    # pole_ <- family@mf$pole_
    
    # initialize response AND inner weights of geometry
      # load response
      family@mf$initialize(data = data_FDboost, formula = obj.formula)
    
    # if(all(sapply(pole_, is.null))) {
      ### create offset from response and obj.formula
      family@pole$fit(obj.formula = obj.formula, data = data_FDboost)
      p <- family@pole$predict()
      a <- try(p <- family@mf$structure(p))
      if(inherits(a, "try-error"))
        warning("Pole could not be structured.")
    # } 
  
    # make sure that pole_ really meets all requirements
    # i.e. it's also normalized etc.
    family@mf$pole_ <- family@mf$register(p)
    
    # specify weights for ngradient fit
    if(!all(sapply(family@mf$weights_, is.null))) {
      ngradient_weights <- family@mf$unstructure_weights(family@mf$weights_)
      if(!is.vector(ngradient_weights) | !is.numeric(ngradient_weights))
        stop("Only numeric vector weights can be specified for the gradient fit, sorry.")
    }
    
    # update timeformula according to family
    formulae_FDboost$timeformula <- family@update_formula(
      formulae_FDboost$timeformula, family@mf$pole_)
    data_FDboost <- c(data_FDboost, 
                      attr(formulae_FDboost$timeformula, "update_formula_vars"))
    attr(formulae_FDboost$timeformula, "update_formula_vars") <- NULL
  }
  
# fit FDboost model -------------------------------------------------------
  yvalname <- lhs.vars(formulae_FDboost$formula)
  ret <- FDboost(formula = formulae_FDboost$formula, 
                 timeformula = formulae_FDboost$timeformula, 
                 id = if(!is.matrix(data_FDboost[[yvalname]])) 
                   formulae_FDboost$id, 
                 numInt = ngradient_weights, data = data_FDboost, 
                 family = family,
                 offset = 0, check0 = FALSE, ...)
  
  ### assign new class (e.g., for specialized predictions)
  class(ret) <- c("mfboost", class(ret))
  # save the call
  ret$call <- match.call()
  ret$obj.formula <- obj.formula
  ret$formula <- formula
  
  # for compatibility with FDboost::predict.FDboost
  attr(ret$yind, "nameyind") <- names(ret$yind)[1]
  
  ret
}

# plot method for mfboost -------------------------------------------------

#' Plot effect estimates of mfboost model
#' 
#' Illustrate estimated effects by comparing predictions with corresponding 
#' reference poles.
#'
#' @param x object of class \code{mfboost} for plotting
#' @param which a subset of base-learners used for plotting. 
#' If which is given (as an integer vector or characters corresponding to base-learners) 
#' predictions only based on the corresponding learners are depicted. 
#' Per default all selected base-learners are plotted.
#' @param ids response ids to be plotted.
#' @param multiplier a number multiplied to the predictor before passing it to the
#' response function. Can be used for re-scaling / amplification.
#' @param ... other arguments passed to \code{mfGeometry} plot method.
#' 
#' @details 
#' Default plot function for \code{mfboost} objects. Example predictions are 
#' plotted using the plot method of \code{x$family@mf}, the geometry 
#' \code{mfGeom} class of the response.
#' 
#' 
#' @export
#'
#' @example tests/mfboost_cells.R
#' 
plot.mfboost <- function(x, which = NULL, ids = 1:6, multiplier = 1, ...) {
  if(is.null(x$family@mf$plot)) {
    warning("No plot function available for geometry.")
  } else {
    e <- environment(x$family@response)
    mf <- x$family@mf$clone(deep = TRUE)
    mf$y_ <- e$response_(
      multiplier * predict.mboost(x, which = which, type = "link"))
    mf$slice(ids)$plot(...)
  }
}


# predict method for mfboost --------------------------------------------

#' Predict with mfboost model
#' 
#' Get predictions of an \code{mfboost} model.
#'
#' @param object object of class \code{mfboost} for plotting
#' @param newdata an optional dataset with new data for prediction. Should be
#' of the same structure as the original dataset 
#' (also including a slot for the response!). The default NULL returns predictions
#' on the data used for fitting.
#' @param which a subset of base-learners to take into account for computing 
#' predictions or coefficients as in \code{predict.mboost}. \code{which = 0} will
#' predict the pole.
#' @param ... other arguments passed to \code{predict.mboost}.
#'
#' @return vector with predictions to plug into \code{newdata}.
#' @export
#' 
#' @example tests/mfboost_cells.R
predict.mfboost <- function(object, newdata = NULL, which = NULL, 
                            type = c("response", "link"), ...) {
  
  type <- match.arg(type)
  
  if(is.null(newdata)) {
    if(identical(which,0)) {
      pred <- rep(0, length(object$response))
      if(type == "link")
        return(pred) else
          return(object$family@response(pred))
    }
    
    return(predict.mboost(object, which = which, 
                           type = type, ...))
  } else {
    # Interpret new response data according to obj.formula ----------------------
    formula <- object$formula
    obj.formula <- object$obj.formula
    
    data_FDboost <- as_FD(newdata, formula = formula, obj.formula = obj.formula, ...)
    
    family <- object$family
    
    # equip manifold geometry in family ------------------------
    
    if(inherits(family, "mfboost_family")){
      
      family <- clone(family)
      
      # initialize response AND inner weights of geometry
      family@mf$initialize(data = data_FDboost, formula = obj.formula)
      mf <- family@mf
      
      ### create offset from response and obj.formula 
      p <- family@pole$predict(newdata = data_FDboost)
      p <- mf$structure(p)
      
      # make sure that pole_ really meets all requirements
      # i.e. it's also normalized etc.
      mf$pole_ <- mf$register(p)
      
      # # update response function to new geometry
      # response_ <- environment(family@response)$response_
      # if(is.function(response_)) environment(response_) <- environment()
      # environment(family@response) <- environment()
      
      # attach new tangent space constraint vectors
      data_FDboost <- c(data_FDboost, 
                        attr(family@update_formula(~0, pole_ = mf$pole_),
                             "update_formula_vars"))

    }
    
    # expand covariates (due to a bug in FDboost - should in fact be done there)
    nameid <- attr(object$id, "nameid")
    if(!is.null(nameid)) {
      excov <- which(lengths(data_FDboost) == length(unique(data_FDboost[[nameid]])))
      for(i in excov) {
        data_FDboost[[i]] <- data_FDboost[[i]][data_FDboost[[nameid]]]
      }
    }
    
    class(object) <- setdiff(class(object), "mfboost")
    
    if(identical(which, 0)) {
      pred <- rep(0, length(data_FDboost[[attr(object$yind, "nameyind")]]))
    } else {
      # arrange variable order as needed for FDboost
      nms <- names(object$yind)
      onms <- setdiff(names(data_FDboost), nms)
      data_FDboost <- data_FDboost[c(onms, nms)]
      
      pred <- predict(object, newdata = data_FDboost, to_FDboost = FALSE,
                      which = if(!is.null(which)) object$which(which), 
                      type = "link", ...)
    }
    if(type == "response") pred <- family@response(pred)
  }
  
  return(pred)  
}

