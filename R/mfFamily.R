

# mfPole class ------------------------------------------------------------

setOldClass("mfPole")

#' R6 Class for mfboost pole model object for use in \code{mfboost_family}
#'
#' A pole (similar to an offset) for a \code{mfboost} model is typically result
#' of some kind of model fit, which might be another simpler mfboost model
#' or some other kind of model procedure estimating the pole from the data.
#' For later model predictions on new data, it is important to store the 
#' not only the resulting pole, but the whole model behind it, to have access
#' to new pole predictions. The abstract R6 class \code{mfPole} implements
#' the basic structure needed to provide a pole for mfboost.
#'
#' @field model A fitted model object.
#'
# #' @name mfPole
# #' @rdname mfPole
#' @import R6
#' @export

mfPole <- R6Class("mfPole", 
                  public = list(
                    #' @description An initialization function determining the pole from the 
                    #' data. The function takes two arguments.
                    #' @param obj.formula Model formula for the pole model.
                    #' @param data Data in \code{FDboost} format.
                    fit = function(obj.formula, data) {
                      stop("No pole model fitting procedure determined.")
                    },
                    model = NULL, # model object
                    #' @description A function taking \code{newdata} in the same format as \code{data}
                    #' and generating a prediction as a simple vector fitting to the data.
                    #' @param newdata in \code{FDboost} format. The default \code{NULL}
                    #' specifies the original data.
                    predict = function(newdata = NULL) {
                      predict(model, newdata = newdata)
                    }
                  )    
) #mfGeometry



# mfboost_family class  ---------------------------------------------------------

#' S4 class for a gradient boosting family with geometry slot.
#'
#' The class \code{mfboost_family} inherits from \code{\link[boost_family]{boost_family}},
#' the \code{mboost} family class. In particular, it's used to define the negative
#' gradient fitted throughout the boosting algorithm, which might now also depend 
#' on the underlying geometry of the response.
#'
#' @slot mf an \code{\link[mfGeometry]{mfGeometry}} object specifying the 
#' underlying geometry.
#' @slot pole an \code{\link[mfPole]{mfPole}} object containing a pole model.
#' @slot ngradient,loss,risk,response,nuisance functions defining the 
#' \code{mboost::boost_family} object.
#' @slot weights,name \code{boost_family} slots indicating the weights allowed and
#' the name of the loss function.
#' @slot offset,fW,rclass,check_y,type other slots of \code{boost_family} objects
#' that typically do not play a role in \code{mfboost_family} objects.
#' 
#' @name mfboost_family-class
#' @rdname mfboost_family-class
#' @export
setClass("mfboost_family", contains = "boost_family",
         representation = representation(
           mf = "ANY",
           pole = "ANY",
           update_formula = "function"
         )) 


# constructor for mfboost_family S4 class --------------------------------

#' Manifold Gradient Boosting Families
#' 
#' \code{mfboost_family} objects extend \code{boost_family} objects to contain,
#' i.a., a geometry slot, which is used to transfer geometric structure information
#' of the response into the functions defined in a \code{boost_family}.
#' It is used to specify loss, risk and other functions for fitting regression 
#' problems with manifold valued response with \code{mfboost}
#' in a flexible and modular way.
#' 
#' @param mf an \code{\link[mfGeometry]{mfGeometry}} object specifying the 
#' underlying geometry.
#' @param initialize a function with arguments \code{formula} (in obj.formula 
#' format) and \code{data} (in a format used by FDboost) initializing \code{mf} 
#' by setting \code{mf$y_} and \code{mf$weights_}. 
#' @param pole a function with arguments \code{formula} and \code{data} 
#' (in a format used by FDboost) returning the pole in the format required by \code{mf}.
#' @param update_formula a function with arguments \code{formula} and \code{pole_}
#' updating the formula in dependence of the \code{pole_}, typically by applying a
#' function to the right side.
#' @param ngradient,loss,risk,response,nuisance functions defining the 
#' \code{mboost::boost_family} object.
#' @param weights,name \code{boost_family} character strings indicating the weights 
#' allowed and the name of the loss function.
#' @param offset,fW,rclass,check_y,type other slots of \code{boost_family} objects
#' that typically do not play a role in \code{mfboost_family} objects.
#' 
#' @export
#' @name mfFamily
#' @rdname mfFamily
mfFamily <- function(mf, pole, ngradient, 
                     update_formula = function(timeformula, pole_) timeformula, 
                     loss = NULL, 
                     risk = NULL, response = function(f) NA, 
                     nuisance = function() return(NA), 
                     weights = c("any", "none", "zeroone", "case"), 
                     name = "user-specified mfFamily", 
                     offset = function(y, w) 0, 
                     fW = NULL) {
  # check geometry
  stopifnot(inherits(mf, "mfGeometry"))
  stopifnot(inherits(pole, "mfPole"))
  
  ret <- Family(ngradient = ngradient, loss = loss, risk = risk, 
         offset = offset, check_y = function(y) y, 
         weights = weights, nuisance = nuisance, 
         name = name, fW = fW, response = response)
  ret <- as(ret, c("mfboost_family"))
  ret@mf <- mf
  ret@pole <- pole
  ret@update_formula <- update_formula
  ret
}

#' (Deep) clone an S4 class object containing R6 slots
#'
#' @param object an S4 class object
#' @param deep should clone of R6 slots be deep
#'
#' @return cloned version of the object 
#' @export
#'
setGeneric("clone", function(object, deep = TRUE) {
  standardGeneric("clone")
})

#' (Deep) clone an S4 class object containing R6 slots
#'
#' @param object mfFamily object
#' @param deep should clone of R6 slots be deep 
#'
#' @return cloned object of class mfFamily
#' @export
#' @importFrom rlang env_clone
#'
setMethod("clone", signature(object = "mfboost_family"), 
          function(object, deep = TRUE) {
  old <- environment(object@ngradient)
  e <- rlang::env_clone(old)
  object@mf <- e$mf <- e$mf$clone(deep = deep)
  object@pole <- e$pole <- e$pole$clone(deep = deep)
  # change all environments from old to e
  for(i in ls(envir = e)) {
    this <- environment(e[[i]])
    if(!is.null(this) & identical(this, old))
      environment(e[[i]]) <- e
  }
  # change all environments of all function fields
  slts <- getSlots("mfboost_family")
  for(i in names(slts)[slts == "function"])
    environment(slot(object, i)) <- e
  
  object
})


# constructor for Riemmanian L2 boosting ---------------------------------
#' 
#' @param mf The response geometry supplied as an \code{mfGeometry} object.
#' @param pole.type one of "RiemannL2"(default) and "Gaussian"
#' @param pole.control a list of parameters controlling the \code{mboost} algorithm. 
#' For more details see \code{\link{boost_control}}.
#' 
#' @export
#' @name RiemannL2
#' @rdname mfFamily
RiemannL2 <- function(mf, pole.type = c("RiemannL2", "Gaussian"), pole.control = boost_control()) {
  stopifnot(inherits(mf, "mfGeometry"))
  pole.type <- match.arg(pole.type)

  response_ <- function(f){
    f_ <- mf$structure(f)
    f_ <- mf$register_v(f_, mf$pole_)
    f_resp <- mf$exp(f_, mf$pole_)
  }
  
  response <- function(f){
    mf$unstructure(
      response_(f)
    )
  }
  
  loss <- function(y, f, w = 1) {
    mu_ <- response_(f)
    eps_ <- mf$log(mf$y_,mu_)
    w * mf$unstructure(eps_)^2
  }
  
  risk <- function(y, f, w = 1) {
    sum(loss(y, f, w))
  }
  
  ngradient <- function(y, f, w = 1) {
    f_resp <- response_(f)
    # Tangent vector in tangent space of respective y[i,]
    eps_ <- mf$log(mf$y_, f_resp)
    # Parallel transport tangent vectors to TpM
    eps_ <- mf$transport(eps_, f_resp, mf$pole_)
    # Return eps
    mf$unstructure( eps_ )
  }
  
  pole <- R6Class("mfPoleRecursive", inherit = mfPole,
                  public = list(
                    fit = function(obj.formula, data) {
                      new_fam <- switch(pole.type, 
                                        Gaussian = Gaussian(),
                                        RiemannL2 = RiemannL2(mf, pole.type = "Gaussian", pole.control = pole.control)
                                        )
                      self$model <- mfboost( formula = ~ 1, obj.formula = obj.formula, 
                                       data = data, family = new_fam, control = pole.control )
                    },
                    predict = function(newdata = NULL) 
                        predict(self$model, type = "response", newdata = newdata)
                      ))$new()
  
  update_formula <- function(timeformula, pole_) {
    
    normal_vecs <- mf$get_normal(y0_ = pole_, weighted = TRUE)
    
    if(!is.null(normal_vecs)) {
      normal_vecs <- as.data.frame(normal_vecs)
      names(normal_vecs) <- paste0("update_formula_var", 1:ncol(normal_vecs))
      
      ## transform obj.formula into appropriate tangent space
      traformula <- as.formula(paste0("~ . %-% bols(", 
                                      paste(names(normal_vecs), collapse = ","), 
                                      ", intercept = FALSE)"))
      timeformula <- update(timeformula, traformula )
      attr(timeformula, "update_formula_vars") <- normal_vecs
    }
    
    timeformula
  }
  
  # update_formula <- function(timeformula, pole_) timeformula
  
  mfFamily(mf, pole = pole, ngradient = ngradient, update_formula = update_formula,
           loss = loss, response = response, 
           weights = "any", risk = risk, 
           name = "Riemmanian L2-Boosting")
}


# Planar shape regression family ------------------------------------------

#' @param weight_fun a function producing inner product weights 
#' taking the arguments \code{arg} (vector of arguments of the function) and 
#' \code{range} (range of the arguments). Passed to \code{mf$initialize}.
#' @param arg_range vector of length 2 specifying the \code{range} argument of 
#' the \code{weight_fun}. The default \code{NULL} will take the minimum and maximum 
#' of the supplied \code{arg}.
#' 
#' @export
#' @name PlanarShapeL2
#' @rdname mfFamily
PlanarShapeL2 <- function(pole.type = "RiemannL2", pole.control = boost_control(), 
                          weight_fun = NULL, arg_range = NULL) {
  mf <- mfGeomProduct$new(
    mfGeom_default = mfGeomPlanarShape$new(weight_fun = weight_fun, arg_range = arg_range))
  
  RiemannL2(mf = mf, pole.type = pole.type, pole.control = pole.control)
}


# Planar size and shape regression family ------------------------------------------

#' 
#' @export
#' @name PlanarSizeShapeL2
#' @rdname mfFamily
PlanarSizeShapeL2 <- function(pole.type = "RiemannL2", pole.control = boost_control(), 
                              weight_fun = NULL, arg_range = NULL) {
  mf <- mfGeomProduct$new(
    mfGeom_default = mfGeomPlanarSizeShape$new(weight_fun = weight_fun, arg_range = arg_range))
  
  RiemannL2(mf = mf, pole.type = pole.type, pole.control = pole.control)
}

# Riemannian L2 Regression Euclidean case ----------------------------------

#' 
#' @export
#' @name EuclideanL2
#' @rdname mfFamily
EuclideanL2 <- function(cyclic = FALSE,
                        arg.grid.len = 5000, weights = weights, 
                        arg.range = NULL,
                        smoothed.cov = NULL, cov.k = 10) {
  mf <- mfGeomEuclidean_irreg$new()
  
  # TODO: replace by proper Euclidean mean function
  pole <- function(formula, data) 
    planarshape_full_proc(formula = formula, dim = dim, id = id, data = data,
                          cyclic = cyclic, smoothed.cov = smoothed.cov,
                          cov.k = cov.k, arg.grid.len = arg.grid.len,
                          mf = mf, weights = weights,
                          arg.range = arg.range,
                          mfboost.return = TRUE)
  warning("Currently, Euclidean family is just for test purposes and lacks a proper pole function.")
  RiemannL2(mf, pole)
}
