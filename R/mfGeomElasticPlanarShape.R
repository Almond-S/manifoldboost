# Shape space of elastic plane curves (single vector) --------------------------

#' Plane Shape Space Geometry with warping alignment
#' 
#' @description Geometry of the elastic shape space identifying planar shapes
#' with a centered and scaled complex vector on the complex sphere / complex
#' projective space.
#' 
#' @param y a numeric vector representing an element of the manifold in an 
#' unstructured way. The position of elements and dimensions is obtained from 
#' the initialization of the geometry.
#' @param v a numeric vector representing an tangent vector in an 
#' unstructured way. Complex interpretation as above. 
#' @param y_ a complex vector on the sphere, orthogonal to the constant.
#' @param y0_ a complex vector on the sphere, orthogonal to the constant.
#' @param y1_ a complex vector on the sphere, orthogonal to the constant.
#' @param v_ a complex tangent vector.
#' @param v0_ a complex tangent vector.
#' @param weights_ numeric vector of inner product weights matching the \code{y_}
#' as a complex vector.
#' @param weights numeric vector of inner product weights matching the \code{y}
#' as a numeric vector. \code{weights_} are duplicated for values of both dimensions.
#'
#' @import elasdics
#' @export
mfGeomWarpPlanarShape <- R6Class("mfGeomWarpPlanarShape", inherit = mfGeomPlanarShape, 
                                 public = list(
                                   #' @description Initialize planar shape geometry with warping 
                                   #' for data given in (long) FDboost format. 
                                   #' Same as in mfGeomPlanarShape but storing also function arguments. 
                                   #' 
                                   #' @param data data in the format used in \code{\link{FDboost}}.
                                   #' @param formula formula describing the internal structure of the data, 
                                   #' intended for the \code{obj.formula} of \code{mfboost} 
                                   #' (interpreted by \code{mfInterpret_objformula}).
                                   #' @param weight_fun a function computing inner product weights for the geometry
                                   #' taking two arguments: \code{arg}, the argument of the functional data 
                                   #' specified in \code{formula}, \code{range} the range of \code{arg}.
                                   #' @param arg_range the \code{range} supplied to the \code{weight_fun}.
                                   #' 
                                   initialize = function(data, formula, 
                                                         weight_fun = NULL, 
                                                         arg_range = NULL) {
                                     # dataless initializations
                                     if(!is.null(weight_fun)) {
                                       stopifnot(is.function(weight_fun))
                                       private$weight_fun <- weight_fun
                                     }
                                     
                                     if(!is.null(arg_range)) {
                                       private$arg_range <- arg_range
                                     }
                                     
                                     # with no data stop here !!!!!!!!!!!!!!!!!!
                                     if(missing(data)) {
                                       return(invisible(self))
                                     } 
                                     
                                     v <- mfInterpret_objformula(formula)
                                     
                                     # check dimension
                                     if(length(unique(data[[v$dim]]))!=2)
                                       stop(paste0("Length of dim variable '", 
                                                   v$dim,"' must be 2 to interpret it as complex."))
                                     n <- length(data[[v$value]])
                                     if(n%%2 != 0)
                                       stop("Number of curve measurements has to be even for complex interpretion.")
                                     if(length(data[[v$dim]]) != n | length(data[[v$arg]]) != n)
                                       stop("value, dim and arg variable have to be of equal length.")
                                     
                                     # order data with respect to arg and dim
                                     structure_index <- order(data[[v$dim]], data[[v$arg]])
                                     private$structure_index <- list(real = head(structure_index,n/2),
                                                                     imaginary = tail(structure_index, n/2)) 
                                     
                                     # structure y (internally using $structure_index just defined above)
                                     private$.y_ <- self$structure(data[[v$value]])
                                     # store argument and add as attribute
                                     attr(private$.y_, "arg") <- private$.arg_ <-
                                       data[[v$arg]][private$structure_index$real]
                                     
                                     # set up index for $unstructure
                                     index <- self$structure(seq_along(data[[v$value]]))
                                     private$unstructure_index <- order(private$.unstructure(index))
                                     
                                     # compute weights
                                     if(!is.null(private$weight_fun)) {
                                       private$.weights_ <- private$weight_fun(
                                         arg = data[[v$arg]][private$structure_index$real], 
                                         private$arg_range )
                                     }
                                     
                                     # registration already depends on weights
                                     self$y_ <- self$register(self$y_)
                                     stopifnot(length(private$.weights_) == length(private$.y_) | 
                                                 is.null(private$.weights_))
                                     
                                     invisible(self)
                                   },
                                   align = function(y_, y0_) {
                                     missy_ <- missing(y_)
                                     missy0_ <- missing(y0_)
                                     
                                     if(missy_) y_ <- private$.y_
                                     if(missy0_) y0_ <- private$.y_
                                     
                                     rot <- private$.innerprod(y_, y0_)
                                     y_ <- rot * y_ / Mod(rot)
                                     warp <- private$.warp(y_, y0_)
                                     
                                     if(missy_) {
                                       attr(private$.y_, "arg") <- warp$y_arg_opt
                                       private$.y_ <- self$register(private$.y_)
                                     }
                                     
                                     y_ <- warp$new_y_
                                     attr(y_, "arg") <- warp$new_y_arg
                                     
                                     self$register(y_)
                                   }
                                 ),
                             private = list(
                               .arg_ = NULL,
                               .warp = function(y_, y0_, closed = FALSE, eps = .01) {
                                 d_ <- data.frame(x = Re(y0_), y = Im(y0_))
                                 if(!is.null(attr(d_, "arg")))
                                   d_$t <- attr(d_, "arg")
                                 d0_ <- data.frame(x = Re(y_), y = Im(y_))
                                 if(!is.null(attr(d0_, "arg")))
                                   d0_$t <- attr(d0_, "arg")
                                 
                                 w <- elasdics::align_curves(
                                   d_,d0_,
                                   eps = eps,
                                   closed = closed
                                 )
                                 
                                 list(
                                   new_y_ = complex(
                                     re = approx(w$data_curve2_aligned$t_optim, 
                                                 Re(y_), 
                                                 xout = w$data_curve1$t)$y,
                                     im = approx(w$data_curve2_aligned$t_optim, 
                                                 Im(y_), 
                                                 xout = w$data_curve1$t)$y
                                   ),
                                   y_arg_opt = w$data_curve2_aligned$t_optim,
                                   new_y_arg = w$data_curve1$t
                                 )
                               },
                               .distance = function(y0_, y1_, squared = FALSE) {
                                 y1_ <- self$align(y1_, y0_)
                                 ip <- private$.innerprod(y0_, y1_)
                                 if(Mod(ip-1) < 1e-15) return(0)
                                 if(Mod(ip+1) < 1e-15) return(pi)
                                 ret <- acos( Mod(ip) ) 
                                 if(squared)
                                   return(ret^2)
                                 ret
                               }
                             ))

# warping aligned planar shape regression family -------------------------------

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
WarpPlanarShapeL2 <- function(pole.type = "RiemannL2", pole.control = boost_control(), 
                          weight_fun = equal_weights, arg_range = NULL) {
  mf <- mfGeomProduct$new(
    mfGeom_default = mfGeomWarpPlanarShape$new(weight_fun = weight_fun, arg_range = arg_range))
  
  RiemannL2(mf = mf, pole.type = pole.type, pole.control = pole.control)
}
