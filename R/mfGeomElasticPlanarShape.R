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
                                     # store also data.frame version of y for warping
                                     private$.y_dat <- data.frame(
                                       x = Re(private$.y_),
                                       y = Im(private$.y_),
                                       t = attr(private$.y_, "arg"))
                                     
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
                                   align = function(y_, y0_, closed = FALSE, eps = 0.01) {
                                     if(missing(y_)) {
                                       y_ <- self$register(
                                         private$.warp(y0_ = y0_, 
                                                       closed = closed, 
                                                       eps = eps))
                                       rot <- private$.innerprod(y_, y0_)
                                       return( rot * y_ / Mod(rot) )
                                     }
                                     
                                     self$register(
                                       private$.warp(y_, y0_, 
                                                     closed = closed, 
                                                     eps = eps))
                                     rot <- private$.innerprod(y_, y0_)
                                     y_ <- rot * y_ / Mod(rot)
                                   },
                                   #' @description Apply Log function of the sphere after rotation alignment
                                   #' @param method alternatives "simple" and "alternative" for the expression 
                                   #' used to compute the sphere Log-map. Passed to parent method.
                                   log = function(y_, y0_ = self$pole_, method = c("simple", "alternative")) {
                                     if(missing(y_))
                                       return(super$log(self$align(y0_ = y0_), 
                                                        y0_, method))
                                     super$log(
                                       private$.align(y_ = y_, y0_ = y0_)
                                       , y0_, method)
                                   },
                                   distance = function(y0_, y1_, squared = FALSE) {
                                     if(missing(y1_)) {
                                       y1_ <- self$align(y0_ = y0_)
                                     } else {
                                       y1_ <- self$align(y_ = y1_, y0_ = y0_)
                                     }
                                     ip <- private$.innerprod(y0_, y1_)
                                     if(Mod(ip-1) < 1e-15) return(0)
                                     if(Mod(ip+1) < 1e-15) return(pi)
                                     ret <- acos( Mod(ip) ) 
                                     if(squared)
                                       return(ret^2)
                                     ret
                                   }
                                 ),
                             private = list(
                               .y_dat = NULL,
                               .arg_ = NULL,
                               .warp = function(y_, y0_, closed = FALSE, eps = .01) {
                                 d0_ <- data.frame(
                                   x = Re(y0_), 
                                   y = Im(y0_),
                                   t = if(is.null(d0_, "arg")) private$.arg_ else
                                     is.null(d0_, "arg") )
                                 
                                 if(missing(y_)) {
                                   w <- elasdics::align_curves(
                                     d0_, private$.y_dat,
                                     eps = eps,
                                     closed = closed
                                   )
                                   # update internal parameterization
                                   private$.y_dat$t <- w$data_curve2_aligned$t_optim
                                 } else {
                                   d_ <- data.frame(
                                     x = Re(y_), 
                                     y = Im(y_),
                                     t = if(is.null(d_, "arg")) private$.arg_ else
                                       is.null(d_, "arg") )
                                   
                                   w <- elasdics::align_curves(
                                     d0_,d_,
                                     eps = eps,
                                     closed = closed
                                   )
                                 }
                                   
                                 structure(
                                   complex(
                                     re = approx(w$data_curve2_aligned$t_optim, 
                                                 Re(y_), 
                                                 xout = w$data_curve1$t)$y,
                                     im = approx(w$data_curve2_aligned$t_optim, 
                                                 Im(y_), 
                                                 xout = w$data_curve1$t)$y
                                   ),
                                   arg = w$data_curve1$t
                                 )}
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
