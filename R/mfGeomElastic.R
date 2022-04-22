


# Square-root-velocity geometry of plane curves (single curve) -----------------

#' Square-root-velocity (SRV) Geometry
#' 
#' @description L2 geometry on SRV-transforms of curves
#' 
#' @param y,v numeric vectors of SRVs in an unstructured way. 
#' @param y_,y0_,y1_,v_,v0_ complex vectors of SRVs in a structured way.
#' @param weights_ numeric vector of inner product weights matching the \code{y_}
#' as a complex vector.
#' @param weights numeric vector of inner product weights matching the \code{y}
#' as a numeric vector. \code{weights_} are duplicated for values of both dimensions.
#'
#' @import elasdics
#' @export
mfGeomSRV_closed <- R6Class("mfGeomSRV", inherit = mfGeomEuclidean, 
                                           public = list(
                                             #' @description Initialize plane SRV geometry.
                                             #' 
                                             #' @param data data in the format used in \code{\link{FDboost}}.
                                             #' @param formula formula describing the internal structure of the data, 
                                             #' intended for the \code{obj.formula} of \code{mfboost} 
                                             #' (interpreted by \code{mfInterpret_objformula}). 
                                             #' If a \code{formula} is already stored in the \code{mfGeometry} object,
                                             #' this is taken as default.
                                             #' @param weight_fun a function computing inner product weights for the geometry
                                             #' taking two arguments: \code{arg}, the argument of the functional data 
                                             #' specified in \code{formula}, \code{range} the range of \code{arg}.
                                             #' @param arg_range the \code{range} supplied to the \code{weight_fun}.
                                             #' 
                                             initialize = function(data, formula, 
                                                                   weight_fun = trapez_weights, 
                                                                   arg_range = c(0,1)) {
                                               # dataless initializations
                                               if(!is.null(weight_fun)) {
                                                 stopifnot(is.function(weight_fun))
                                                 private$weight_fun <- weight_fun
                                               }
                                               
                                               if(!is.null(arg_range)) {
                                                 private$arg_range <- arg_range
                                               }
                                               
                                               if(missing(formula)) {
                                                 formula <- private$.formula
                                               } else {
                                                 stopifnot(inherits(formula, "formula"))
                                                 private$.formula <- formula
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
                                               private$.curve_ <- private$close(
                                                 self$structure(data[[v$value]]))
                                               # store argument and add as attribute
                                               private$.arg_ <- data[[v$arg]][private$structure_index$real]
                                               attr(private$.curve_, "arg") <- c(private$.arg_, 
                                                                                 1-private$.arg_[1])
                                               
                                               # store SRV trafo of curve 
                                               private$.y_ <- self$srv_trafo(private$.curve_)
                                               
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
                                               self$register()
                                               stopifnot(length(private$.weights_) == length(private$.y_) |
                                                           is.null(private$.weights_))
                                               
                                               invisible(self)
                                             },
                                             #' @description convert 2D vectors (given in long format with dimension 
                                             #' indicator) to complex
                                             structure = function(y) {
                                               if(is.null(private$structure_index)) 
                                                 stop("Please, initialize geometry before using $structure.")
                                               complex(
                                                 real = y[private$structure_index$real], 
                                                 imag = y[private$structure_index$imaginary])
                                             },
                                             #' @description convert complex back to long format numeric
                                             unstructure =  function(y_) {
                                               if(is.null(private$unstructure_index)) 
                                                 stop("Please, initialize geometry before using $unstructure.")
                                               private$.unstructure(y_)[private$unstructure_index]
                                             },
                                             #' @description double inner product weights for complex vectors 
                                             #' to appear in both dimensions in numeric long format
                                             unstructure_weights = function(weights_) {
                                               self$unstructure(complex(re = weights_, im = weights_))
                                             },
                                             align = function(y_, y0_ = self$pole_, closed = TRUE) {
                                               if(missing(y_))
                                                 y_ <- private$.y_
                                               if(is.null(attr(y_, "arg")))
                                                 attr(y_, "arg") <- private$.arg_
                                               
                                                 private$interpolate(y_, 
                                                                     arg = if(is.null(attr(y0_, "arg"))) 
                                                                       private$.arg_ else 
                                                                         attr(y0_, "arg"), 
                                                                     closed = closed)
                                             },
                                             #' @description Perform SRV-transform of a curve and its inverse.
                                             #' @param inverse logical, should the inverse SRV-transform be conducted.
                                             #' @param center logical, if \code{inverse = TRUE}, should the output be centered.
                                             srv_trafo = function(y_, inverse = FALSE, center = FALSE) {
                                               # if(inverse & center) 
                                               #   warning("Sorry, centering not implemented, yet.")
                                               arg <- attr(y_, "arg")
                                               if(is.null(arg))
                                                 arg <- private$.arg_
                                               long_arg <- length(arg) == length(y_) + 1 & 
                                                 arg[1] %% 1 == tail(arg, 1) %% 1
                                               if(inverse) {
                                                 ret <- elasdics::get_points_from_srv(data.frame(
                                                   t = if(long_arg) arg[-length(arg)] else 
                                                     arg,
                                                   x = Re(y_), 
                                                   y = Im(y_)))
                                               } else {
                                                 ret <- elasdics::get_srv_from_points(data.frame(
                                                   t = arg,
                                                   x = if(long_arg) 
                                                     private$close(Re(y_)) else Re(y_), 
                                                   y = if(long_arg) 
                                                     private$close(Im(y_)) else Im(y_)))
                                               }
                                               ret <- structure(complex(
                                                 re = ret$x, im = ret$y
                                               ), arg = attr(y_, "arg"))
                                               
                                               if(inverse & center) 
                                                 ret <- private$center(ret)
                                               
                                               ret
                                             },
                                             #' @description default plotting function for shapes of plane curves, plotting \code{y_}
                                             #' in front of \code{y0_} (after alignment).
                                             #' @param level one of 'curve' (default) and 'SRV' indicating
                                             #' whether curves are to be plotted or their SRV-transforms
                                             #' @param closed logical, should the plotted functions be closed?
                                             #' @param center logical, if \code{level='curve'}, should curves be centered?
                                             #' @param add_original logical, should original curve y_ be displayed besides its aligned version. 
                                             #' @param col,pch,type graphical parameters passed to \code{base::plot} referring to \code{y_}.
                                             #' @param ylab,xlab,xlim,ylim,xaxt,yaxt,asp graphical parameters passed to \code{base::plot} 
                                             #' with modified defaults.
                                             #' @param ... other arguments passed to \code{base::plot}.
                                             #' @param y0_par graphical parameters for \code{y0_}.
                                             #' @param seg_par graphical parameters for line segments connecting \code{y_} 
                                             #' and \code{y0_}.
                                             plot = function(y_, y0_ = self$pole_, 
                                                             level = c("curve", "SRV"), 
                                                             closed = TRUE,
                                                             center = TRUE,
                                                             ylab = NA, xlab = NA, 
                                                             col = "black",
                                                             pch = 19, asp = 1,
                                                             xlim, ylim,
                                                             yaxt='n',
                                                             xaxt='n',
                                                             type = "l",
                                                             y0_par = list(col = "darkgrey", type = type), 
                                                             seg_par = list(col = "grey"), ...) {
                                               level <- match.arg(level)
                                               stopifnot(is.logical(closed))
                                               
                                               missy_ <- missing(y_)
                                               
                                               if(!is.null(y0_)) {
                                                 if(level == "curve") {
                                                   y_ <- if(missy_) private$.curve_ else 
                                                     self$srv_trafo(y_, inverse = TRUE, center = center)
                                                   y0_ <- self$srv_trafo(y0_, inverse = TRUE, center = center)
                                                   
                                                   if(missy_&center)
                                                     y_ <- private$center(y_)
                                                 } else {
                                                   if(missy_) y_ <- private$.y_
                                                 }
                                                 
                                                 if(closed) {
                                                   y_ <- private$close(y_)
                                                   y0_ <- private$close(y0_)
                                                 }
                                                 
                                                 if(is.null(y0_par)) 
                                                   y0_par <- list()
                                                 if(is.list(y0_par)) {
                                                   .y0_par <- c(private$.y0_par, type = type)
                                                   nm <- setdiff(names(.y0_par), names(y0_par))
                                                   y0_par[nm] <- .y0_par[nm]
                                                   y0_par$x <- Re(y0_)
                                                   y0_par$y <- Im(y0_)
                                                 } else {
                                                   stop("y0_par has to be supplied as a list (or NULL).")
                                                 }
                                                 
                                                 if(!is.na(seg_par)) {
                                                   if(is.null(seg_par)) 
                                                     seg_par <- list()
                                                   if(is.list(seg_par)) {
                                                     nm <- setdiff(names(private$.seg_par), names(seg_par))
                                                     seg_par[nm] <- private$.seg_par[nm]
                                                     seg_par[c("x0", "y0", "x1", "y1")] <- 
                                                       list(x0 = Re(y0_), y0 = Im(y0_), 
                                                            x1 = Re(y_), y1 = Im(y_))
                                                   } else {
                                                     stop("seg_par has to be supplied as a list (or NULL).")
                                                   }
                                                 } 
                                               } else {
                                                 if(level == "curve") {
                                                   y_ <- if(missy_) 
                                                     private$.curve_ else 
                                                       self$srv_trafo(y_, inverse = TRUE, center = center)
                                                 } else {
                                                   if(missy_) y_ <- private$.y_
                                                 }
                                                 if(closed) {
                                                   y_ <- private$close(y_)
                                                 }
                                               }
                                               
                                               if(missing(xlim))  
                                                 xlim <- range(Re(c(y_, y0_)))
                                               if(missing(ylim))
                                                 ylim <- range(Im(c(y_, y0_)))
                                               
                                               plot(x = Re(y_), y = Im(y_), ylab = ylab, xlab = xlab, 
                                                    yaxt = yaxt, xaxt = xaxt, type = type,
                                                    pch = pch, col = if(is.null(y0_)) col,
                                                    xlim = xlim, ylim = ylim, asp = asp, ...)
                                               
                                               if(!is.null(y0_)) {
                                                 if(!is.na(seg_par[1])) 
                                                   do.call(segments, seg_par)
                                                 do.call(points, y0_par)
                                                 points(x = Re(y_), y = Im(y_), pch = pch, col = col, type = type, ...)
                                               }
                                             }
                                           ), active = list(
                                             y_ = function(value) {
                                               if(missing(value))
                                                 private$.y_ else {
                                                   private$.y_ <- self$validate(value)
                                                   private$.curve_ <- self$srv_trafo(value, inverse = TRUE)
                                                 } }
                                           ),
                                           private = list(
                                             .curve_ = NULL,
                                             .arg_ = NULL,
                                             close = function(x) x[c(seq_len(length(x)), 1)],
                                             center = function(curve_) {
                                               ONE <- rep(1, length(curve_) -1)
                                               curve_ <- curve_ - self$innerprod(ONE, head(curve_, -1)) / self$innerprod(ONE)
                                             },
                                             interpolate = function(y_, arg, closed = TRUE) {
                                               if(missing(y_))
                                                 curve_ <- private$.curve_ else
                                                 curve_ <- self$srv_trafo(y_, inverse = TRUE, center = FALSE)
                                               
                                               curve_ <- structure(
                                                 complex(
                                                   re = private$approx(attr(curve_, "arg"), Re(curve_), xout = arg, closed = closed)$y,
                                                   im = private$approx(attr(curve_, "arg"), Im(curve_), xout = arg, closed = closed)$y
                                                 ),
                                                 arg = arg
                                               )
                                               self$srv_trafo(curve_)
                                             },
                                             approx = function(x, y = NULL, xout, closed = TRUE, xleft = 0, xright = 1, ...) {
                                               if(closed) {
                                                 x <- (x - xleft) %% (xright - xleft) + xleft
                                                 xout <- (xout - xleft) %% (xright - xleft) + xleft
                                                 xorder <- order(x)
                                                 x <- x[xorder]
                                                 y <- y[xorder]
                                                 xrange <- c(x[1], x[length(x)])
                                                 where <- range(findInterval(xout, xrange))
                                                 if(any(where != 1)) {
                                                   y <- c(if(0 %in% where) tail(y, 1),
                                                          y,
                                                          if(2 %in% where) head(y, 1))
                                                   x <- c(if(0 %in% where) xleft + tail(x, 1) -xright,
                                                          x,
                                                          if(2 %in% where) xright + head(x, 1)-xleft)
                                                 }
                                               }
                                               approx(x = x, y = y, xout = xout, ...)
                                             },
                                             .y0_par = list(col = "darkgrey"),
                                             .seg_par = list(col = "grey"),
                                             .unstructure = function(y_) c(Re(y_), Im(y_)),
                                             # needed for structuring and unstructuring
                                             structure_index = NULL,
                                             unstructure_index = NULL,
                                             weight_fun = NULL,
                                             arg_range = NULL
                                           ))

# Family for regression on SRV level -------------------------------

#' @param weight_fun a function producing inner product weights 
#' taking the arguments \code{arg} (vector of arguments of the function) and 
#' \code{range} (range of the arguments). Passed to \code{mf$initialize}.
#' @param arg_range vector of length 2 specifying the \code{range} argument of 
#' the \code{weight_fun}. The default \code{NULL} will take the minimum and maximum 
#' of the supplied \code{arg}.
#' 
#' @export
#' @name SquareRootVelocityL2
#' @rdname mfFamily
SquareRootVelocityL2 <- function(
  weight_fun = trapez_weights, arg_range = c(0,1), closed = TRUE) {
  if(!closed) stop("Sorry, only closed case implemented so far.")
  mf <- mfGeomProduct$new(
    mfGeom_default = mfGeomSRV_closed$new(weight_fun = weight_fun, 
                                                         arg_range = arg_range))
  
  RiemannL2(mf = mf, pole = mfPoleZero$new())
}




# Shape space of elastic plane curves (single vector) --------------------------

#' Elastic Plane Shape Space Geometry
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
mfGeomElasticPlanarShape_closed <- R6Class("mfGeomElasticPlanarShape_closed", inherit = mfGeomPlanarShape, 
                                           public = list(
                                             #' @description Initialize planar shape geometry with warping 
                                             #' for data given in (long) FDboost format. 
                                             #' Same as in mfGeomPlanarShape but storing also function arguments. 
                                             #' 
                                             #' @param data data in the format used in \code{\link{FDboost}}.
                                             #' @param formula formula describing the internal structure of the data, 
                                             #' intended for the \code{obj.formula} of \code{mfboost} 
                                             #' (interpreted by \code{mfInterpret_objformula}). 
                                             #' If a \code{formula} is already stored in the \code{mfGeometry} object,
                                             #' this is taken as default.
                                             #' @param weight_fun a function computing inner product weights for the geometry
                                             #' taking two arguments: \code{arg}, the argument of the functional data 
                                             #' specified in \code{formula}, \code{range} the range of \code{arg}.
                                             #' @param arg_range the \code{range} supplied to the \code{weight_fun}.
                                             #' 
                                             initialize = function(data, formula, 
                                                                   weight_fun = trapez_weights, 
                                                                   arg_range = c(0,1), warp_update) {
                                               # dataless initializations
                                               if(!is.null(weight_fun)) {
                                                 stopifnot(is.function(weight_fun))
                                                 private$weight_fun <- weight_fun
                                               }
                                               
                                               if(!is.null(arg_range)) {
                                                 private$arg_range <- arg_range
                                               }
                                               
                                               if(!missing(warp_update)) {
                                                 stopifnot(is.function(warp_update))
                                                 private$warp_update <- warp_update
                                               }
                                               
                                               if(missing(formula)) {
                                                 formula <- private$.formula
                                               } else {
                                                 stopifnot(inherits(formula, "formula"))
                                                 private$.formula <- formula
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
                                               private$.curve_ <- private$close(
                                                 self$structure(data[[v$value]]))
                                               # store argument and add as attribute
                                               private$.arg_ <- data[[v$arg]][private$structure_index$real]
                                               attr(private$.curve_, "arg") <- 
                                                 private$.arg_ <- c(private$.arg_, 1-private$.arg_[1])
                                               
                                               # store SRV trafo of curve 
                                               private$.y_ <- self$srv_trafo(private$.curve_)
                                               
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
                                               self$register()
                                               stopifnot(length(private$.weights_) == length(private$.y_) |
                                                           is.null(private$.weights_))
                                               
                                               invisible(self)
                                             },
                                             #' @description convert 2D vectors (given in long format with dimension 
                                             #' indicator) to complex
                                             structure = function(y) {
                                               if(is.null(private$structure_index)) 
                                                 stop("Please, initialize geometry before using $structure.")
                                               structure( complex(
                                                 real = y[private$structure_index$real], 
                                                 imag = y[private$structure_index$imaginary]),
                                                 arg = private$.arg_ )
                                             },
                                             #' @description scale \code{y_} to unit length or
                                             #' an SRV-trnasform \code{q_} to unit norm, respectively.
                                             register = function(y_) {
                                               if(missing(y_)) {
                                                 L <- self$innerprod(private$.y_)
                                                 private$.curve_ <- private$.curve_ / L
                                                 return((private$.y_ <- private$.y_ / sqrt(L)))
                                               }
                                               y_ / sqrt(self$innerprod(y_))
                                             },
                                             align = function(y_, y0_, 
                                                              warp = self$warp, 
                                                              eps = 0.01) {
                                               y_ <- self$register(
                                                 if(missing(y_)) {
                                                   private$.warp(y0_ = y0_,
                                                                 optimize = warp,
                                                                 eps = eps) } else {
                                                                   private$.warp(y_, y0_,
                                                                                 optimize = warp,
                                                                                 eps = eps)
                                                                 })
                                               if(!is.null(private$warp_update)) 
                                                 self$warp <- private$warp_update(
                                                   warp_memory = self$warp_memory, warp = warp)
                                               self$warp_memory <- self$warp_memory + warp
                                               super$align(y_, y0_)
                                             },
                                             distance = function(y0_, y1_, squared = FALSE) {
                                               super$distance(y0_ = y0_, 
                                                              y1_ = if(missing(y1_)) 
                                                                self$align(y0_ = y0_) else 
                                                                  self$align(y_ = y1_, y0_ = y0_),
                                                              squared = squared)
                                             },
                                             #' @description Perform SRV-transform of a curve and its inverse.
                                             #' @param inverse logical, should the inverse SRV-transform be conducted.
                                             #' @param center logical, if \code{inverse = TRUE}, should the output be centered.
                                             srv_trafo = function(y_, inverse = FALSE, center = FALSE) {
                                               # if(inverse & center) 
                                               #   warning("Sorry, centering not implemented, yet.")
                                               arg <- attr(y_, "arg")
                                               long_arg <- length(arg) == length(y_) + 1 & 
                                                 arg[1] %% 1 == tail(arg, 1) %% 1
                                               if(inverse) {
                                                 ret <- elasdics::get_points_from_srv(data.frame(
                                                   t = if(long_arg) arg[-length(arg)] else 
                                                     arg,
                                                   x = Re(y_), 
                                                   y = Im(y_)))
                                               } else {
                                                 ret <- elasdics::get_srv_from_points(data.frame(
                                                   t = arg,
                                                   x = if(long_arg) 
                                                     private$close(Re(y_)) else Re(y_), 
                                                   y = if(long_arg) 
                                                     private$close(Im(y_)) else Im(y_)))
                                               }
                                               ret <- structure(complex(
                                                 re = ret$x, im = ret$y
                                               ), arg = attr(y_, "arg"))
                                               
                                               if(inverse & center) 
                                                 ret <- private$center(ret)
                                               
                                               ret
                                             },
                                             #' @description default plotting function for shapes of plane curves, plotting \code{y_}
                                             #' in front of \code{y0_} (after alignment).
                                             #' @param level one of 'curve' (default) and 'SRV' indicating
                                             #' whether curves are to be plotted or their SRV-transforms
                                             #' @param closed logical, should the plotted functions be closed?
                                             #' @param center logical, if \code{level='curve'}, should curves be centered?
                                             #' @param add_original logical, should original curve y_ be displayed besides its aligned version. 
                                             #' @param col,pch,type graphical parameters passed to \code{base::plot} referring to \code{y_}.
                                             #' @param ylab,xlab,xlim,ylim,xaxt,yaxt,asp graphical parameters passed to \code{base::plot} 
                                             #' with modified defaults.
                                             #' @param ... other arguments passed to \code{base::plot}.
                                             #' @param y0_par graphical parameters for \code{y0_}.
                                             #' @param seg_par graphical parameters for line segments connecting \code{y_} 
                                             #' and \code{y0_}.
                                             plot = function(y_, y0_ = self$pole_, 
                                                             level = c("curve", "SRV"), 
                                                             closed = TRUE,
                                                             center = TRUE,
                                                             add_original = FALSE,
                                                             ylab = NA, xlab = NA, 
                                                             col = "black",
                                                             pch = 19, asp = 1,
                                                             xlim, ylim,
                                                             yaxt='n',
                                                             xaxt='n',
                                                             type = "l",
                                                             y0_par = list(col = "darkgrey", type = type), 
                                                             seg_par = list(col = "grey"), ...) {
                                               level <- match.arg(level)
                                               stopifnot(is.logical(closed))
                                               
                                               missy_ <- missing(y_)
                                               
                                               if(!is.null(y0_)) {
                                                 # align y_
                                                 y_ <- if(missy_) self$align(y0_ = y0_) else 
                                                   self$align(y_, y0_)
                                                 
                                                 if(level == "curve") {
                                                   y_ <- self$srv_trafo(y_, inverse = TRUE, center = center)
                                                   y0_ <- self$srv_trafo(y0_, inverse = TRUE, center = center)
                                                 }
                                                 
                                                 if(is.null(y0_par)) 
                                                   y0_par <- list()
                                                 if(is.list(y0_par)) {
                                                   .y0_par <- c(private$.y0_par, type = type)
                                                   nm <- setdiff(names(.y0_par), names(y0_par))
                                                   y0_par[nm] <- .y0_par[nm]
                                                   y0_par$x <- Re(y0_)
                                                   y0_par$y <- Im(y0_)
                                                 } else {
                                                   stop("y0_par has to be supplied as a list (or NULL).")
                                                 }
                                                 
                                                 if(!is.na(seg_par)) {
                                                   if(is.null(seg_par)) 
                                                     seg_par <- list()
                                                   if(is.list(seg_par)) {
                                                     nm <- setdiff(names(private$.seg_par), names(seg_par))
                                                     seg_par[nm] <- private$.seg_par[nm]
                                                     seg_par[c("x0", "y0", "x1", "y1")] <- 
                                                       list(x0 = Re(y0_), y0 = Im(y0_), 
                                                            x1 = Re(y_), y1 = Im(y_))
                                                   } else {
                                                     stop("seg_par has to be supplied as a list (or NULL).")
                                                   }
                                                 } 
                                               } else {
                                                 if(level == "curve") {
                                                   y_ <- if(missy_) 
                                                     private$.curve_ else 
                                                       self$srv_trafo(y_, inverse = TRUE, center = center)
                                                 } else {
                                                   if(missy_) y_ <- private$.y_
                                                 }
                                               }
                                               
                                               if(closed) {
                                                 y_ <- private$close(y_)
                                                 y0_ <- private$close(y0_)
                                               } 
                                               
                                               if(missing(xlim))  
                                                 xlim <- range(Re(c(y_, y0_)))
                                               if(missing(ylim))
                                                 ylim <- range(Im(c(y_, y0_)))
                                               
                                               plot(x = Re(y_), y = Im(y_), ylab = ylab, xlab = xlab, 
                                                    yaxt = yaxt, xaxt = xaxt, type = type,
                                                    pch = pch, col = if(is.null(y0_)) col,
                                                    xlim = xlim, ylim = ylim, asp = asp, ...)
                                               
                                               if(!is.null(y0_)) {
                                                 if(add_original) {
                                                   if(level == "curve") {
                                                     y_original <- if(missy_) private$center(private$.curve_) else 
                                                       self$srv_trafo(y_original, inverse = TRUE, center = center)
                                                   } else {
                                                     y_original <- if(missy_) 
                                                       private$.y_ else y_
                                                   }
                                                   
                                                   if(closed)
                                                     y_original <- private$close(y_original)
                                                   points(x = Re(y_original), y = Im(y_original), pch = pch, col = col, type = type, ...)
                                                 }
                                                 do.call(segments, seg_par)
                                                 do.call(points, y0_par)
                                                 points(x = Re(y_), y = Im(y_), pch = pch, col = col, type = type, ...)
                                               }
                                             },
                                             warp = TRUE,
                                             warp_memory = 0
                                           ), active = list(
                                             y_ = function(value) {
                                               if(missing(value))
                                                 private$.y_ else {
                                                   private$.y_ <- self$validate(value)
                                                   private$.curve_ <- self$srv_trafo(value, inverse = TRUE)
                                                 } }
                                           ),
                                           private = list(
                                             .curve_ = NULL,
                                             .arg_ = NULL,
                                             close = function(x) x[c(seq_len(length(x)), 1)],
                                             warp_update = NULL,
                                             center = function(curve_) {
                                               ONE <- rep(1, length(curve_) -1)
                                               curve_ <- curve_ - self$innerprod(ONE, head(curve_, -1)) / self$innerprod(ONE)
                                             },
                                             .warp = function(y_, y0_, closed = TRUE, optimize = TRUE, eps = .01) {
                                               find_t_args <- list()
                                               # align to
                                               find_t_args$p <- rbind(
                                                 x = Re(y0_), 
                                                 y = Im(y0_))
                                               # at time points
                                               find_t_args$r = if(is.null(attr(y0_, "arg"))) {
                                                 if(closed) c(private$.arg_, 1+private$.arg_[1])  else 
                                                   private$.arg_
                                               } else
                                                   attr(y0_, "arg")
                                               
                                               if(missing(y_)) {
                                                 # align this
                                                 find_t_args$q <- rbind(
                                                   x = Re(private$.y_),
                                                   y = Im(private$.y_)
                                                 )
                                                 # at time points
                                                 find_t_args$s <- 
                                                   if(is.null(attr(private$.y_, "arg"))) 
                                                     private$.arg_ else
                                                       attr(private$.y_, "arg")
                                                 find_t_args$initial_t <- if(is.null(attr(private$.y_, "arg_optim"))) 
                                                   find_t_args$s else attr(private$.y_, "arg_optim")
                                               } else {
                                                 find_t_args$q <- rbind(
                                                   x = Re(y_), 
                                                   y = Im(y_) )
                                                 # at time points
                                                 find_t_args$s <- if(is.null(attr(y_, "arg"))) 
                                                   private$.arg_ else
                                                     attr(y_, "arg")
                                                 find_t_args$initial_t <- if(is.null(attr(y_, "arg_optim"))) 
                                                   find_t_args$s else attr(y_, "arg_optim")
                                               }
                                               
                                               if(optimize) {
                                                 find_t_args$eps <- eps
                                                 t_optim <- suppressMessages(do.call(
                                                   if(closed) 
                                                     elasdics:::find_optimal_t_discrete_closed else
                                                       elasdics:::find_optimal_t_discrete,
                                                   find_t_args
                                                 ))
                                                 
                                                 if(missing(y_)) {
                                                   attr(private$.y_, "arg_optim") <- t_optim
                                                 }
                                               } else {
                                                 t_optim <- find_t_args$s
                                               }
                                               
                                               if(missing(y_))
                                                 curve_ <- private$.curve_ else
                                                   curve_ <- self$srv_trafo(y_, inverse = TRUE, center = FALSE)
                                               
                                               curve_ <- structure(
                                                 complex(
                                                   re = private$approx(t_optim, Re(curve_), xout = find_t_args$r)$y,
                                                   im = private$approx(t_optim, Im(curve_), xout = find_t_args$r)$y
                                                 ),
                                                 arg = find_t_args$r
                                               )
                                               self$srv_trafo(curve_)
                                             },
                                             approx = function(x, y = NULL, xout, closed = TRUE, xleft = 0, xright = 1, ...) {
                                               if(closed) {
                                                 x <- (x - xleft) %% (xright - xleft) + xleft
                                                 xout <- (xout - xleft) %% (xright - xleft) + xleft
                                                 xorder <- order(x)
                                                 x <- x[xorder]
                                                 y <- y[xorder]
                                                 xrange <- c(x[1], x[length(x)])
                                                 where <- range(findInterval(xout, xrange))
                                                 if(any(where != 1)) {
                                                   y <- c(if(0 %in% where) tail(y, 1),
                                                          y,
                                                          if(2 %in% where) head(y, 1))
                                                   x <- c(if(0 %in% where) xleft + tail(x, 1) -xright,
                                                          x,
                                                          if(2 %in% where) xright + head(x, 1)-xleft)
                                                 }
                                               }
                                               approx(x = x, y = y, xout = xout, ...)
                                             }
                                           ))

# elastic planar shape regression family -------------------------------

#' @param weight_fun a function producing inner product weights 
#' taking the arguments \code{arg} (vector of arguments of the function) and 
#' \code{range} (range of the arguments). Passed to \code{mf$initialize}.
#' @param ... further arguments passed to \code{mfboost} in the pole fit.
#' 
#' @export
#' @name ElasticPlanarShapeL2
#' @rdname mfFamily
ElasticPlanarShapeL2 <- function(weight_fun = trapez_weights, closed = TRUE, formula, warp_update = function(warp_memory, warp) TRUE, ...) {
  if(!closed) stop("Sorry, only closed case implemented so far.")
  mf <- mfGeomProduct$new(
    mfGeom_default = mfGeomElasticPlanarShape_closed$new(weight_fun = weight_fun, 
                                                         arg_range = c(0,1), warp_update = warp_update))
  
  pole <- mfPoleRiemannL2$new(
    mfGeom = mfGeomProduct$new(
      mfGeom_default = mfGeomSRV_closed$new(weight_fun = weight_fun, 
                                            arg_range = c(0,1))),
    mfPole = mfPoleZero$new(mfGeom = mf), ...)
  
  if(!missing(formula)) {
    mf$formula <- pole$mf$formula <- formula
  }
  
  RiemannL2(mf = mf, pole = pole)
}
