# Planar shape space (single vector) ------------------------------------------

#' Planar Shape Space Geometry with simple warping alignment
#' 
#' @description Geometry of the 'classic' shape space identifying planar shapes
#' with a centered and scaled complex vector on the complex sphere / complex
#' projective space, yet also with warping alignment.
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
#' @export
mfGeomWarpPlanarShapeSimple <- R6Class("mfGeomWarpPlanarShapeSimple", inherit = mfGeomPlanarShape, 
                               private = list(
                               .distance = function(y0_, y1_, squared = FALSE) {
                                 if(missing(y1_)) {
                                   y1_ <- self$align(private$.y_, y0_) } else 
                                     y1_ <- self$align(y1_, y0_)
                                   
                                 ip <- private$.innerprod(y0_, y1_)
                                 if(Mod(ip-1) < 1e-15) return(0)
                                 if(Mod(ip+1) < 1e-15) return(pi)
                                 ret <- acos( Mod(ip) ) 
                                 if(squared)
                                   return(ret^2)
                                 ret
                               },
                               .align = function(y_, y0_) {
                                 if(missing(y_)) {
                                   y_ <- private$.warp(private$.y_, y0_) } else 
                                     y_ <- private$.warp(y_, y0_)
                                 y_ <- self$register(y_)
                                 rot <- private$.innerprod(y_, y0_)
                                 return( rot * y_ / Mod(rot) )
                               },
                               .warp = function(y_, y0_, closed = FALSE, eps = .01) {
                                 w <- elasdics::align_curves(
                                     data.frame(
                                       x = Re(y0_),
                                       y = Im(y0_)
                                     ), data.frame(
                                       x = Re(y_),
                                       y = Im(y_)
                                     ),
                                     eps = eps,
                                     closed = closed
                                   )
                                 
                                 thist <- which(
                                   names(w$data_curve2_aligned)%in%
                                     c("t", "t_optim"))
                                 structure(
                                   complex(
                                     re = approx(w$data_curve2_aligned$t_optim, 
                                                 w$data_curve2_aligned[-thist][[1]], 
                                                 xout = w$data_curve1$t)$y,
                                     im = approx(w$data_curve2_aligned$t_optim, 
                                                 w$data_curve2_aligned[-thist][[2]], 
                                                 xout = w$data_curve1$t)$y
                                   ),
                                   arg = w$data_curve1$t
                                 )}
                             ))

# warping aligned planar shape regression family (simple) ----------------------

#' @param weight_fun a function producing inner product weights 
#' taking the arguments \code{arg} (vector of arguments of the function) and 
#' \code{range} (range of the arguments). Passed to \code{mf$initialize}.
#' @param arg_range vector of length 2 specifying the \code{range} argument of 
#' the \code{weight_fun}. The default \code{NULL} will take the minimum and maximum 
#' of the supplied \code{arg}.
#' 
#' @export
#' @name WarpPlanarShapeL2Simple
#' @rdname mfFamily
WarpPlanarShapeL2Simple <- function(pole.type = "RiemannL2", pole.control = boost_control(), 
                              weight_fun = equal_weights, arg_range = NULL) {
  mf <- mfGeomProduct$new(
    mfGeom_default = mfGeomWarpPlanarShapeSimple$new(weight_fun = weight_fun, arg_range = arg_range))
  
  RiemannL2(mf = mf, pole.type = pole.type, pole.control = pole.control)
}


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
                                                         arg_range = NULL,
                                                         closed, warp_update) {
                                     # dataless initializations
                                     if(!is.null(weight_fun)) {
                                       stopifnot(is.function(weight_fun))
                                       private$weight_fun <- weight_fun
                                     }
                                     
                                     if(!is.null(arg_range)) {
                                       private$arg_range <- arg_range
                                     }
                                     
                                     if(!missing(closed))
                                       private$closed <- closed
                                     
                                     if(!missing(warp_update)) {
                                       stopifnot(is.function(warp_update))
                                       private$warp_update <- warp_update
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
                                   align = function(y_, y0_, 
                                                    closed = private$closed, 
                                                    warp = self$warp, 
                                                    eps = 0.01) {
                                     y_ <- self$register(
                                         if(missing(y_)) {
                                           private$.warp(y0_ = y0_, 
                                                       closed = closed, 
                                                       optimize = warp,
                                                       eps = eps) } else {
                                                           private$.warp(y_, y0_, 
                                                                         closed = closed,
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
                                   warp = TRUE,
                                   warp_memory = 0
                                 ),
                             private = list(
                               .y_dat = NULL,
                               .arg_ = NULL,
                               closed = FALSE,
                               warp_update = NULL,
                               .warp = function(y_, y0_, closed = FALSE, optimize = TRUE, eps = .01) {
                                 d0_ <- data.frame(
                                   x = Re(y0_), 
                                   y = Im(y0_),
                                   t = if(is.null(attr(y0_, "arg"))) 
                                     private$.arg_ else
                                     attr(y0_, "arg") )
                                 
                                 # set up output structure
                                 d_ <- if(missing(y_)) private$.y_dat else 
                                   data.frame(
                                     x = Re(y_), 
                                     y = Im(y_),
                                     t = if(is.null(attr(y_, "arg"))) 
                                       private$.arg_ else
                                         attr(y_, "arg") )
                                 
                                 if(optimize) {
                                   d_$t[1] <- 0
                                   d_$t <- elasdics::align_curves(
                                     d0_, d_,
                                     eps = eps,
                                     closed = closed
                                   )$data_curve2_aligned$t_optim[seq_along(d_$t)]
                                   
                                   if(closed) {
                                     # sort if necessary
                                     new_order <- if(any(diff(d_$t) < 0)) 
                                       order(d_$t) else NA
                                     if(!is.na(new_order[1])) 
                                       d_ <- d_[new_order, ]
                                   }
                                   
                                   if(missing(y_)) {
                                     # update internal parameterization
                                     if(is.na(new_order)[1]) {
                                       private$.y_dat$t <- d_$t
                                     } else {
                                       private$.y_dat <- d_
                                       private$.y_ <- private$.y_[new_order]
                                     }
                                     attr(private$.y_, "arg") <- d_$t
                                   } 
                                 }
                                 
                                 structure(
                                   complex(
                                     re = private$approx(d_$t, d_$x, xout = d0_$t, 
                                                         closed = closed)$y,
                                     im = private$approx(d_$t, d_$y, xout = d0_$t, 
                                                         closed = closed)$y
                                   ),
                                   arg = d0_$t
                                 )},
                               approx = function(x, y = NULL, xout, closed = FALSE, xleft = 0, xright = 1, ...) {
                                 if(closed) {
                                   rx <- range(x)
                                   where <- range(findInterval(xout, rx))
                                   if(any(where != 1)) {
                                     if(!is.null(y))
                                       y <- c(if(0 %in% where) y[which.max(x)], 
                                              y, 
                                              if(2 %in% where) y[which.min(x)])
                                     x <- c(if(0 %in% where) xleft + rx[2]-xright, 
                                            x, 
                                            if(2 %in% where) xright + rx[1]-xleft)
                                   }
                                 }
                                 approx(x = x, y = y, xout = xout, ...)
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
#' @name WarpPlanarShapeL2
#' @rdname mfFamily
WarpPlanarShapeL2 <- function(pole.type = "RiemannL2", pole.control = boost_control(), 
                          weight_fun = equal_weights, closed = FALSE, warp_update = function(warp_memory, warp) TRUE) {
  mf <- mfGeomProduct$new(
    mfGeom_default = mfGeomWarpPlanarShape$new(weight_fun = weight_fun, 
                                               arg_range = c(0,1), closed = closed, warp_update = warp_update))
  
  RiemannL2(mf = mf, pole.type = pole.type, pole.control = pole.control)
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
mfGeomElasticClosedPlanarShape <- R6Class("mfGeomElasticClosedPlanarShape", inherit = mfGeomPlanarShape, 
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
                                       y_original <- if(missy_) 
                                         private$.y_ else y_
                                       y_ <- self$align(y_original, y0_)
                                       
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
                                         if(level == "curve") 
                                           y_original <- if(missy_) private$center(private$.curve_) else 
                                             self$srv_trafo(y_original, inverse = TRUE, center = center)
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
                                     find_t_args$r = if(is.null(attr(y0_, "arg"))) 
                                         private$.arg_ else
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
                                       t_optim <- do.call(
                                         if(closed) 
                                           elasdics:::find_optimal_t_discrete_closed else
                                             elasdics:::find_optimal_t_discrete,
                                           find_t_args
                                         )
                                       
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
#' @param arg_range vector of length 2 specifying the \code{range} argument of 
#' the \code{weight_fun}. The default \code{NULL} will take the minimum and maximum 
#' of the supplied \code{arg}.
#' 
#' @export
#' @name WarpPlanarShapeL2
#' @rdname mfFamily
ElasticPlanarShapeL2 <- function(pole.type = "RiemannL2", pole.control = boost_control(), 
                              weight_fun = equal_weights, closed = FALSE, warp_update = function(warp_memory, warp) TRUE) {
  mf <- mfGeomProduct$new(
    mfGeom_default = mfGeomElasticPlanarShape$new(weight_fun = weight_fun, 
                                               arg_range = c(0,1), closed = closed, warp_update = warp_update))
  
  RiemannL2(mf = mf, pole.type = pole.type, pole.control = pole.control)
}
