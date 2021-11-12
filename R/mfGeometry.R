setOldClass("mfGeometry")

#' R6 Class implementing manifold geometries
#'
#' An abstract \code{mfGeometry} class implementing a set of methods, including 
#' e.g. Exp- and Log-mappings, that are needed to define \code{mfboost} families
#' in a concise and modular way. 
#' The arguments of the methods are vectors, matrices
#' or arrays possibly holding additional attributes. Besides these methods, objects
#' of class \code{mfGeometry} offer an active slot for the 'pole' which serves 
#' as prototype for elements of the manifold. Children of this class 
#' specify the geometry of interest, such as \code{\link{mfGeomEuclidean}}, 
#' \code{\link{mfGeomSphere}}, \code{\link{mfGeomPlanarShape}}, and
#' \code{\link{mfGeomPlanarSizeShape}}. 
#'
#' @param y a numeric vector containing the elements of an object y in long format
#' @param v a numeric vector containing the elements of a tangent vector in long format
#' @param y_ an object in the internal format of the geometry 
#' @param y0_ an object in the internal format of the geometry
#' @param y1_ an object in the internal format of the geometry
#' @param v_ an object repesenting a tangent vector (internal format)
#' @param v0_ an object repesenting a tangent vector (internal format)
#' @param v1_ an object repesenting a tangent vector (internal format)
#' @param weights_ inner product weights in the internal format
#' @param ... other arguments.
#'
#' @field y_ slot for a data object in the internal format
#' @field pole_ slot for a data object in the internal format
#' @field weights_ slot for inner product weights in the internal format 
#'
#'
# #' @name mfGeometry
# #' @rdname mfGeometry
#' @import R6
#' @export

mfGeometry <- R6Class("mfGeometry", 
  public = list(
    
    ## data communicators with long format
    
    #' @description structure a vector containing the elements of the response 
    #' \code{y} in long format in the internal format of the geometry.
   structure = function(y) {
     fun <- as.character(match.call()[[1]])[3]
     private$.undefined(fun) },
    #' @description inverse operation of \code{structure}.
   unstructure = function(y_) {
     fun <- as.character(match.call()[[1]])[3]
     private$.undefined(fun) },
    #' @description convert internal inner product weights to a weight vector
    #' matching the un-structered response \code{y} in long format.
   unstructure_weights = function(weights_) NULL,
   
   ## intrinsic functions
    
    #' @description a function aligning \code{y_} to \code{y0_} primarly intended
    #' for quotient geometries. 
   align = function(y_, y0_ = private$.pole_) {
     fun <- as.character(match.call()[[1]])[3]
     private$.undefined(fun) },
    #' @description a function registering \code{y_} to a proper subset, say to
    #' a subspace or an embedded submanifold.
   register = function(y_) {
     fun <- as.character(match.call()[[1]])[3]
     private$.undefined(fun) },
   
    #' @description computes the distance or a vector of distances between objects.
   distance = function(y0_, y1_)
     private$.distance(y0_, y1_, ...),
    #' @description the Riemannian Exp map mapping a tangent vector \code{v_} 
    #' at \code{y_} to the manifold.
   exp = function(v_, y0_ = private$.pole_, ...) {
     fun <- as.character(match.call()[[1]])[3]
     private$.undefined(fun) },
    #' @description the inverse of the Riemannian Exp map.
   log = function(y_, y0_ = private$.pole_, ...) {
     fun <- as.character(match.call()[[1]])[3]
     private$.undefined(fun) },
    #' @description parallel transport of a tangent vector \code{v0_} 
    #' at \code{y0_} to the tangent space at \code{y1_}.
   transport = function(v0_, y0_, y1_, ...) {
     fun <- as.character(match.call()[[1]])[3]
     private$.undefined(fun) },
    #' @description computes the inner product or a vector of inner products 
    #' between objects.
   innerprod = function(v0_, v1_ = v0_, ...)
     private$.innerprod(v0_, v1_, ...),
    #' @description returns the normal vector (in unstructured format) to the 
    #' tangent space at a point \code{y0_} considered as a linear subspace.
   get_normal = function(y0_ = private$.pole_, ...) {
     fun <- as.character(match.call()[[1]])[3]
     private$.undefined(fun) },
   #' @description function for validating structure of internal object \code{y_}.
   validate = function(y_, ...) {
     stopifnot(TRUE)
     y_
   }
   
   ), 
  private = list(
    .pole_ = NULL,
    .y_ = NULL,
    .weights_ = NULL,
    .undefined = function(fun) {
      stop(paste0("The method '", fun, "' is not defined in this geometry.
          To define it, create a new one as a subclass of the current mfGeometry (subclass)."))},
    # methods used by other methods:
    .distance = function(y0_, y1_, ...) {         
      fun <- as.character(match.call()[[1]])[3]
      private$.undefined(fun) }, 
    .innerprod = function(v0_, v1_ = v0_, ...) {
      fun <- as.character(match.call()[[1]])[3]
      private$.undefined(fun) }
    ),
  
  active = list(
    pole_ = function(value) {
    if(missing(value))
      private$.pole_ else {
        private$.pole_ <- self$validate(value)
      } },
    y_ = function(value) {
      if(missing(value))
        private$.y_ else {
          private$.y_ <- self$validate(value)
        } }, 
    weights_ = function(value) {
          if(missing(value))
            private$.weights_ else {
              private$.weights_ <- value
            } })    
  ) #mfGeometry


# Standard euclidean geometry for a single vector  ------------------------

#' Euclidean Geometry Class
#' 
#' An R6 class defining the Euclidean geometry for a single vector. 
#' This class is typically not too much of practical interest, as in this case
#' the \code{Family} could typically by defined without it. However, it serves
#' as parent for most other geometries.
#' 
#' @param y a numeric vector
#' @param v a numeric vector
#' @param y_ a numeric/complex vector
#' @param y0_ a numeric/complex vector
#' @param y1_ a numeric/complex vector
#' @param v_ a (tangent) numeric/complex vector
#' @param v0_ a (tangent) numeric/complex vector
#' @param v1_ a (tangent) numeric/complex vector
#' @param weights a numeric vector
#' @param weights_ a numeric vector of the same length as \code{y_} / \code{pole_}
#' 
#' @export
mfGeomEuclidean <- R6Class("mfGeomEuclidean", inherit = mfGeometry,
  public = list(
    #' @description method for structuring either objects \code{y}, 
    #' tangent vectors \code{v}, or weight vectors \code{weights}. Usually the 
    #' identity, except for complex vectors. FOR COMPLEX VECTORS USE WITH CARE:
    #' STILL ORDER DEPENDEND.
    structure = function(y, v, weights, y0_ = private$.y_) {
      if(missing(y) + missing(v) + missing(weights) != 2)
        stop("One and only one of 'y', 'v' and 'weights' have to be provided.")
      if(is.null(y0_)) stop("Please specify a reference y0_ or a pole beforehand.")
      if(!missing(weights)) {
        stopifnot(is.numeric(weights))
        if(is.complex(y0_)) {
          weights <- weights[1L:length(y0_)]
        }
        return(weights)
      } 
     if(missing(y)) y <- v
     stopifnot(is.numeric(y))
     if(is.complex(y0_)) {
       l <- length(y)/2
       y_ <- complex(real = y[seq_len(l)], imaginary = y[(l+1L):(2L*l)])
     }
     if(!is.null(y0_)) stopifnot(length(y_) == length(y0_))
     y_
   },
   #' @description inverse of \code{$structure}. FOR COMPLEX VECTORS USE WITH CARE:
   #' STILL ORDER DEPENDEND.
   unstructure = function(y_, v_, weights_, y0_ = private$.y_) {
     if(missing(y_) + missing(v_) + missing(weights_) != 2)
       stop("One and only one of 'y_', 'v_' and 'weights_' has to be provided.") 
     if(is.null(y0_)) stop("Please specify a reference y0_ or a pole beforehand.")
     if(missing(weights_)) {
       if(missing(y_)) y_ <- v_ 
       if(is.complex(y0_)) c(Re(y_), Im(y_)) else
         y_
     } else {
         if(is.complex(y0_)) {
           if(is.null(dim(weights_))) 
             c(weights_, weights_) else 
               complex2realmat(weights_)
         } else weights_
       }
   },
   #' @description The identity.
   align = function(y_, y0_) private$.align(y_, y0_),
   #' @description The identity.
   register = function(y_) y_,
   #' @description Simple addition.
   exp = function(v_, y0_ = private$.pole_) y0_ + v_,
   #' @description Simple substraction.
   log = function(y_, y0_ = private$.pole_) y_ - y0_,
   #' @description The identity.
   transport = function(v0_, y0_, y1_) v0_,
   #' @description The weighted scalar product.
   innerprod = function(v0_, v1_ = v0_, weights_ = private$.weights_)
     private$.innerprod(v0_, v1_, weights_),
   #' @description Always returning \code{NULL}.
   get_normal = function(y0_ = private$.pole_) NULL,
   #' @description Check whether numeric or complex.
   validate = function(y_) {
     stopifnot(is.numeric(y_) | is.complex(y_))
     y0_
   }
), private = list(
  .innerprod = function(v0_, v1_ = v0_, weights_ = private$.weights_) {
    # allow for potentially weighted inner product
    if(is.null(weights_)) return(as.complex(sum(Conj(v0_)*v1_)))
    if(is.null(dim(weights_)))  return(as.complex(sum(weights_*Conj(v0_)*v1_)))
    return( as.complex(crossprod(Conj(v0_), weights_)%*%v1_) )
  },
  .distance = function(y0_, y1_, squared = FALSE) {
    d2 <- private$.innerprod(y1_-y0_)
    if(squared) d2 else sqrt(d2)
  },
  .align = function(y_, y0_) y_
))

# Euclidean sphere for a single vector --------------------------------------

#' Geometry of the Unit Sphere
#' 
#' @description The geometry of the n-dimensional unit sphere, treating the complex sphere is, 
#' in case, as a real manifold.
#' 
#' @param y a numeric vector representing an element of the manifold in an 
#' unstructured way
#' @param v a numeric vector representing an tangent vector in an 
#' unstructured way
#' @param y_ a numeric/complex vector on the sphere
#' @param y0_ a numeric/complex vector on the sphere
#' @param y1_ a numeric/complex vector on the sphere
#' @param v_ a numeric/complex tangent vector
#' @param v0_ a numeric/complex tangent vector
#'
#' @export
mfGeomUnitSphere <- R6Class("mfGeomUnitSphere", inherit = mfGeomEuclidean,
 public = list(
  #' @description \code{y} vectors are scaled to unit norm 
  register = function(y_) {
    y_ / sqrt(private$.innerprod(y_)) },
  #' @description vectors \code{v} 
  #' are orthogonalized with respect to \code{y0_} guaranteeing proper tangent 
  #' vectors.
  register_v = function(v_, y0_ = private$.pole_) {
    v_ - y0_ * Re( private$.innerprod(y0_, v_) ) # complex sphere as real manifold!
  },
  #' @param method character string, for choosing one of two slightly different 
  #' methods to implement log-method ("simple" or "alternative").
  log = function(y_, y0_ = private$.pole_, method = c("simple", "alternative")) {
    method = match.arg(method)
    switch(method, 
           simple = {
             ### See for example Dryden&Mardia2012, page 76 -> described for complex vectors, 
             ### slightly adjusted with x0 * innerprod(x0, x) instead of x0 * Mod(innerprod(x0, x))
             rho <- private$.distance(y_, y0_)
             if(!is.na(rho)) if(rho < 1e-15) return(rep(0, length(y))) 
             return( rho / sin(rho) * 
                       ( y_ - y0_ * private$.innerprod(y0_, y_) )) 
           },
           alternative = {
             ### See for example Cornea et al. 2017
             v_ <- y_ - y0_ * private$.innerprod(y0_, y_)
             if(sqrt(private$.innerprod(v_)) == 0) return(v_)
             else return( v_ * private$.distance(x_, x0_) / 
                            sqrt(private$.innerprod(v_)))  
           }) },
  #' @description moving \code{v_} along an arc on the sphere.
  exp = function(v_, y0_) {
    if(private$.innerprod(v_) == 0) return(y0_)
    s <- sqrt(private$.innerprod(v_))
    return( y0_ * cos(s) + v_ / s * sin(s) )   
  },
  #' @param method expression used for parallel transport: "general" corresponds to the one of
  #' Cornea et al. 2017, "horizontal" to the one of Dryden & Mardia 2012 p. 76, 
  #' and "simple" to a simplified version of "general" where \code{v0_} 
  #' is horizontal to \code{y0_}.
  transport = function(v0_, y0_, y1_, method = c("general", "horizontal", "simple")) {
    method <- match.arg(method)
    # if( Im(private$.innerprod(y0_, y1_)) > tol & mode(v0_) == "complex") stop("For complex vector methods only work in the horizontal case, yet.")
    if( method == "horizontal" & mode(v0_) == "numeric") stop("horizontal method only works for the horizontal complex case")
    switch(method,
           simple = {
             # assuming v0_ tangent to y0_
             y1tilde_ <- (y1_ - private$.innerprod(y0_, y1_)*y0_) / sqrt(1- private$.innerprod(y0_, y1_)^2)
             return( v0_ - private$.innerprod(y1tilde_, v0_) * y1tilde_ - sqrt(1 - private$.innerprod(y0_, y1_)^2) * private$.innerprod(y1tilde_, v0_) * y0_ + 
                       Conj(private$.innerprod(y0_, y1_))*private$.innerprod(y1tilde_, v0_) * y1tilde_ )
           },
           general = { ### See for example Cornea et al. 2017
             y1tilde_ <- (y1_ - private$.innerprod(y0_, y1_)*y0_) / sqrt(1- private$.innerprod(y0_, y1_)^2)
             return( v0_ - private$.innerprod(y0_, v0_) * y0_ - private$.innerprod(y1tilde_, v0_) * y1tilde_ + 
                       (private$.innerprod(y0_, y1_)*private$.innerprod(y0_, v0_) - 
                          sqrt(1 - private$.innerprod(y0_, y1_)^2) * private$.innerprod(y1tilde_, v0_)) * y0_ + 
                       (sqrt(1- private$.innerprod(y0_,y1_)^2) * private$.innerprod(y0_,v0_) + 
                          Conj(private$.innerprod(y0_, y1_))*private$.innerprod(y1tilde_, v0_)) * y1tilde_ )
           }, 
           horizontal = { 
             ### See for example Dryden&Mardia 2012 page 76, described for complex vectors
             return( v0_ - ( private$.innerprod(y1_, v0_) * (y0_ + y1_) ) / 
                       ( 1 + Mod(private$.innerprod(y0_, y1_)) ) )
           })
  },
  #' @description normal vector just corresponds to the pole itself
  get_normal = function(y0_ = private$.pole_) cbind(Re(y0_), Im(y0_))
 ), private = list(
   .distance = function(y0_, y1_, squared = FALSE) {
     ip <- private$.innerprod(y0_, y1_)
     if(Mod(ip-1) < 1e-15) return(0)
     if(Mod(ip+1) < 1e-15) return(pi)
     ret <- acos( Re(ip) )
     if(squared)
       return(ret^2)
     ret
   }
 ) )


# Planar shape space (single vector) ------------------------------------------

#' Planar Shape Space Geometry
#' 
#' @description Geometry of the 'classic' shape space identifying planar shapes
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
#' @export
mfGeomPlanarShape <- R6Class("mfGeomPlanarShape", inherit = mfGeomUnitSphere, 
 public = list(
   #' @description Initialize planar shape geometry for data given 
   #' in (long) FDboost format. 
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
     
     # with no data stop here !!!!!!!!!!!!!!!!!
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
   #' @description remove double weights to obtain 
   #' inner product weights for complex vectors 
   structure_weights = function(weights) {
     if(is.null(private$structure_index)) 
       stop("Please, initialize geometry before using $structure.")
     weights[private$structure_index$real]
     },
   #' @description center \code{y_} and scale it to unit norm.
   register = function(y_) {
     ONE <- rep(1, length(y_))
     y_ <- y_ - private$.innerprod(ONE, y_) / private$.innerprod(ONE)
     y_ / sqrt(private$.innerprod(y_)) 
   },
   #' @description orthogonally project \code{v_} into the tangent space of \code{y0_}.
   register_v = function(v_, y0_ = self$pole_) {
     v_ - y0_ * private$.innerprod(y0_, v_) 
   },
   #' @description Apply Log function of the sphere after rotation alignment
   #' @param method alternatives "simple" and "alternative" for the expression 
   #' used to compute the sphere Log-map. Passed to parent method.
   log = function(y_, y0_ = self$pole_, method = c("simple", "alternative")) {
    super$log(
      private$.align(y_, y0_)
      , y0_, method)
  },
  #' @description for the parallel transport previous alignment is assumed,
  #' such that the transport of the sphere is directly applied.
  #' @param method expression used for parallel transport: "general" corresponds to the one of
  #' Cornea et al. 2017, "horizontal" to the one of Dryden & Mardia 2012 p. 76, 
  #' and "simple" to a simplified version of "general" where \code{v0_} 
  #' is horizontal to \code{y0_}. Passed to parent method.
  transport = function(v0_, y0_, y1_, method = c("horizontal", "simple", "general")) {
    method <- match.arg(method)
    super$transport(v0_, y0_, y1_, method = method)
  },
  #' @description default plotting function for planar shapes, plotting \code{y_}
  #' in front of \code{y0_} (after alignment).
  #' @param col,pch,type graphical parameters passed to \code{base::plot} referring to \code{y_}.
  #' @param ylab,xlab,xlim,ylim,xaxt,yaxt,asp graphical parameters passed to \code{base::plot} 
  #' with modified defaults.
  #' @param ... other arguments passed to \code{base::plot}.
  #' @param y0_par graphical parameters for \code{y0_}.
  #' @param seg_par graphical parameters for line segments connecting \code{y_} 
  #' and \code{y0_}.
  plot = function(y_ = self$y_, y0_ = self$pole_, 
                  ylab = NA, xlab = NA, 
                  col = "black",
                  ylim = range(Im(c(y_, y0_))), 
                  xlim = range(Re(c(y_, y0_))),
                  pch = 19, asp = 1,
                  yaxt='n',
                  xaxt='n',
                  type = "p",
                  y0_par = list(col = "darkgrey", type = type), 
                  seg_par = list(col = "grey"), ...) {
    if(!is.null(y0_)) {
      
      # align y_
      y_ <- self$align(y_, y0_)
      
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
    }
    plot(x = Re(y_), y = Im(y_), ylab = ylab, xlab = xlab, 
         yaxt = yaxt, xaxt = xaxt, type = type,
         pch = pch, col = if(is.null(y0_)) col,
         xlim = xlim, ylim = ylim, asp = asp, ...)
    if(!is.null(y0_)) {
      do.call(segments, seg_par)
      do.call(points, y0_par)
      points(x = Re(y_), y = Im(y_), pch = pch, col = col, type = type, ...)
    }
  },
  #' @description check whether \code{is.complex(y_)}. 
  validate = function(y_) {
    stopifnot(is.complex(y_))
    y0_
  },
  #' @description Obtain "design matrix" of tangent space normal vectors in 
  #' unstructured long format.
  #' @param weighted logical, should inner product weights be pre-multiplied to 
  #' normal vectors?
  get_normal = function(y0_ = self$pole_, weighted = FALSE) {
    if(is.null(private$unstructure_index)) 
      stop("Please, initialize geometry before using $get_normal.")
    one <- rep(1, length(y0_))
    if(weighted & !is.null(private$.weights_)) {
      y0_ <- private$.weights_ * y0_
      one <- private$.weights_ * one
    }
    if(!is.null(y0_)) 
      cbind(complex2realmat(y0_), 
            complex2realmat(one))[private$unstructure_index, ]
  }
), private = list(
  .y0_par = list(col = "darkgrey"),
  .seg_par = list(col = "grey"),
  .distance = function(y0_, y1_, squared = FALSE) {
    ip <- private$.innerprod(y0_, y1_)
    if(Mod(ip-1) < 1e-15) return(0)
    if(Mod(ip+1) < 1e-15) return(pi)
    ret <- acos( Mod(ip) ) 
    if(squared)
      return(ret^2)
    ret
  },
  .align = function(y_, y0_) {
    rot <- private$.innerprod(y_, y0_)
    return( rot * y_ / Mod(rot) )
  },
  .unstructure = function(y_) c(Re(y_), Im(y_)),
  # needed for structuring and unstructuring
  structure_index = NULL,
  unstructure_index = NULL,
  weight_fun = NULL,
  arg_range = NULL
))


# Planar size and shape space ------------------------------------------

#' Planar Size and Shape Space Geometry
#' 
#' @description Geometry of the 'classic' size-and-shape space identifying planar 
#' configurations with a centered and complex vector and treating them as rotation
#' invariant.
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
mfGeomPlanarSizeShape <- R6Class("mfGeomPlanarSizeShape", inherit = mfGeomPlanarShape, 
                             public = list(
                               #' @description center \code{y_}.
                               register = function(y_) {
                                 ONE <- rep(1, length(y_))
                                 y_ <- y_ - private$.innerprod(ONE, y_) / private$.innerprod(ONE)
                               },
                               #' @description orthogonally project \code{v_} into the tangent space of \code{y0_}.
                               register_v = function(v_, y0_ = self$pole_) {
                                 v_ - y0_ * complex(im=1) * Im( private$.innerprod(y0_, v_) / private$.innerprod(y0_) )
                               },
                               #' @description \code{y_-y0_} after aligning \code{y_} to \code{y0_}.
                               log = function(y_, y0_ = self$pole_) {
                                 private$.align(y_, y0_) - y0_
                               },
                               #' @description \code{y0_+v_} assuming \code{v_} in the tangent space of \code{y0_}.
                               exp = function(v_, y0_ = private$.pole_) y0_ + v_,
                               #' @description size-and-shape parallel transport
                               #' @param method currently, only "horizontal" assuming \code{y0_, y1_} aligned
                               #' and \code{v0_} a proper horizontal tangent vector at \code{y0_}.
                               transport = function(v0_, y0_, y1_, method = c("horizontal")) {
                                 method <- match.arg(method)
                                 # if( Im(private$.innerprod(y0_, y1_)) > tol & mode(v0_) == "complex") stop("For complex vector methods only work in the horizontal case, yet.")
                                 if( method == "horizontal" & mode(v0_) == "numeric") 
                                   stop("horizontal method only works for the horizontal complex case")
                                 switch(method,
                                        horizontal = { 
                                          y0_ <- y0_ / sqrt(private$.innerprod(y0_))
                                          y1_ <- y1_ / sqrt(private$.innerprod(y1_))
                                          return( 
                                            v0_ - complex(re = 0, im = 1) * 
                                              Im( private$.innerprod(y1_, v0_) ) * (y0_ + y1_) / 
                                                    ( 1 + Mod(private$.innerprod(y0_, y1_)) ) )
                                        })
                               },
                               #' @description default plotting function for planar shapes, plotting \code{y_}
                               #' in front of \code{y0_} (after alignment).
                               #' @param col,pch,type graphical parameters passed to \code{base::plot} referring to \code{y_}.
                               #' @param ylab,xlab,xlim,ylim,xaxt,yaxt,asp,bty graphical parameters passed to \code{base::plot} 
                               #' with modified defaults.
                               #' @param ... other arguments passed to \code{base::plot}.
                               #' @param y0_par graphical parameters for \code{y0_}.
                               #' @param seg_par graphical parameters for line segments connecting \code{y_} 
                               #' and \code{y0_}.
                               plot = function(y_ = self$y_, y0_ = self$pole_, 
                                               yaxt='s', bty='n', ...) {
                                 super$plot(y_ = y_, y0_ = y0_, 
                                            yaxt=yaxt, bty=bty, ...)
                               },
                               #' @description Obtain "design matrix" of tangent space normal vectors in 
                               #' unstructured long format.
                               #' @param weighted logical, should inner product weights be pre-multiplied to 
                               #' normal vectors?
                               get_normal = function(y0_ = self$pole_, weighted = FALSE) {
                                 if(is.null(private$unstructure_index)) 
                                   stop("Please, initialize geometry before using $get_normal.")
                                 one <- rep(1, length(y0_))
                                 if(weighted) {
                                   y0_ <- private$.weights_ * y0_
                                   one <- private$.weights_ * one
                                 }
                                 if(!is.null(y0_)) 
                                   cbind(c(Im(y0_), -Re(y0_)), 
                                         complex2realmat(one))[private$unstructure_index, ]
                               }
                             ), private = list(
                               .distance = function(y0_, y1_, squared = FALSE) {
                                 d2 <- private$.innerprod(y0_ - private$.align(y1_, y0_))
                                 if(squared)
                                   return(d2)
                                 sqrt(d2)
                               },
                               .align = function(y_, y0_) {
                                 ip <- c(
                                   y_ = private$.innerprod(y_),
                                   y0_ = private$.innerprod(y0_)
                                 )
                                 if(any(ip==0))
                                   return(y_) 
                                 rot <- private$.innerprod(y_, y0_) / sqrt(prod(ip))
                                 return( rot * y_ / Mod(rot) )
                               }
                             ))



# Product geometries ------------------------------------------------------


#' A general product geometry
#'
#' An \code{mfGeometry} class implementing a product geometry in an abstract
#' general way. While there might be more efficient implementations using the 
#' particular structure of your geometry of interest, this should work for
#' products of any mfGeometries.
#'
#' @param y a list of `y` elements of the component manifolds in an 
#' unstructured 'flattened' way
#' @param v a list of `v` vectors representing tangent vectors in the component
#' manifolds in an unstructured 'flattened' way
#' @param y_ a list of objects living in the component manifolds
#' @param y0_ a list of objects living in the component manifolds
#' @param y1_ a list of objects living in the component manifolds
#' @param v_ a list of tangent vectors living in the component manifolds
#' @param v0_ a list of tangent vectors living in the component manifolds
#' @param v1_ a list of tangent vectors living in the component manifolds
#' @param ... other arguments passed to underlying \code{mfGeometry}.
#' 
#' @field mfGeom_default an \code{mfGeometry} object implementing the default 
#' geometry of a component of the product space.
#' 
#' @field y_ return lists of respective elements in the components.
#' \code{y_} can be specified as a list of \code{mfGeometry} objects constituting 
#' the components of the product geometry.
#' @field pole_ return lists of respective elements in the components.
#' @field weights_ return lists of respective elements in the components.
#' @field mfGeom_classes return a list of the classes of \code{y_} (read only).
#'
# #' @name mfGeometry
# #' @rdname mfGeometry
#' @import R6
#' @export
#' 
#' @rdname mfGeomProduct
mfGeomProduct <- R6Class("mfGeomProduct", inherit = mfGeometry, 
                      public = list(
                        mfGeom_default = NA,
                        #' @description Initialize product geometry with each 
                        #' component the \code{mfGeom_default} based on the 
                        #' \code{data} containing the variables for \code{y_} 
                        #' specified in the \code{formula}.
                        #' @param mfGeom_default \code{mfGeometry} object carrying the
                        #' component geometry.
                        #' @param formula formula containing interpreted 
                        #' via \code{mfInterpret_formula}) and indicating, in particular,
                        #' the ID variable for splitting observations in the product 
                        #' components.
                        #' @param data containing required variables.
                        #' 
                        initialize = function(mfGeom_default, data, formula) {
                          if(!missing(mfGeom_default)) 
                            self$mfGeom_default <- mfGeom_default
                          if(missing(data) & missing(formula)) 
                            return(invisible(self))
                          
                          if(!inherits(self$mfGeom_default, "mfGeometry"))
                            stop("Default geometry $mfGeom_default has to be set 
                                   before initliaizing data.")
                          
                          v <- mfInterpret_objformula(formula)
                          
                          if(is.matrix(data[[v$value]])) {
                            ## regular FDboost data format
                            # take first data point as template
                            dat1 <- data[names(data) != v$value]
                            dat1[[v$value]] <- data[[v$value]][1, ]
                            self$mfGeom_default$initialize(data = dat1, 
                                                           formula = formula)
                            # multiplicate first to set structure
                            n <- nrow(data[[v$value]])
                            
                            # fill in true values
                            private$.y_ <- lapply(1L:n, function(i) {
                              e <- self$mfGeom_default$clone()
                              e$y_ <- self$mfGeom_default$register(
                                self$mfGeom_default$structure(data[[v$value]][i, ]))
                              e
                            })
                            # prepare skeleton for reconstruction
                            private$y_skeleton <- asplit(
                              matrix(1L:length(data[[v$value]]), nrow = n), 1)
                            
                            names(private$.y_) <- data[[v$id]]
                            
                            # only for $get_normal() and $unstructure_weights()
                            private$regular_design <- TRUE
                          } else {
                            ## irregular FDboost data format
                            data <- as.data.frame(data[
                              sapply(data, length) == length(data[[v$id]])])
                            
                            datalist <- split(data, data[[v$id]])
                            
                            private$.y_ <- lapply(
                              datalist,
                              function(data) {
                                e <- self$mfGeom_default$clone()
                                e$initialize(data = data, formula = formula)
                                e
                              })
        
                            # prepare skeleton for reconstruction
                            private$y_skeleton <- split(1L:length(data[[v$value]]), data[[v$id]])
                            
                            names(private$.y_) <- sort(unique(data[[v$id]]))
                          }
                          # make sure order is preserved
                          index <- self$structure(1L:length(data[[v$value]]))
                          private$unstructure_index <- order(private$.unstructure(index))
                        },
                        #' @description subset product geometry to specified components.
                        #' Make sure to clone the product geometry first if 
                        #' complete geometry should be preserved.
                        #' @param which integers or character strings indicating
                        #' which elements of \code{y_} should be selected. 
                        #'  
                        slice = function(which) {
                          private$.y_ <- private$.y_[which]
                          y_skeleton <- private$y_skeleton[which]
                          private$y_skeleton <- relist(
                            as.numeric(as.factor(unlist(y_skeleton))), y_skeleton)
                          # make sure order is preserved
                          index <- self$structure(1L:length(unlist(y_skeleton)))
                          private$unstructure_index <- order(private$.unstructure(index))
                          invisible(self)
                        },
                        #' @description bring \code{y} into the right order to
                        #' split and pass it to the component geometries.
                        structure = function(y) {
                          if(is.null(private$y_skeleton))
                            stop("Please, initialize geometry before using $structure.")
                          private$.structure(y[unlist(private$y_skeleton)])
                        },
                        #' @description unstructure \code{y_} in the component geometries
                        #' and arrange it into the right order.
                        unstructure = function(y_) {
                          if(is.null(private$unstructure_index))
                            stop("Please, initialize geometry before using $unstructure.")
                          private$.unstructure(y_)[private$unstructure_index]
                        },
                        #' @description unstructure \code{weights_} in the component geometries
                        #' and arrange it into the right order.
                        unstructure_weights = function(weights_) {
                          uw <- Map(function(e, weights_) e$unstructure_weights(weights_), 
                                    private$.y_, weights_)
                          if(private$regular_design) {
                            uw[[1]]
                            } else {
                              c(unlist(uw))[private$unstructure_index]
                            }
                          },
                        #' @description bring \code{weights} into the right order to
                        #' split and pass it to the component geometries.
                        structure_weights = function(weights) {
                          self$structure(weights)
                        },
                        #' @description component-wise alignment.
                        align = function(y_, y0_, ...) {
                          mapply( function(e, y_, y0_) e$align(y_, y0_, ...), 
                                  private$.y_, y_, y0_, SIMPLIFY = FALSE ) 
                        },
                        #' @description component-wise registration.
                        register = function(y_, ...) {
                          mapply( function(e, y_) e$register(y_, ...), 
                                  private$.y_, y_, SIMPLIFY = FALSE ) 
                        },
                        #' @description component-wise registration.
                        register_v = function(v_, y0_, ...) {
                          mapply( function(e, v_, y0_) e$register_v(v_, y0_, ...), 
                                  private$.y_, v_, y0_, SIMPLIFY = FALSE ) 
                        },
                        #' @description compute vector of distances in component
                        #' geometries.
                        distance = function(y0_, y1_, ...) {
                          mapply( function(e, y0_, y1_) e$distance(y0_, y1_, ...), 
                                  private$.y_, y0_, y1_, SIMPLIFY = TRUE )
                        }, 
                        #' @description loop over Exp-maps in the component
                        #' geometries.
                        exp = function(v_, y0_, ...) {
                          mapply( function(e, v_, y0_) e$exp(v_, y0_, ...), 
                                  private$.y_, v_, y0_, SIMPLIFY = FALSE )
                        }, 
                        #' @description loop over Log-maps in the component
                        #' geometries.
                        log = function(y_, y0_, ...) {
                          mapply( function(e, y_, y0_) e$log(y_, y0_, ...), 
                                  private$.y_, y_, y0_, SIMPLIFY = FALSE )
                        },
                        #' @description loop over parallel transports in the component
                        #' geometries.
                        transport = function(v0_, y0_, y1_, ...) {
                          mapply( function(e, v0_, y0_, y1_) e$transport(v0_, y0_, y1_, ...), 
                                  private$.y_, v0_, y0_, y1_, SIMPLIFY = FALSE )
                        },
                        #' @description compute vector of inner products in the 
                        #' component geometries. An inner product on the product
                        #' space can be obtained as a scalar product of the 
                        #' returned vector.
                        innerprod = function(v0_, v1_ = v0_, ...) {
                          mapply( function(e, v0_, v1_) e$innerprod(v0_, v1_, ...), 
                                  private$.y_, v0_, v1_, SIMPLIFY = FALSE )
                        },
                        #' @description loop over individual plots of the components.
                        plot = function(y_ = self$y_, y0_ = self$pole_, main = names(y_), ...) {
                          stopifnot(is.list(y0_)|is.null(y0_))
                          stopifnot(is.list(y_))
                          if(is.null(y0_))
                            y0_ <- list(NULL)[rep(1, length(y_))]
                          warn <- TRUE
                          for(i in seq_along(y_)) {
                            if(is.null(private$.y_[[i]]$plot)) {
                              if(warn) {
                                warning("No plotting method for geometry specified.")
                                warn <- FALSE
                              }
                            } else {
                              private$.y_[[i]]$plot(y_ = y_[[i]], 
                                                    y0_ = y0_[[i]], main = main[i], ...)
                            }
                          }
                        },
                        get_normal = function(y0_ = self$pole_, weighted = FALSE) {
                          if(private$regular_design) {
                            if(!all(sapply(y0_, identical, y0_[[1]])))
                              warning("$get_normal for regular data, but not all elements of y0_ are identical.
                                      Only the first one will be used for all.")
                            #only return normal vector for first component
                            nvecs <- private$.y_[[1]]$get_normal(y0_ = y0_[[1]], weighted = weighted)
                          } else {
                            nvecs <- mapply( function(e, y0_) e$get_normal(y0_, weighted = weighted),
                                             private$.y_, y0_, SIMPLIFY = FALSE)
                            nvecs <- do.call(rbind, nvecs)[private$unstructure_index, ]
                          }
                          nvecs 
                        },
                        validate = function(y0_) { 
                          if(is.null(private$.y_)) 
                            is.list(y0_) else 
                          mapply(function(e, y0_) e$validate(y0_),
                                 private$.y_, y0_)
                        y0_
                        }
                      ),
                      private = list(
                        deep_clone = function(name, value) {
                          if(name == ".y_") {
                            # also clone R6 objects
                            return(
                              lapply(value, function(y) {
                              if(inherits(y, "R6"))
                                y$clone(deep = TRUE) else
                                  y
                            })
                            )
                          }
                          if(name == "mfGeom_default") {
                            if(inherits(value, "R6"))
                              return(value$clone(deep = TRUE))
                          }
                            # if no special case
                            value
                        },
                        regular_design = FALSE,
                        y_skeleton = NULL,
                        unstructure_index = NULL,
                        .unstructure = function(y_) c(unlist(
                          Map(function(e, y_) e$unstructure(y_), 
                              private$.y_, y_))),
                        .structure = function(y) {
                          Map(function(e, y) e$structure(y),
                              private$.y_, 
                              relist(y, private$y_skeleton))
                        }
                      ),
                      active = list(
                        y_ = function(value) {
                          if(missing(value))
                            lapply(private$.y_, function(e) e$y_) else {
                              private$.y_ <- Map(function(e, .y_) {
                                if(inherits(e, "mfGeometry")) return(e) else {
                                  if(is.null(.y_)) {
                                    if(!inherits(self$mfGeom_default, "mfGeometry")) 
                                      stop("Product geometry component is not specified as
                                         mfGeometry, but no mfGeometry in $mfGeom_default.")
                                    e_ <- self$mfGeom_default$clone()
                                    e_$y_ <- e
                                    e_
                                  } else {
                                    .y_$y_ <- e
                                    .y_
                                  }
                                }
                              }, value, private$.y_)
                              }
                            },
                        pole_ = function(value) {
                          if(missing(value)) {
                            if(is.null(private$.y_)) return(NULL)
                            lapply(private$.y_, function(e) e$pole_)
                          } else {
                              for(i in seq_along(value)) private$.y_[[i]]$pole_ <- value[[i]]
                            } },
                        weights_ = function(value) {
                          if(missing(value)) {
                            if(is.null(private$.y_)) return(NULL)
                            lapply(private$.y_, function(e) e$weights_)
                          } else {
                              if(is.null(private$.y_)) stop("Please initialize $y_ first.")
                              for(i in seq_along(value)) private$.y_[[i]]$weights_ <- value[[i]]
                            } },
                        mfGeom_classes = function(value) {
                          if(missing(value)) {
                            if(is.null(private$.y_)) 
                              stop("Set $y_ first please.")
                            lapply(y_, class)
                          } else {
                            stop("This is read only.")
                          }
                        })    
) #mfGeomProduct

# 
# # Product geometries regular --------------------------------------------
# 
# #' @export
# #' @rdname mfGeom_make_product
# mfGeom_make_powering <- function( inherit, classname = paste0(inherit$classname, "_reg"),
#                                   public = NULL, private = NULL, active = NULL) {
#   stopifnot(is.R6Class(inherit))
#   ## check whether is a mfGeometry class generator
#   # function from wch at https://stackoverflow.com/questions/37303552/r-r6-get-full-class-name-from-r6generator-object
#   findClasses <- function(x) {
#     if (is.null(x))
#       return(NULL)
#     parent <- x$get_inherit()
#     c(x$classname, findClasses(parent))
#   }
#   if(!("mfGeometry" %in% findClasses(inherit))) 
#     stop("'inherit' is not inheriting from mfGeometry.")
#   
#   # set up default public slots
#   public_ <- list(
#     structure = function(y, y0_ = private$.y_, v, ...) {
#       if( missing(y) + missing(v) != 1) 
#         stop("One and only one of y and v have to be supplied.")
#       
#       n <- if(is.null(dim(y0_))) 
#         length(y0_) else
#           nrow(y0_)
#       
#       if(missing(v)) {
#         y <- matrix(y, nrow = n)
#         y <- lapply(seq_len(n), function(i) y[i, ])
#         mapply(super$structure, 
#                y = y,
#                y0_ = y0_, MoreArgs = list(...), 
#                SIMPLIFY = FALSE)
#       } else {
#         v <- matrix(v, nrow = n)
#         v <- lapply(seq_len(n), function(i) v[i, ])
#         mapply(super$structure, 
#                v = v, 
#                y0_ = y0_, MoreArgs = list(...),
#                SIMPLIFY = FALSE)
#       }
#       
#     },
#     unstructure = function(y_, v_, weights_, y0_ = private$.y_, ...) {
#       if(missing(y_) + missing(v_) + missing(weights_) != 2)
#       stop("One and only one of 'y_', 'v_' and 'weights_' has to be provided.") 
#       if(is.null(y0_)) stop("Please specify a reference y0_ or a pole beforehand.")
#       
#       if(!missing(v_)) ret <- unlist(mapply(super$unstructure, v_ = v_, y0_ = y0_, 
#                                      MoreArgs = ...)) else {
#           if(!missing(weights_)) ret <- unlist(mapply(super$unstructure, 
#                                                          weights_ = weights_, 
#                                                          y0_ = y0_, MoreArgs = ...)) else
#                   ret <- unlist(mapply(super$unstructure, y_ = y_, y0_ = y0_, 
#                                                                 MoreArgs = ...))
#                                      }
#       n <- length(y0_)
#       ret <- matrix(unlist(ret), nrow = n, byrow = TRUE)
#       c(ret)
#     }, 
#     align = function(y_, y0_, ...) {
#       mapply( super$align, y_ = y_, y0_ = y0_, 
#               MoreArgs = list(...), SIMPLIFY = FALSE ) 
#     }, 
#     distance = function(y0_, y1_, ...) {
#       mapply( super$distance, y0_ = y0_, y1_ = y1_, 
#               MoreArgs = list(...), SIMPLIFY = TRUE )
#     }, 
#     exp = function(v_, y0_, ...) {
#       mapply( super$exp, v_ = v_, y0_ = y0_, 
#               MoreArgs = list(...), SIMPLIFY = FALSE ) 
#     }, 
#     log = function(y_, y0_, ...) {
#       mapply( super$log, y_ = y_, y0_ = y0_, 
#               MoreArgs = list(...), SIMPLIFY = FALSE ) 
#     }, 
#     transport = function(v0_, y0_, y1_, ...) {
#       mapply( super$transport, v0_ = v0_, y0_ = y0_, y1_ = y1_, 
#               MoreArgs = list(...), SIMPLIFY = FALSE )
#     },
#     innerprod = function(v0_, v1_ = v0_) {
#       if(is.null(weights_)) 
#         mapply( super$innerprod, v0_ = v0_, v1_ = v1_, 
#                 SIMPLIFY = TRUE ) else
#                   mapply( super$innerprod, v0_ = v0_, v1_ = v1_, 
#                           SIMPLIFY = TRUE )
#     },
#     get_normal = function(y0_ = private$.pole_, ...) {
#       nvecs <- lapply( y0_, super$get_normal ) 
#       matrix(unlist(nvecs), ncol = ncol(nvecs[[1]]))
#     },
#     validate = function(y0_, ...) {
#       stopifnot(is.list(y0_))
#       mapply( super$validate, y0_ = y0_, 
#               MoreArgs = list(...) )
#       if(is.null(attr(y0_, "mf.id"))) {
#         # generate an .id vector for the product
#         prod.id <- paste0(".id.", seq_along(y0_), "_")
#         temy0_ <- y0_
#         names(temy0_) <- prod.id
#         prod.id <- substr(prod.id, 1, nchar(prod.id)-1)
#         get_id <- function(a) regmatches(a, regexpr("^([^_]+)", a))
#         attr(y0_, "mf.id") <- factor(get_id(names(self$unstructure(y_ = temy0_, y0_ = temy0_))),
#                                      levels = prod.id)
#       }
#       y0_
#     }
#     
#   )
#   
#   # modify public slots
#   public_[names(public)] <- public
#   
#   # return new R6 class generator
#   R6Class(classname = classname, 
#           inherit = inherit, 
#           public = public_, private = private, 
#           active = active
#   )
# } # mfGeom_make_powering
# 
# # Create Euclidean family for multiple (irreg.) vectors -------------------
# 
# # #' @name mfGeomEuclidean_reg
# # #' @rdname mfGeometry-class
# # #' @export
# mfGeomEuclidean_reg <- mfGeom_make_powering(mfGeomEuclidean, 
#                                             "mfGeomEuclidean_reg")
# 
# # Create general planar shape family for multiple (irreg.) shapes ---------
# 
# # #' @name mfGeomPlanarShapes_reg
# # #' @rdname mfGeometry-class
# # #' @export
# mfGeomPlanarShapes_reg <- mfGeom_make_powering(mfGeomPlanarShape,
#                                                "mfGeomPlanarShapes_reg")