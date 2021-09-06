


# inner product computation -------------------------------------------------------------

#' Compute an inner product
#' 
#' @description 
#' 
#' `innerprod()` computes the inner product of two `complex` or `numeric` vectors.
#' 
#' @param x,y `complex` or `numeric` vectors.
#'  
#' @export
#' 
#'
innerprod <- function(x, y, A = NULL) {
  UseMethod("innerprod")
}

innerprod.complex <- function(x, y = x, A = NULL) {
  if(is.null(A)) return(as.complex(sum(Conj(x)*y)))
  if(is.null(dim(A)))  return(as.complex(sum(A*Conj(x)*y)))
    return( as.complex(crossprod(Conj(x),A)%*%y) )
}

innerprod.numeric <- function(x, y = x, A = NULL) {
  if(is.null(A)) return(sum(x*y))
  if(is.null(dim(A)))  return(sum(A*x*y))
  return( sum(crossprod(x,A)%*%y) )
}

# rotation operators ------------------------------------------------------

#' Rotate shape configurations
#' 
#' @description 
#' 
#' Rotate shape either by some angle alpha, which might be
#' the angle aligning it to a template
#' 
#' @param x The shape object(s) to be rotated. Might be a 
#' single shape given as `complex` vector or `numeric` matrix, 
#' or an entire shape dataset provided as
#' a k x 2 x n `array`.
#' @param theta The rotation angle provided in radiance.
#' @param x0 a template shape `x` is rotation aligned to 
#' if `x0` is provided instead of `theta`. 
#'  
#' @export
#' @example R/preshape_example.R   
#' 
rotate <- function(x, theta, x0, ...) {
  if(!missing(theta) & !missing(x0)) 
    stop("Either theta or x0 have to be specified.")
  UseMethod("rotate")
}

rotate.complex <- function(x = 1, theta, x0, A = NULL, type.theta = c("radial", "degree")) {
  if(missing(x0)) {
    return( x * exp(complex(im = 1) * theta) )
  }
  
  if(missing(theta)) {
    rot <- innerprod(x, x0, A = A)
    return( rot * x / Mod(rot) )
  }
} 

rotate.matrix <- function(x, theta, x0) {
  dimx <- dim(x)
  if(mode(x) == "numeric") {
        if(missing(x0)) {
        rotationMatrix <- if(ncol(x) == 2) {
          matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), nrow = 2, ncol = 2) 
        } else stop( "Only implemented for dim(x)[2], the dimension of the values, is 2." )
        return(x%*%rotationMatrix)
      }
      
      if(missing(theta)) {
        stopifnot(dim(x0) == dim(x))
        SVD <- La.svd(crossprod(xi, x0))
        rot <- SVD$u %*% SVD$vt
        return( x %*% rot )
      }
  }
  if(mode(x) == "complex") stop("For complex matrices not implemented, yet.")
  if(!(mode(x) %in% c("numeric", "complex"))) stop("mode(x) must either be numeric or complex.")
}

rotate.array <- function(x, theta, x0) {
  dimx <- dim(x)
  if(missing(x0)) {
    rotationMatrix <- if(dim(x)[2] == 2) {
      matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), nrow = 2, ncol = 2) 
    } else stop( "Only implemented for dim(x)[2], the dimension of the values, is 2." )
    x <- apply(x, 3, function(x) x%*%rotationMatrix)
    return(structure(x, dim = dimx))
  }
  
  if(missing(theta)) {
    stopifnot(dim(x0) == dim(x)[1:2])
    rot <- apply(x, 3, function(xi) {
      SVD <- La.svd(crossprod(xi, x0))
      SVD$u %*% SVD$vt } )
    dim(rot) <- rep(dim(x),c(0,2,1))
    return( structure(sapply(1:dim(x)[3], 
                             function(i) x[,,i]%*%rot[,,i] ), 
                      dim = dimx ) )
  }
}



# orthogonalization operator ----------------------------------------------

#' Orthogonalize a matrix with respect to another  
#' 
#' @description 
#' 
#' `orthogonalize()` determines a linear transformation matrix Z for a matrix x 
#' such that the columns of xZ are orthogonal to the columns of the matrix x0.
#' 
#' @param x matrix to orthogonalize.
#' @param x0 matrix with columns to be orthogonal to. 
#'  
#' 
#'
orthogonalize <- function(x, x0) {
  UseMethod("orthogonalize")
}

# Not tested in this context, yet:
orthogonalize.matrix <- function(x, x0) {
  # From FDboost %Xc% (respectively analogously):
  C <- crossprod( x, x0 )  # X1[ , -1], X2[ , -1])
  qr_C <- qr(C, tol = 1e-20)
  if( any(class(qr_C) == "sparseQR") ){
    rank_C <- qr_C@Dim[2]
  }else{
    rank_C <- qr_C$rank 
  } 
  Q <- qr.Q(qr_C, complete=TRUE)        # orthonormal matrix of QR decomposition
  Z <- Q[ , (rank_C+1):ncol(Q) ]        
  return(Z)
}


# compute geodesic distance -------------------------------------------------------

#' Geodesic distance
#' 
#' @description 
#' 
#' Function computing the geodesic distance between two shapes.
#' 
#' @param x,y Two shapes given as complex vector, configuration matrix or in long format.
#' @param type The type of shape considered.
#'  
#' @export
#' 
geodist <- function(x, y, A = NULL, type = c("center")) {
  type <- match.arg(type)
  switch(type,
    center = {
      ip <- innerprod(x,y, A = A)
      if(Mod(ip-1) < 1e-15) return(0)
      if(Mod(ip+1) < 1e-15) return(pi)
      acos( Mod(ip) ) })
}


# log mapping -------------------------------------------------------------

#' Manifold Log mapping on shape space.
#' 
#' @description 
#' 
#' Function computing the Log mapping of a shape into the tangent space at a reference shape.
#' 
#' @param x A shape given as complex vector, configuration matrix or in long vector form.
#' @param x0 Reference shape, the pole of the tangent space.
#' @param method Method for computing the Log mapping; one of 'simple' and 'general' (see details).
#'  
#' @details 
#' The $Log: M \mapsto T_{x_0}M$ maps a point $x \in M$ to a tangent vector $v \in T_{x_0}M$ keeping 
#' the length $||v||_2 = geodist(x_0, x)$. It is the inverse of the $Exp$-map and uniquely defined in
#' a neighbourhood of $x_0$.
#' Method `simple` corresponds widely to a formula given by Dryden and Mardia 2012 (p.76) with a small adjustment.
#' Method `general` is a formula used also by Cornea et al. 2017.
#'  
#' 
Log <- function(x, x0, method = c("simple", "alternative")) {
  method = match.arg(method)
  switch(method, 
         simple = { 
           rho <- geodist(x,x0)
           if(!is.na(rho)) if(rho < 1e-15) return(rep(0, length(x))) 
           return( rho / sin(rho) * ( x - x0 * innerprod(x0, x) )) ### See for example Dryden&Mardia2012, page 76 -> described for complex vectors, slightly adjusted with x0 * innerprod(x0, x) instead of x0 * Mod(innerprod(x0, x))
         },
         alternative = {
           v <- x - x0 * innerprod(x0, x)
           if(sqrt(innerprod(v)) == 0) return(v)
           else return( v * geodist(x, x0) / sqrt(innerprod(v)))  ### See for example Cornea et al. 2017
         })
}


# Exp mapping -------------------------------------------------------------

#' Manifold Exp mapping on shape space.
#' 
#' @description 
#' 
#' Function computing the Exp mapping of a tangent vector at a reference shape into the shape space.
#' 
#' @param x A shape given as complex vector, configuration matrix or in long vector form.
#' @param x0 Reference shape, the pole of the tangent space.
#' @param method Method for computing the Log mapping; one of 'simple' and 'general' (see details).
#'  
#' @details 
#' The $Exp: T_{x_0}M \mapsto M$ maps a tangent vector $v \in T_{x_0}M$ to a point $M \ni x = c(1)$
#' with $c: (0,1) \mapsto M$ a geodesic path with $c(0) = x_0$ and constant velocity $dc/dt = v$.
#'  
#' 

Exp <- function(v, x0, A = NULL) {
  if(innerprod(v) == 0) return(x0)
  
  s <- sqrt(innerprod(v, A = A))
  return( x0 * cos(s) + v / s * sin(s) )   ### see 
}


# Parallel transport ------------------------------------------------------

#' Parallel transport in shape space
#' 
#' @description 
#' 
#' Parallel transport of a tangent vector in one tangent space into another.
#' 
#' @param v A complex/numeric tangent vector in the tangent space at `x_0`.
#' @param x0 The reference shape given as complex vector, configuration matrix or in long vector form.
#' @param x1 Reference shape, the pole of the tangent space.
#' @param method Method for carrying out the parallel transport; 
#' one of 'general' and 'horizontal' (see details).
#' @param tol Tolerance when checking whether x0 and x1 are rotation aligned 
#' (i.e., the geodesic is horizontal). 
#'  
#' @details 
#' The $Exp: T_{x_0}M \mapsto M$ maps a tangent vector $v \in T_{x_0}M$ to a point $M \ni x = c(1)$
#' with $c: (0,1) \mapsto M$ a geodesic path with $c(0) = x_0$ and constant velocity $dc/dt = v$.
#' Method 'general' corresponds to the formulation in Cornea et al. 2017; 
#' it generally performs parallel transport in the preshape space. Method 'horizontal',
#' by contrast, implements the formula in Dryden and Mardia 2012 (p.76) performing parallel transport
#' along horizontal geodesics for complex vectors.
#'
Transport <- function(v, x0, x1, method = c("general", "horizontal", "simple"), tol = 1e-15) {
  method <- match.arg(method)
  if( Im(innerprod(x0, x1)) > tol & mode(v) == "complex") stop("For complex vector methods only work in the horizontal case, yet.")
  if( method == "horizontal" & mode(v) == "numeric") stop("horizontal method only works for the horizontal complex case")
  switch(method,
         simple = {
           x1tilde <- (x1 - innerprod(x0, x1)*x0) / sqrt(1- innerprod(x0, x1)^2)
           return( v - innerprod(x1tilde, v) * x1tilde - sqrt(1 - innerprod(x0, x1)^2) * innerprod(x1tilde, v) * x0 + Conj(innerprod(x0, x1))*innerprod(x1tilde, v) * x1tilde )
         },
         general = { ### See for example Cornea et al. 2017
           x1tilde <- (x1 - innerprod(x0, x1)*x0) / sqrt(1- innerprod(x0, x1)^2)
           return( v - innerprod(x0, v) * x0 - innerprod(x1tilde, v) * x1tilde + (innerprod(x0, x1)*innerprod(x0, v) - sqrt(1 - innerprod(x0, x1)^2) * innerprod(x1tilde, v)) * x0 + (sqrt(1- innerprod(x0,x1)^2) * innerprod(x0,v) + Conj(innerprod(x0, x1))*innerprod(x1tilde, v)) * x1tilde )
         }, 
         horizontal = { ### See for example Dryden&Mardia 2012 page 76, described for complex vectors
           return( v - ( innerprod(x1, v) * (x0 + x1) ) / ( 1 + Mod(innerprod(x0, x1)) ) )
         })
}

