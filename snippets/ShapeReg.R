
# ShapeReg mboost Family --------------------------------------------------

# Let y a numeric vector containing the preshapes and structured as produced by as_long() in BaseFunctions.R, i.e. in particular with an attribute 'Dim' containing c(n.landmarks, n.dimensions, n.individuals)
# Correspondingly, y has to be structured such that array(y, dim = Dim(x)) is an array of the corresponding design.
# p is the pole for the response, typically p = meanShape(y)

### TODO: specify Dim(y) as a ShapeReg-Argument, because sadly, attributes of y are lost somewhere

#' Shape regression family
#' 
#' @description 
#' 
#' Family for Shape Regression for planar shapes and a fix 'intercept' p.
#' 
#' @param p The shape pole of the tangent space the linear predictor lives in. 
#' Can be understood as an offset.
#'  
#' @export
#' @examples
#' #TODO: Add examples.
#' 
ShapeReg <- function(p) {
  
  #TODO: check for data being balanced
  
  p_cplx <- as_complex(p)
  k <- length(p_cplx)
  
  response <- function(f){
    f <- matrix(f, nrow = 2*k)
    return( apply( f, 2, Exp, x0 = p) ) 
  }
  
  loss <- function(y, f) {
    dim(f) <- dim(y) <- c( 2*k, length(y)/(2*k) )
    return( vapply(1:ncol(y), function(i) geodist(y[, i], f[, i] )^2, numeric(1)) )
  }
  
  ngradient <- function(y, f, w = 1) {
    dimy <- c(k,2, length(y)/(2*k))
    f_resp <- response(f)
    dim(f_resp) <- dim(y) <- dimy
    f_resp <- as_complex(f_resp)
    y <- as_complex(y)
    # Rotation align y and resp
    y <- vapply(1:ncol(y), function(i) rotate(x = y[, i], x0 = f_resp[, i]), complex( l = dimy[1] ) )
    # Tangent vector in tangent space of respective y[i,]
    v <- vapply(1:ncol(y), function(i) Log(x = y[, i], x0 = f_resp[, i]), complex( l = dimy[1] ) )
    # return tangent vectors parallel transported to offset
    # rotation align y to p
    f_resp <- apply(f_resp, 2, rotate, x0 = p_cplx)
    # Parallel transport tangent vectors to TpM
    v <- vapply(1:ncol(y), function(i) Transport(v[, i], x0 = f_resp[, i], x1 = p_cplx, method = "horizontal"), complex(dimy[1]) )
    # Return v after basis tranform and in appropriate n.individuals x (n.landmarks-1) format
    v <- structure( c(Re(v), Im(v)), dim = c(dim(v), 2) )
    return( aperm(v, c(1,3,2)) )
  }
  
  offset <- weighted.mean
  
  Family(ngradient = ngradient, loss = loss, response = response, weights = "none",
         offset = offset, name = "Planar Shape Regression with fixed base point")
}
