
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
ShapeReg <- function(pole) {
  
  if(!inherits(pole, "shape_long")) pole <- as_shape_long(pole)
  pole_cplx <- as_shape_complex(pole) 
  long.id <- pole$.id
  p <- with(pole, split(.value, .id))
  p_cplx <- with(pole_cplx, split(.value, .id))
    
  response <- function(f){
    return( unlist( mapply( Exp, split(f, long.id), p, SIMPLIFY = FALSE) ) )
  }
  
  loss <- function(y, f) {
    return( unlist( mapply(function(x,y) geodist(x,y)^2, 
                           split(y, long.id), split(f, long.id), SIMPLIFY = FALSE) ) )
  }
  
  # helper functions
  vec2cplx <- function(x) {
    k <- length(x)/2
    complex(re = head(x,k), im = tail(x,k))
  }
  cplx2vec <- function(x) {
    c(Re(x), Im(x))
  }
  proj <- function(x, p) {
    x - p * crossprod(Conj(p), x)
  }
  
  ngradient <- function(y, f, w = 1) {
    # project f into TpM
    f <- tapply(f, long.id, vec2cplx, simplify = FALSE)
    f <- mapply(proj, f, p_cplx, SIMPLIFY = FALSE)
    # calc predictions
    f_resp <- mapply(Exp, f, p_cplx, SIMPLIFY = FALSE) 
    y <- tapply(y, long.id, vec2cplx, simplify = FALSE)
    # Rotation align y and resp
    y <- mapply(rotate, x = y, x0 = f_resp, SIMPLIFY = FALSE )
    # Tangent vector in tangent space of respective y[i,]
    v <- mapply(Log, x = y, x0 = f_resp, SIMPLIFY = FALSE)
    # return tangent vectors parallel transported to offset
    # rotation align f_resp to p
    f_resp <- mapply(rotate, x = f_resp, x0 = p_cplx, SIMPLIFY = FALSE)
    # Parallel transport tangent vectors to TpM
    v <- mapply(Transport, v = v, x0 = f_resp, x1 = p_cplx, 
                MoreArgs = list(method = "horizontal"), SIMPLIFY = FALSE )
    # Return v after basis tranform and in appropriate format
    v <- unlist(lapply(v, cplx2vec), use.names = FALSE)
    return( v )
  }
  
  offset <- function(y, w) 0 #weighted.mean
  
  Family(ngradient = ngradient, loss = loss, response = response, weights = "any",
         offset = offset, name = "Planar Shape Regression with fixed base point")
}
