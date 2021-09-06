
# IrregularShapeReg mboost Family --------------------------------------------------

#' Shape regression family
#' 
#' @description 
#' 
#' Family for Shape Regression for planar shapes and a fixed 'intercept', the pole.
#' 
#' @param pole an vactor with shape pole in long format,
#' where N is total amount of measurements over all response curves.
#' Can be understood as an offset.
#' @param param a vector of length N containing the 'time-points' of the measurements of all curves, 
#' i.e. the argument of the observed curves. Numeric for functional shapes, factor for landmark shapes  
#' @param id a factor vector of length N (or coerced to one) with the IDs identifying measurements 
#' belonging to the same shape.
#'  
#' @export
#' @examples
#' #TODO: Add examples.
#' 
IrrShapeReg <- function( pole, pole.arg = NULL, pole.dim = NULL, pole.id = NULL, pole.arg.range = NULL) {
  
  if(!"shape_long" %in% class(pole)) pole <- shape_long(pole, pole.arg, pole.dim, pole.id)
  if(is.null(pole.arg.range) & is.numeric(pole$.arg)) pole.arg.range <- range(pole$.arg)
  
  p <- with(pole, split(.value, .id))
  
  # generate integration matrix / vector
  if(is.factor(pole$.arg)) A <- NULL else {
    # trapezoidal rule
    A <- with(pole, tapply(.arg, list(.id, .dim), trapez_weights, pole.arg.range))
    A <- mapply(c, A[,1], A[,2], SIMPLIFY = FALSE)
  }
  
  response <- function(f){
    f <- split(f, pole$.id)
    c(unlist(mapply( Exp, f, p, A)))
  }
  
  loss <- function(y, f) {
    f <- split(f, pole$.id)
    f_resp <- mapply( Exp, f, p, A)
    f_resp <- lapply(f_resp, function(x) as_complex(matrix(x, ncol = 2)))
    y <- lapply( split(y, pole$.id), function(x) as_complex(matrix(x, ncol = 2)))
    return( mapply(geodist, y, f_resp, A[length(unique(pole$.id))]) )
  }
  
  ngradient_one <- function(y, f, A, p) {
    f_resp <- Exp(f, p, A)
    f_resp_c <- as_complex(matrix(f_resp, ncol = 2))
    y_c <- as_complex(matrix(y, ncol = 2))
    y_c <- rotate(y_c, x0 = f_resp_c)
    y <- c(Re(y_c), Im(y_c))
    
    norm_f <- sqrt(innerprod(f))
    if(norm_f>1e-5) {
      # original formula
      dmu <- sin(norm_f)/norm_f * ( - tcrossprod(p, f) + diag(nrow = length(f)) ) +
      (norm_f * cos(norm_f) - sin(norm_f)) / norm_f^3 * tcrossprod(f) } else {
        # Taylor series expensions at 0
        dmu <- (1 - norm_f/6 + norm_f^4/120) * ( - tcrossprod(p, f) + diag(nrow = length(f)) ) +
          ( -1/3 + norm_f^2/20 - norm_f^4/840 ) * tcrossprod(f)
      }
    
    geodist(f_resp, y, A) / sqrt(1 - innerprod(f_resp, y, A)) * innerprod(y, dmu, A)
  }
  
  ngradient <- function(y, f, w = 1) {
    
    f <- split(f, pole$.id)
    y <- split(y, pole$.id)
    
    unlist(mapply(ngradient_one, y, f, A, p))
    
  }
  
  offset <- function(y, w = 1) 0
  
  Family(ngradient = ngradient, loss = loss, response = response, weights = "none",
         offset = offset, name = "Planar Shape Boosting with fixed base point")
}
