
#' Periodic / centered B-splines   
#' 
#' @description 
#' 
#' `hspline()` implements B-splines with the option of removing the intercept from the 
#' span of the basis and/or applying a periodicity constraint.
#' 
#' @param x B-spline argument.
#' @param knots if of length 1, the number of equally spaced inner knots, 
#' but also a vector of inner knots can be supplied. (For an unconstrained 
#' B-spline, the number of the B-spline basis functions is then given by the 
#' number of inner knots + degree + 1.)
#' @param boundary.knots the boundary knots of the B-spline. Defaults to 
#' range(x). 
#' @param degree the B-spline degree, i.e. the degree of the piece-wise 
#' polynomials; one less than the B-spline order. 
#' @param remove.constant logical, should the constant function be removed 
#' from the span of the B-spline functions. Induces an integral to zero constraint.
#' @param periodic logical, should the B-splines be constrained to be periodic.
#' @param derivs return the `derivs`th derivative of the B-spline basis functions.
#' Can be supplied as a vector of length(x). See `splines::splineDesign`.
#' @param sparse return design matrix as sparse matrix. See `splines::splineDesign`.
#'
#' @details 
#' NOTE: due to machine accuracy setting boundary.knots = range(x) can lead to 
#' an spurious error claiming x is outside the boundary knots. In this case, it
#' is advised to substract/add a very small constant to them.
#'

hspline <- function(x, knots = 10, boundary.knots = NULL, degree = 3, 
                    remove.constant = FALSE, periodic = FALSE, derivs = 0, sparse = FALSE)
{
  if(is.null(boundary.knots)) {
    if(length(x) == 1) stop("If only one x value supplied, boundary.knots must
                            be given.")
    boundary.knots <- range(x) 
  } else {
    if(length(boundary.knots) != 2) stop("boundary.knots has to be of length 2.")
  }
   
  # set up homogeneous basis by continuing the knots outside of the boundary
  if(length(knots) == 1) {
    h <- diff(boundary.knots)/(knots+1) 
    knots <- seq(min(boundary.knots) - degree*h, 
               max(boundary.knots) + degree*h, len = knots+2+2*degree)
  } else {
    knots <- sort(knots)
    knots <- c( tail(knots, degree) - diff(range(boundary.knots)), 
                min(boundary.knots), knots, max(boundary.knots), 
                head(knots, degree) + diff(range(boundary.knots)))
  }
  
  B <- splines::splineDesign(knots = knots, x = x, ord = degree+1, 
                             derivs = derivs, sparse = sparse, 
                             outer.ok = FALSE) # error if x ouside of boundary.knots
  
  # apply constraints if supplied
  if(periodic) {
    # identify the spline functions on right boundary with those on the left
    B[, 1:degree] <- B[, 1:degree] + B[, ncol(B) + 1 - degree:1]
    B <- B[, 1:(ncol(B)-degree), drop = FALSE]
  }
  
  if(remove.constant) {
    # motivated from implementation in FDboost:::X_bbsc
    # for homogeneous B-splines the constant function corresponds to constant
    # basis coefficients. This is also pertains for periodic B-splines.
    # Thus, we remove the constant function from the B-spline-space by
    # using a coefficient basis orthogonal to the constant coefficient vector.
    C <- rep(1, ncol(B))
    # use QR decomposition only to obtain desired orthonormal basis Q
    qr_C <- qr(C) 
    Q <- qr.Q(qr_C, complete=TRUE)  # orthonormal matrix of QR decomposition
    Z <- Q[ , 2:ncol(Q), drop = FALSE ] # first element is multiple of C itself
    B <- B %*% Z
  }
  
  B
}


