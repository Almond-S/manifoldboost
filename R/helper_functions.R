

# function transforming complex linear operator matrix to real lin --------

#' Convert complex to real liner operator matrix 
#' 
#' @description 
#' 
#' `complex2realmat()` converts a complex \eqn{n\times m}{(n x m)}-matrix \eqn{A} 
#' representing a 
#' complex linear operator in \eqn{GL(\mathds{C}^m, \mathds{C}^n)}{GL(C^m, C^n)} to a
#' real \eqn{2n\times 2m}{(2n x 2m)} matrix \eqn{A'} representing the corresponding 
#' real linear 
#' operator in \eqn{GL(\mathds{R}^{2m}, \mathds{R}^{2n})}{GL(R^2m, R^2n)}, 
#' identifying \eqn{z \in \mathds{C}^n}{z in C^n} with \eqn{(Re(z^\top), Im(z^\top))^\top}{(Re(z),Im(z)}.
#' 
#' @param A a `complex` (or `numeric`) matrix.
#'  
#'
complex2realmat <- function(A) {
  A <- as.matrix(A)
  colnames(A) <- NULL
  n <- nrow(A)
  m <- ncol(A)
  A_ <- matrix(nrow = 2*n, ncol = 2*m)
  A_[1:n, 1:m] <- A_[n + 1:n, m + 1:m] <- Re(A)
  A_[n + 1:n, 1:m] <- Im(A)
  A_[1:n, m + 1:m] <- -Im(A)
  A_
}


# casting two dims into complex -------------------------------------------

#' Converting to complex by casting dim into real/imaginary
#'
#' @param x an object of class `matrix` or `array`.
#' @param .margin the dimension to interpret as real and imaginary numbers
#'
#' @return complex vector, matrix, or array
#'
as_complex <- function(x, .margin = 2) {
  stopifnot(length(dim(x)[.margin]) == 2)
  UseMethod()
}

as_complex.matrix <- function(x, .margin = 2) {
  if(.margin == 1) 
    complex(real = x[1, ], imag = x[2, ]) else
        complex(real = x[, 1], imag = x[, 2])
}

as_complex.array <- function(x, .margin = 2) {
  .dim <- dim(x)
  .dimnames <- dimnames(x)
  x <- matrix(
    aperm(x, c(setdiff(seq_len(dim(x)), .margin), .margin)),
    ncol = .margin )
  x <- as_complex.matrix(x, .margin = 2)
  array(x, dim = .dim, dimnames = .dimnames)
}


# Trapezoidal rule --------------------------------------------------------

#' Trapezoidal rule numerical integration weights
#' 
#' @description 
#' 
#' Constructs a vector with numerical integration weights for one single function.
#' 
#' @param arg A numeric vector with 'time-points' (argument values) the function is evaluated at. 
#' Can be in non-increasing order.
#' @param range the range of the interval integrated over. Must include all arg values.
#' 
#' @export

trapez_weights <- function(arg, range = NULL) {
  if(is.null(range)) range <- range(arg) else {
    if(any(arg < range[1] | arg > range[2])) 
      stop("All arg values have to lie within the range.") }
  t_diffs <- diff( c(range[1], sort(arg), range[2]) )
  t_diffs[-c(1,length(t_diffs))] <- t_diffs[-c(1,length(t_diffs))]/2
  weights_sorted <- t_diffs[-1] + t_diffs[-length(t_diffs)]
  # weights_sorted <- 1/2 * diff( c(range[1], sort(arg), range[2]), lag = 2 )
  return(weights_sorted[order(order(arg))])
}

simpson_weights <- function(arg, range = NULL, sparseMat = FALSE) {
  if(is.null(range)) range <- range(arg)
  M <- bandSparse(n = length(arg), k = c(0,1), symmetric = TRUE, 
                  diagonals = list("0" = 2/3*trapez_weights(arg = arg, range = range),
                                   "1" = 1/3*diff(sort(arg))[order(order(arg))] ))
  if(sparseMat) return(M)
  as.matrix(M)
}

# Simple equal weighting --------------------------------------------------------

#' Simple 1/n_grid weights 
#' 
#' @description 
#' 
#' Constructs a vector with equal weights of 1/n_grid instead of numerical 
#' integration weights for one single function.
#' 
#' @param arg A numeric vector with 'time-points' (argument values) the function is evaluated at. 
#' Can be in non-increasing order.
#' @param range the range of the interval integrated over. Must include all arg values.
#' 
#' @export

equal_weights <- function(arg, range = NULL) 
  rep(1, length(arg)) / length(arg)


