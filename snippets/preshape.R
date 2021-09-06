
# classical preshapes ---------------------------------------------------

#' Make preshapes from shape configurations
#' 
#' @description 
#' 
#' `preshape()` transforms shape configurations into preshapes
#' for further processing and analysis.
#' 
#' @param x The shape object. If `x` is a `complex` vector or 
#' `numeric` vector or matrix, it will be considered as one single shape;
#' entire shape datasets should be provided as `shape_long`/`shape_wide`
#' object or as an `array`.
#' @param translate Logical, should translation be removed.
#' @param scale Logical, should scale be removed.
#' @param type The type of the preshape. Currently, only 'center'
#' is possible.
#'  
#' @export
#' @example R/preshape_example.R   
#' 
preshape <- function(x, translate = TRUE, scale = TRUE, int.weights = NULL, arg.range = NULL, A = NULL, type = "center") {
  UseMethod("preshape")
}

preshape.complex <- function(x, translate = TRUE, scale = TRUE, A = NULL, type = "center") {
  type <- match.arg(type)
  switch(type,
              center = {
                if(translate) {
                  ONE <- rep(1, length(x))
                  x <- x - innerprod(ONE, x, A = A) * ONE / innerprod(ONE, A = A)
                }
                if(scale) x <- x / sqrt(innerprod(x, A = A))
              }
  )
  return(x)
}

preshape.matrix <- function(x, translate = TRUE, scale = TRUE, A = NULL, type = "center") {
  type <- match.arg(type)
  switch(type,
              center = {
                if(translate) {
                  ONE <- rep(1, len = nrow(x))
                  x <- x - apply(x, 2, function(x) innerprod(ONE, x, A = A) * ONE / innerprod(ONE, A = A))
                }
                if(scale) x <- x / sqrt(innerprod(x, A = A))
              }
  )
  return(x)
}

preshape.array <- function(x, translate = TRUE, scale = TRUE, A = NULL, type = "center") {
  type <- match.arg(type)
  switch(type,
         center = {
           if(translate) {
             if(is.null(A)) x <- sweep(x, 2:3, apply(x, 2:3, mean), `-`)
             if(!is.null(A) & is.null(dim(A))) x <- sweep(x, 2:3, apply(x, 2:3, weighted.mean, w = A), `-`)
             if(!is.null(A) & !is.null(dim(A))) x <- sweep(x, 2:3, apply(x, 2:3, weighted.mean, w = colSums(A)), `-`)
           }
           if(scale) x <- sweep(x, 3, 
                                apply(x, 3, function(x) sqrt(innerprod(x, A = A))), 
                                `/`)
         }
  )
  return(x)
}

preshape.shape_long <- function(x, translate = TRUE, scale = TRUE, int.weights = NULL, arg.range = NULL, type = "center") {
  as_shape_long( 
    preshape.shape_complex( as_shape_complex(x), 
              translate = translate, scale = scale, 
              int.weights = int.weights, arg.range = arg.range, 
              type = type))
}

preshape.shape_complex <- function(x, translate = TRUE, scale = TRUE, int.weights = NULL, arg.range = NULL, type = "center") {
  if(is.null(int.weights)) {switch(type, 
         center = {
           value <- split(x$.value, x$.id)
           pre <- sapply(value, preshape.complex, 
                         translate = translate, scale = scale, type = type)
           x$.value <- c(unlist(pre))
         })
  return(x)}
  if(is.null(arg.range)) arg.range <- range(x$.arg)
  switch(type, 
         center = {
           pre <- lapply(split(x, x$.id), function(x) {
             A <- switch(int.weights,
                         empirical = NULL,
                         trapezoidal = trapez_weights(x$.arg, arg.range),
                         simpson = simpson_weights(x$.arg, arg.range) )
             preshape.complex( x$.value, translate = translate, scale = scale, 
                              A = A, type = type)
           })
           x$.value <- c(unlist(pre))
         })
  return(x)
}

preshape.shape_wide <- function(x, translate = TRUE, scale = TRUE, int.weights = NULL, arg.range = NULL, type = "center") {
  if(is.null(int.weights)) {switch(type, 
         center = {
           value <- split(x[, c(".value1", ".value2")], x$.id)
           pre <- lapply(value, function(x) preshape.matrix(as.matrix(x), 
                         translate = translate, scale = scale, type = type) )
           x[, c(".value1", ".value2")] <- do.call(rbind, pre)
         })
  return(x)
  }
  if(is.null(arg.range)) arg.range <- range(x$.arg)
  switch(type, 
         center = {
           pre <- lapply(split(x, x$.id), function(x) {
             A <- switch(int.weights,
                         empirical = NULL,
                         trapezoidal = trapez_weights(x$.arg, arg.range),
                         simpson = simpson_weights(x$.arg, arg.range) )
             preshape.matrix( as.matrix(x[, c(".value1", ".value2")]), translate = translate, scale = scale, 
                               A = A, type = type)
           })
           x[, c(".value1", ".value2")] <- do.call(rbind, pre)
         })
  return(x)
}

