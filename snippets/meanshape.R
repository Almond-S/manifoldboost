# compute mean shape ------------------------------------------------------

#' Compute shape mean
#' 
#' @description 
#' 
#' Function computing a Full Procrustres mean of a collection of shapes.
#' 
#' 
#' @param x A shape dataset provided as `shape_long`/`shape_wide`/`shape_complex` object
#' or as a k x 2 x n `array` (corresponding to the format in the R package shapes).
#' @param distance Underlying distance for mean computation. So far, only `"fullprocrustres"` is available.
#' @param type Type of shape space. Defaults to `"center"`.
#' @param return_preshape Logical, if `TRUE` (default) the mean is returned as preshape; 
#' if `FALSE` it is returned as shape configuration.
#'  
#' @export
#' @example R/meanshape_example.R   
#' 

meanshape <- function(x, distance = "fullprocrustres", type = c("center"), return_preshape = TRUE) {
  UseMethod("meanshape")
}

meanshape.matrix <- function(x, distance = "fullprocrustres", type = c("center"), return_preshape = TRUE) {
  if(mode(x) != "complex") stop("If x is a matrix, mode(x) must be complex.")
  x <- preshape(x, type = type)
  S <- matrix( rowSums( apply(x, 2, function(x) tcrossprod(x, Conj(x))) ), nrow = nrow(x))
  ei <- eigen( S )$vectors[,1]
  if(!return_preshape) ei <- crossprod( Conj(shapes::defh( length(ei) )), ei )
  return( rotate(ei, x0 = x[,1]) )
}

meanshape.default <- function(x, distance = "fullprocrustres", type = c("center"), return_preshape = TRUE) {
  type <- match.arg(type)
    x <- as_complex(x)
    xmean <- meanshape(x, distance = distance, type = type, return_preshape = return_preshape)
    xmean <- cbind(Re(xmean), Im(xmean)) 
    rownames(xmean) <- rownames(x)
  xmean
}


# B-spline mean shapes ----------------------------------------------------

# depends on packages FDboost / mboost, Matrix

bshape <- function(x,
                   bfun = bbsc, smoothed.cov = FALSE, make.preshapes = TRUE, 
                   cyclic = FALSE, differences = c(0, 2), lambda0 = 0.1, cv.seed = NULL, 
                   cv.k = NULL, cov.k = -1, int.weights = NULL, arg.grid = NULL, 
                   return.grid = 200, return.aligned = TRUE, ...) {
  lambda0 <- matrix(lambda0, nrow = length(differences), ncol = 1)
  if(make.preshapes) 
    x <- preshape( as_shape_complex(x), int.weights = int.weights ) else
      x <- as_shape_complex(x)
  if(is.null(arg.grid)) {
    arg.grid <- unique(x$.arg)
  }
  
  # compute bspline normalization ... 
  Bfun <- bfun(arg.grid, differences = differences[1], cyclic = cyclic, ...)
  B <- extract(Bfun, "design")
  int.range <- environment(Bfun$dpp)$args$knots[[1]]$boundary.knots
  # compute basis integration weights
  if(is.null(int.weights)) A <- 1/length(arg.grid) else {
    A <- switch (int.weights,
            empirical = {
              if(is.null(arg.grid)) {
                freq <- table(x$.arg)
                as.numeric(freq) / length(freq) } else 
                  freq <- density(x$.arg)
                approx(freq$x, freq$y, arg.grid)$y 
                },
            trapezoidal = trapez_weights(arg.grid, range = int.range),
            simpson = simpson_weights(arg.grid, range = int.range)
    )
  }
  if(is.null(dim(A))) 
    R <- crossprod(B, A * B) else {
    R <- crossprod(B, A %*% B)
    }
  R <- as.matrix(R)
  
  if(!smoothed.cov) {
    # perform penalized PCA approach similar to the one described by 
    # Silverman (1996): smoothed functional principal component analysis by choice of norm or 
    # by Ramsay and Silverman (2012)
    
    # ...and penalization matrices
    P <- vapply(1:length(differences), function(i) {
      if(i == 1) return(as.matrix(extract(Bfun, "penalty")))
      as.matrix(
        extract(bfun(arg.grid, differences = differences[i], cyclic = cyclic, ...), "penalty"))
    }, numeric(length = ncol(B)^2 ), 
    USE.NAMES = FALSE)
    
    # compute transformed covariance operator
    if(is.null(int.weights)) {
      By <- lapply(split(x, x$.id), 
                   function(x) {
                     crossprod( as.matrix(
                       extract(bfun(x$.arg, differences = differences[1], cyclic = cyclic, ...), "design")), 
                       x$.value )
                   })
    } else { 
      By <- switch(int.weights,
                   empirical = lapply(split(x, x$.id), 
                                      function(x) {
                                        crossprod( as.matrix(
                                          extract(bfun(x$.arg, differences = differences[1], cyclic = cyclic, ...), "design")), 
                                          x$.value )
                                      }),
                   trapezoidal = lapply(split(x, x$.id), 
                                        function(x) {
                                          crossprod( as.matrix(
                                            extract(bfun(x$.arg, differences = differences[1], cylic = cyclic, ...), "design")), 
                                            trapez_weights(x$.arg, range = int.range) * x$.value )
                                        }),
                   simpson = lapply(split(x, x$.id), 
                                    function(x) {
                                      crossprod( as.matrix(
                                        extract(bfun(x$.arg, differences = differences[1], cyclic = cyclic, ...), "design")), 
                                        simpson_weights(x$.arg, range = int.range) %*% x$.value )
                                    })
      )
    }
    
    
    ByyB <- sapply(By, function(x) tcrossprod(x, Conj(x)))
    
    calc_pca <- function(lambdas, which_fold) {
      
      # compute penalty matrix
      P <- matrix( P %*% lambdas, ncol = ncol(B) )
      
      sumByyB <- matrix( rowSums( ByyB[, which_fold] ), ncol = ncol(B) )
      
      # compute penalized mean coefficients
      if(qr(R+P)$rank < nrow(R)) 
        stop("Singular basis covariance matrix:\n increase regularization param lambda0 or try (longer) arg.grid.")
      eigen(as.matrix(solve(R + P)) %*% sumByyB)
    }
    
    # set up shape-wise cross-validation folds
    if(!is.null(cv.k)) {
      if(!is.null(cv.seed)) set.seed(cv.seed)
      ids <- unique(x$.id)
      cv_folds <- matrix(0, nrow = cv.k, ncol = length(ids))
      cv_folds[1:length(ids)] <- sample(ids, length(ids))
      
      cv_risk <- function(lambdas) {
        cv_loss <- apply(cv_folds, 1, function(fo) {
          which_fold <- ids %in% ids[fo]
          
          mean_estimate <- calc_pca(lambdas, !which_fold)$vectors[,1]
          sumByyB <- matrix( rowSums( ByyB[, which_fold] ), ncol = ncol(B) )
          
          Mod( crossprod(Conj(mean_estimate), sumByyB) %*% mean_estimate  /
                 crossprod(Conj(mean_estimate), R + matrix( P %*% lambdas, ncol = ncol(B) ) ) %*% mean_estimate )
        })
        sum(cv_loss)
      }
      
      lambda0 <- optim(par = lambda0, fn = cv_risk)$par
    }
    
    # perform pca with optimal lambda0
    ei <- calc_pca(lambda0, TRUE)
  }
  
  if(smoothed.cov) {
    # use a PACE-type procedure generalizing the 
    # Cederbaum, Scheipl, Greven (2016): Fast symmetric additive covariance smoothing
    # to complex covariance matrices
    require(mgcv)
    
    # we model the part of E[Conj(y(t1))y(t2)] where t1 < t2 
    # assuming that the time points are ordered
    cov_dat <- lapply(split(x, x$.id), function(x) {
                          combs <- combn(1:nrow(x), 2)
                          with(x, data.frame(
                            yyt = .value[combs[1,]] * Conj(.value[combs[2,]]), 
                            arg1 = .arg[combs[1,]],
                            arg2 = .arg[combs[2,]]
                          ))
                        })
    cov_dat <- do.call(rbind, cov_dat)
    cov_fit_re <- bam( Re(yyt) ~ s(arg1, arg2, bs = "sps", k = cov.k, 
                                   xt = list(cyclic = cyclic)), 
                       data = cov_dat, knots = list(arg1 = int.range) )
    cov_fit_im <- bam( Im(yyt) ~ -1 + s(arg1, arg2, bs = "sps", k = cov.k,  
                                        xt = list(skew = TRUE, cyclic = cyclic)), 
                       data = cov_dat, knots = list(arg1 = int.range ) )
    # predict smoothed covariance
    cov_dat <- expand.grid(arg1 = arg.grid, arg2 = arg.grid)
    yyt <- matrix( complex(
      real = predict(cov_fit_re, newdata = cov_dat),
      imaginary = predict(cov_fit_im, newdata = cov_dat) ), 
      ncol = length(arg.grid))
      
    if(is.null(dim(A))) 
               sumByyB <- crossprod( A * B, yyt) %*% (A * B) else {
                 sumByyB <- crossprod( A %*% B, yyt) %*% A %*% B 
               }
    
    if(qr(R)$rank < nrow(R)) 
      stop("Singular basis covariance matrix:\n try (longer) arg.grid.")
    ei <- eigen(as.matrix(solve(R)) %*% sumByyB)
  }
  
  #### prepare results 
  
  btangent <- bfun(x$.arg, differences = differences[1], cyclic = cyclic, ...)
  
  # function converting complex matrix operators to vecctorized real analoga
  cmat2vmat <- function(A) {
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

  # transform and vectorize design & penalty matrix
  e <- environment(btangent$dpp)
  # id-wise trafo 
  e$X <- do.call(rbind, lapply(
    split(as.data.frame(as.matrix(e$X) %*% ei$vectors[,-1]), x$.id), 
    cmat2vmat))
  e$K <- cmat2vmat(
    crossprod( ei$vectors[,-1]
               , as.matrix(e$K) ) %*% ei$vectors[,-1]
  )
  
  mean_coefs <- ei$vectors[,1]
  
  if(return.aligned) mean_coefs <- mean_coefs / innerprod(mean_coefs, A = R)
  
  if(is.null(return.grid) | return.aligned) {
    shp_mean <- c(as.matrix(
      extract(bfun(x$.arg, differences = differences[1], cyclic = cyclic, ...), "design")) %*% mean_coefs)
    if(return.aligned) {
      rot <- c( crossprod(Conj(shp_mean), x$.value) / crossprod(Conj(shp_mean), shp_mean) )
    }
    if(is.null(return.grid)) {
      x$.value <- if(return.aligned) rot * shp_mean else shp_mean
      
      ret <- list( meanshape = as_shape_long(x), btangent = btangent )
    }
  } 
  # return on grid
  if(length(return.grid)==1) return.grid <- seq(int.range[1], int.range[2], len = return.grid) 
  
  if(!is.null(return.grid)) ret <- list( 
    meanshape = as_shape_long( shape_complex(
      value = ifelse(return.aligned, rot, 1) * c(as.matrix(
        extract(bfun(return.grid, differences = differences[1], cyclic = cyclic, ...), "design")) %*% mean_coefs) , 
      arg = return.grid,
      id = "mean_fun")
    ), 
    btangent = btangent )
  
  ret
}





# REML utilities for specification of penalization ------------------------

## Probably not needed anymore

# normalization constant of the complex Bingham distribution (for )
c_Bingham <- function(A) {
  ei <- eigen(A)
  s <- vapply(1:nrow(A), function(j) {
    aj <- prod(eigen(A)$value[j] - eigen(A)$value[-j])
    if(any(aj==0)) stop("A can only have distinct eigenvalues.\n
                        One other case is implemented in c_Watson.")
    aj * exp(eigen(A)$value[j])
             }, 1)
  2*pi^nrow(A)*sum(s)
}

c_Watson <- function(kappa, k) {
  2*pi^(k-1)*kappa^(2-k)*exp(kappa)* 
    (1- sum(sapply(0:(k-3), function(r) 
      kappa^r * exp(-kappa) / gamma(r+1))))
}



