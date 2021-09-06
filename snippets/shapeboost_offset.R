
shapeboost_offset <- function(
  formula,
  dim = NULL,
  id = NULL,
  data,
  cyclic = FALSE,
  smoothed.cov = NULL,
  smoothed.cov_control = NULL,
  weights = NULL,
  numInt = "equal"
) {
  
  ### check formulas
  stopifnot(class(formula) == "formula")
  if(class(try(id)) == "try-error") stop("id must either be NULL or a formula object.")
  if(class(try(dim)) == "try-error") stop("dim must either be NULL or a formula object.")
  if(!is.null(id)) stopifnot(class(id) == "formula")
  if(!is.null(dim)) stopifnot(class(dim) == "formula")
  
  ### initialize base-learners
  g2 <- Gaussian()
  g2@check_y <- function(y) y # just to avoid any errors
  g2@offset <- function(y, weights) 0
  mod0 <- mboost(formula, data, 
                 family = g2, control = boost_control(mstop = 0))
  if(length(mod0$baselearner) != 1) 
    stop("Exactly one base-learner has to be specified for offset.")
  bl <- mod0$baselearner[[1]]
  
  ### extract id if supplied
  if(!is.null(id)) {
    stopifnot(class(id) == "formula")
    nameid <- all.vars(id)[[1]]
    stopifnot(length(nameid) == 1)
    .id <- data[[nameid]]
    data[[nameid]] <- NULL
  } else {
    nameid <- NULL
    .id <- NULL
  }
  
  ### extract dim if supplied
  if(!is.null(dim)) {
    stopifnot(class(dim) == "formula")
    ydimvar <- all.vars(dim)[[1]]
    nameydim <- ydimvar
    stopifnot(length(ydimvar) == 1)
    .dim <- data[[ydimvar]]
    if(length(unique(.dim))!=2) 
      stop("Only planar 2-dimensional shapes provided, yet.")
    data[[ydimvar]] <- NULL
  } else {
    nameydim <- NULL
    .dim <- NULL
  }
  
  ### prepare response
  response <- mod0$response
  if(!is.array(response)) {
    ## basic checks
    if(!is.null(dim(response))) 
      stop("If not supplied as an array, the response has to be a vector.")
    if(is.null(id)) 
      stop("If response is not an array, id has to be specified.")
    if(length(.id) != length(response))
      stop("id must be of same length as response.")
  }
  
  ### prepare arg
  nameyarg <- bl$get_names()
  .arg <- data[[nameyarg]]
  if(is.array(response))
    stopifnot(length(.arg) == dim(response)[1]) else
      stopifnot(length(.arg) == length(response))
  if(!is.null(weights)) stopifnot(length(weights) == length(.arg))
  arg.range <- environment(bl$dpp)$args$knots[[1]]$boundary.knots
  # make arg.grid for numerical integrations
  arg.grid <- if(is.array(response)) .arg else
    seq(arg.range[1], arg.range[2], len = 40)
  ### make numerical integration weights for arg.grid
  ### function setting up integration weights
  make_integration_weights <- function(numInt, arg = NULL, range = NULL) {
    switch (numInt,
            equal = NULL,
            trapezoidal = trapez_weights(arg, range = arg),
            simpson = simpson_weights(arg, range = arg)
    )
  } 
  A <- make_integration_weights(numInt, arg.grid, arg.range)
  
  ### Built shape dataset
  if(is.array(response)) dresponse <- as_shape_long.array(
    x = response, 
    arg = .arg, 
    id = .id) else
      dresponse <- shape_long(value = response, arg = .arg, id = .id)
  
  B <- extract(bl, "design")
  
  if(is.null(smoothed.cov)) {
    smoothed.cov <- !is.array(response) 
  }
  
  ### get estimate of ByyB, the complex covariance matrix of the basis coefs
  ### without covariance smoothing (regular / irregular case) or with cov smooth
  if(!smoothed.cov) {
    # could be extended to regularized FPCA
    # a la Silverman (1996): smoothed functional principal component analysis 
    # by choice of norm
    # but didn't work so well and is, thus, not implemented here
    # => if penalization necessary always use smoothed.cov
    
    
    ## function returning ByyB
    compute_ByyB <- function(value, B, arg = NULL, range = NULL, int.weights = NULL) {
      yy <- tcrossprod(value, Conj(value))
      ByyB <- if(is.null(int.weights))
        crossprod(B, yy) %*% B else {
          if(is.null(dim(int.weights)))
            crossprod(int.weights*B, yy) %*% (int.weights*B) else
              crossprod(int.weights %*% B, yy) %*% int.weights %*% B
        }
      ByyB
    }
    
    if(is.array(response)) {
      ByyB <- compute_ByyB(as_complex(response), B, .arg, arg.range, A)
        } else { 
      warning("For irregular data it is often not advisable to do 
              non-regularized PCA. You might want to set smoothed.cov = TRUE.")
      ## combine relevant info in one matrix
      B_ <- cbind(response, .arg, B)
      B_ <- split(B_, .id)
      ByyB <- sapply(B_, function(B_id) {
        A <- make_integration_weights(numInt, 
                                      arg = B_id[1:nrow(B_id), 2], arg.range)
        compute_ByyB(
          value = complex(matrix(B_id[,1], ncol = 2)),
          B = B_id[1:nrow(B_id), -(1:2)], 
          arg = B_id[1:nrow(B_id), 2], arg.range, numInt)
      })
      ByyB <- matrix( rowSums(ByyB), ncol = ncol(B) )
    }
  } else {
    # -> smoothed.cov == TRUE
    # use a PACE-type procedure generalizing the 
    # Cederbaum, Scheipl, Greven (2016): Fast symmetric additive covariance smoothing
    # to complex covariance matrices
    require(mgcv)
    
    # we model the part of E[Conj(y(t1))y(t2)] where t1 < t2 
    # assuming that the time points are ordered
    cresponse <- as_shape_complex(dresponse)
    cov_dat <- lapply(split(cresponse, cresponse$.id), function(x) {
      combs <- combn(1:nrow(x), 2)
      with(x, data.frame(
        yy = .value[combs[1,]] * Conj(.value[combs[2,]]), 
        arg1 = .arg[combs[1,]],
        arg2 = .arg[combs[2,]]
      ))
    })
    cov_dat <- do.call(rbind, cov_dat)
    
    cov_fit_re <- bam( Re(yy) ~ s(arg1, arg2, bs = "sps", k = cov.k, 
                                   xt = list(cyclic = cyclic)), 
                       data = cov_dat, knots = list(arg1 = int.range) )
    cov_fit_im <- bam( Im(yy) ~ -1 + s(arg1, arg2, bs = "sps", k = cov.k,  
                                        xt = list(skew = TRUE, cyclic = cyclic)), 
                       data = cov_dat, knots = list(arg1 = int.range ) )
    # predict smoothed covariance
    cov_dat <- expand.grid(arg1 = arg.grid, arg2 = arg.grid)
    yy <- matrix( complex(
      real = predict(cov_fit_re, newdata = cov_dat),
      imaginary = predict(cov_fit_im, newdata = cov_dat) ), 
      ncol = length(arg.grid))
    
    ByyB <- if(is.null(A))
      crossprod(B, yy) %*% B else {
        if(is.null(dim(A)))
          crossprod(A*B, yy) %*% (A*B) else
            crossprod(A %*% B, yy) %*% A %*% B
    }
  }
  
  ### compute matrix R = BB
  if(is.array(response)) {
    if(is.null(A)) 
      R <- crossprod(B) else {
        if(is.null(dim(A))) 
          R <- crossprod(B, A * B) else {
            R <- crossprod(B, A %*% B)
          }
      }
    R <- as.matrix(R)
  } else { # irregular case
    newdat <- data.frame(.arg = arg.grid)
    names(newdat) <- nameyarg
    # evaluate B on integration grid 
    B.grid <- with(newdat, 
                   extract(eval(as.expression(bl$get_call())), "design")
                   )
    if(is.null(A)) 
      R <- crossprod(B.grid) else {
        if(is.null(dim(A))) 
          R <- crossprod(B.grid, A * B.grid) else {
            R <- crossprod(B.grid, A %*% B.grid)
          }
      }
    R <- as.matrix(R)
  }
  
  ### Compute leading eigenvector
  ei <- eigen(as.matrix(solve(R)) %*% ByyB)
  
  ### return offset (and offset coefs)
  ret <- list(
    offset_coefs = ei$vectors[,1],
    offset = numeric(length = ncol(B)),
    offset_tangent_coefs = ei$vectors[,-1]
  )
  if(is.array(response)) {
    ret$offset <- B %*% cbind(Re(ret$offset_coefs), Im(ret$offset_coefs))
  } else {
    dim1 <- .dim == unique(.dim)[1]
    ret$offset[dim1] <- B[dim1, ] %*% Re(ret$offset_coefs)
    ret$offset[!dim1] <- B[!dim1, ] %*% Im(ret$offset_coefs)
  }
  
  ret
}
