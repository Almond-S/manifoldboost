# function essentially copied from sparseFLMM:::smooth.construct.symm.smooth.spec
# only different constraint

smooth.construct.lc.smooth.spec <- function (object, data, knots) 
{
  if (length(object$term) != 2) 
    stop("basis only handels 2D smooths")
  x <- data[[object$term[1]]]
  y <- data[[object$term[2]]]
  if (length(unique(x)) < object$bs.dim) 
    warning("basis dimension is larger than number of unique covariates")
  if (is.null(object$xt)) 
    object$xt <- list(bsmargin = "ps", kroneckersum = TRUE)
  if (is.null(object$xt$kroneckersum)) 
    object$xt$kroneckersum <- TRUE
  if (is.null(object$xt$bsmargin)) 
    object$xt$bsmargin <- "ps"
  if (object$xt$bsmargin != "ps") 
    stop("marginal smooth class need to be 'ps'")
  if (length(object$p.order) == 1) {
    m <- rep(object$p.order, 2)
  }
  else {
    m <- object$p.order
  }
  m[is.na(m)] <- 2
  object$p.order <- m
  if (object$bs.dim < 0) 
    object$bs.dim <- max(10, m[1])
  nk <- object$bs.dim - m[1]
  if (nk <= 0) 
    stop("basis dimension too small for b-spline order")
  k1 <- knots[[object$term[1]]]
  k2 <- knots[[object$term[2]]]
  if (!is.null(k1) & !is.null(k2)) {
    if ((k1 != k2)) 
      stop("number of specified knots is not equal for both margins")
  }
  object$line.con <- object$point.con
  object$point.con <- NULL #list(0, 0)
  # names(object$point.con) <- object$term
  
  Sm <- list()
  # create marginal smooths
  object$margin <- lapply(object$term, function(ter) {
    smooth1 <- smooth.construct(eval(as.call(list(as.symbol("s"), 
                                       as.symbol(ter), bs = object$xt$bsmargin, 
                                       pc = if(!is.na(object$line.con[[ter]])) object$line.con[[ter]],
                                       k = object$bs.dim, m = object$p.order))), 
                     data = data, 
                     knots = knots)
    # apply line constraint(s)
    if(!is.na(object$line.con[[ter]])) {
      # for centering
      dat0 <- data.frame(t0 = object$line.con[[ter]])
      smooth0 <- suppressWarnings(smooth.construct(eval(as.call(list(as.symbol("s"), 
                                                as.symbol("t0"), bs = object$xt$bsmargin, 
                                                k = object$bs.dim, m = object$p.order))),
                                                data = dat0, knots = list(t0 = smooth1$knots)))
      # adapted form FDboost:::X_bbsc
      C <- c(smooth0$X)
      qr_C <- qr(C)
      Q <- qr.Q(qr_C, complete = TRUE)
      Z <- Q[, (2):ncol(Q)]
    
      smooth1$X <- smooth1$X %*% Z
      smooth1$S[[1]] <- crossprod(Z, smooth1$S[[1]]) %*% Z
      smooth1$bs.dim <- smooth1$bs.dim - 1
      smooth1$D <- smooth1$D %*% Z
      r <- qr(smooth1$S[[1]])$rank
      nsd <- nrow(smooth1$S[[1]]) - r
      smooth1$rank <- r
      smooth1$null.space.dim <- nsd
      smooth1$Z <- Z
      # smooth1$side.constrain <- TRUE
    }
    smooth1
  })
  
  X <- tensor.prod.model.matrix(X = list(object$margin[[1]]$X, object$margin[[2]]$X))
  Sm[[1]] <- object$margin[[1]]$S[[1]]
  Sm[[2]] <- object$margin[[2]]$S[[1]]
  if (object$xt$kroneckersum) {
    S <- tensor.prod.penalties(list(Sm[[1]], Sm[[2]]))
    S <- S[[1]] + S[[2]]
  }
  else {
    S <- Sm[[1]] %x% Sm[[2]]
  }
  r <- qr(S)$rank
  nsd <- nrow(S) - r
  object$S <- list(S)
  object$X <- X
  object$rank <- r
  object$null.space.dim <- nsd
  object$m <- m
  object$knots <- k1
  class(object) <- "lc.smooth"
  # avoid fürther identifiability constraints
  object$C <- matrix(0, 0, ncol(object$X)) 
  object$side.constrain <- FALSE
  object
}


# predict function --------------------------------------------------------

Predict.matrix.lc.smooth <- function (object, data) 
{
  m <- length(object$margin)
  X <- list()
  for (i in 1:m) {
    term <- object$margin[[i]]$term
    dat <- list()
    for (j in 1:length(term)) {
      dat[[term[j]]] <- data[[term[j]]]
    }
    X[[i]] <- PredictMat(object$margin[[i]], dat, n = length(dat[[1]]))
    if(!is.null(object$margin[[i]]$Z)) {
      X[[i]] <- X[[i]] %*% object$margin[[i]]$Z
    }
  }
  tensor.prod.model.matrix(X)
}


