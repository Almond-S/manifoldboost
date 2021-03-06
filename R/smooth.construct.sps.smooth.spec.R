# adapted from sparseFLMM::make_summation_matrix

#' Construct basis transformation matrix for (skew-)symmetry constraints
#' 
#' Construct basis transformation matrix for (skew-)symmetry constraints for 
#' bivariate P-spline smooths, which can also be combined with periodicity constraints
#' for cyclic marginals. The function is adapted from the 
#' \code{\link[sparseFLMM]{sparseFLMM}} with the same name.
#'
#' @param F number of marginal basis functions.
#' @param skew logical, should the basis be constraint to skew-symmetry instead 
#' of symmetry.
#' @param cyclic.degree integer, specifying the number of basis functions identified
#' with each other at the boundaries in order to implement periodicity. Should 
#' be specified to match the degree of the B-spline basis.
#' 
#' @details This function is used to implement the respective constraints in the
#' function \code{\link[smooth.construct.sps.smooth.spec]{smooth.construct.sps.smooth.spec}}.
#'
#' @return A basis transformation matrix of dimension \eqn{F^2 \times G} with 
#' \eqn{G<F^2} depending on the specified constraint.
#'
# #' @examples
make_summation_matrix <- function (F, skew = FALSE, cyclic.degree = 0) {
  
  ind_mat <- matrix(seq_len(F^2), ncol = F, nrow = F)
  pairs <- cbind(c(ind_mat), c(t(ind_mat)))
  C <- diag(F^2)
  cons <- pairs[pairs[, 1] < pairs[, 2], , drop = FALSE]
  
  if(skew) {
    C[, cons[, 1]] <- C[, cons[, 1]] - C[, cons[, 2]] 
  } else {
    C[, cons[, 1]] <- C[, cons[, 1]] + C[, cons[, 2]]
  }
  if(skew) ind_vec <- cons[,1] else 
    ind_vec <- pairs[pairs[, 1] <= pairs[, 2], 1, drop = FALSE]
  C <- C[, ind_vec]
  
  if(cyclic.degree>0) {
    if(F < 2*cyclic.degree) stop("For F<2*cyclic.degree not implemented, yet.")
    # helper function for mapping values of ind_mat to column idx of C
    match_ind <- function(ind) 
      sapply(ind, function(x) which.max(x == ind_vec))
    
    # match edges
    if(F > 2*cyclic.degree) {
      edge_pairs <- cbind(
        c(ind_mat[(cyclic.degree+1):(F-cyclic.degree), 1:cyclic.degree]), 
        c(t(ind_mat[F-cyclic.degree + 1:cyclic.degree, (cyclic.degree+1):(F-cyclic.degree)]))
        )
      edge_pairs <- apply(edge_pairs, 2, match_ind)
      if(skew) {
        C[, edge_pairs[, 1]] <- C[, edge_pairs[, 1]] - C[, edge_pairs[, 2]]
      } else {
        C[, edge_pairs[, 1]] <- C[, edge_pairs[, 1]] + C[, edge_pairs[, 2]]
      }
      C <- C[, -edge_pairs[, 2]]
      ind_vec <- ind_vec[-edge_pairs[, 2]]
    }
    # match corners
    corner_loc <- subset(expand.grid(row = 1:cyclic.degree, 
                                     col = 1:cyclic.degree), 
                         if(skew) row > col else row >= col)
    if(nrow(corner_loc)>0) {
      corner_pairs <- cbind(
        mapply(function(x,y) ind_mat[x,y], corner_loc$row + F-cyclic.degree, corner_loc$col),
        mapply(function(x,y) ind_mat[x,y], corner_loc$row, corner_loc$col),
        mapply(function(x,y) ind_mat[x,y], corner_loc$row + F-cyclic.degree, corner_loc$col + F-cyclic.degree)
      )
      
      corner_pairs <- matrix(apply(corner_pairs, 2, match_ind), ncol = ncol(corner_pairs))
      
      C[, corner_pairs[, 1]] <- C[, corner_pairs[, 1]] + 
        C[, corner_pairs[, 2]] + C[, corner_pairs[, 3]]
      C <- C[, -c(corner_pairs[, 2:3]), drop = FALSE]
      ind_vec <- ind_vec[-c(corner_pairs[, 2:3])]
    }
    
    # symmetrize lower corner square
    lower_square <- ind_mat[F-cyclic.degree + 1:cyclic.degree,
                            1:cyclic.degree ]
    if(length(lower_square)>1) {
      lower_pairs <- cbind(c(lower_square), c(t(lower_square)))
      lower_pairs <- lower_pairs[lower_pairs[, 1] < lower_pairs[, 2], , drop = FALSE]
      lower_pairs <- matrix(apply(lower_pairs, 2, match_ind), ncol = ncol(lower_pairs))
      if(skew) {
        C[, lower_pairs[, 1]] <- C[, lower_pairs[, 1]] - C[, lower_pairs[, 2]]
      } else { 
        C[, lower_pairs[, 1]] <- C[, lower_pairs[, 1]] + C[, lower_pairs[, 2]]
      }
      C <- C[,  -lower_pairs[,2], drop = FALSE]
      ind_vec <- ind_vec[-lower_pairs[,2]]
    }
    if(skew) C <- C[, -match_ind(diag(as.matrix(lower_square)))]
  }
  
  C
}


# construct (skew-)symmetic smooths ---------------------------------------

# adapted from sparseFLMM:::smooth.construct.symm.smooth.spec

#' (Skew-)Symmetric smooths constructor
#' 
#' The \code{sps} smoother class for \code{mgcv} is adapted from the \code{symm}
#' class in the R package \code{sparseFLMM}. Besides symmetric bivariate P-spline
#' tensor product smooths, it also offers the opportunity to specify skew-symmetric
#' smooths, the symmetry constraints also in combination with cyclic (periodic) 
#' marginal P-spline smooths, and the corresponding constraints for univariate 
#' smooths. 
#'
#' @param object a smooth specfication object, usually generated by a term of 
#' the form \code{s(x, bs = "sps", ...)}.
#' @param data a list containing the data.
#' @param knots a list containing any knots for basis setup. Can be \code{NULL}.
#'
#' @return An object of class "\code{sympspline.smooth}". See 
#' \code{\link[smooth.construct]{smooth.construct}} for its elements.
#' @importFrom mgcv smooth.construct
#' @export
#'
# #' @examples
smooth.construct.sps.smooth.spec <- function (object, data, knots) {
  
  if (length(object$term) > 2) 
    stop("basis only handels 1D and 2D smooths")
  
  if (is.null(object$xt)) 
    object$xt <- list(skew = FALSE, cyclic = FALSE)
  if(is.null(object$xt$skew))
    object$xt$skew <- FALSE
  if(is.null(object$xt$cyclic))
    object$xt$cyclic <- FALSE
  if(is.null(object$xt$bsmargin))
    object$xt$bsmargin <- "ps"
  
  # __ 1D case _____________________________________________________
  # determine designmat X, penaltymat S and Z transformation matrix 
  
  if(length(object$term) == 1) {
    
    if(object$xt$cyclic) 
      warning("Only 2D splines can be cross-cyclic.
              Hence, cyclic = TRUE is ignored.
              You might want to specify bsmargin = 'cp' instead
              to get cyclic B-splines.")
    
    # borrow form pspline smooth
    object <- smooth.construct(eval(as.call(list(as.symbol("s"), 
                                                 as.symbol(object$term[1]), 
                                                 bs = object$xt$bsmargin, 
                                                 pc = object$point.con, xt = object$xt,
                                                 k = object$bs.dim, m = object$p.order))), 
                               data = data, 
                               knots = knots)
    # make (skew)-symmetric coefficient basis
    if(object$xt$skew) {
      bs.dim <- floor(object$bs.dim/2)
      Z <- rbind( diag(nrow = bs.dim), 
                if(object$bs.dim %% 2) 0,
                - diag(nrow = bs.dim)[, bs.dim:1] )
    } else {
      bs.dim <- ceiling(object$bs.dim/2)
      Z <- rbind( diag(nrow = bs.dim),
                  diag(nrow = bs.dim)[,bs.dim:1] )
      if(object$bs.dim %% 2) Z <- Z[-bs.dim, ]
    }
    S <- object$S[[1]]
  } 
  
  # __ 2D case _____________________________________________________
  # determine designmat X, penaltymat S and Z transformation matrix 


  if(length(object$term) == 2) {
    
    x <- data[[object$term[1]]]
    y <- data[[object$term[2]]]
    if (length(unique(x)) < object$bs.dim) 
      warning("basis dimension is larger than number of unique covariates")
    if (is.null(object$xt$kroneckersum)) 
      object$xt$kroneckersum <- TRUE
    if (!all(sapply(object$xt$bsmargin, '%in%', c("ps", "cp")))) 
      stop("marginal smooth classes need to be 'ps' or 'cp'.")
    if(length(object$xt$bsmargin) == 1)
      object$xt$bsmargin <- rep(object$xt$bsmargin, 2)
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
    k1 <- if(is.null(knots[[object$term[1]]])) 
                 knots[[object$term[2]]] else knots[[object$term[1]]]
    k2 <- knots[[object$term[2]]]
    if(!is.null(k2)) {
      if(!identical(k1, k2)) 
        stop("number of specified knots is not equal for both margins")
    }
    if(is.null(k1)) k1 <- range(data[object$term])
    object$knots <- list(k1, k1) 
    names(object$knots) <- object$term
    
    Sm <- list()
    smooth1 <- smooth.construct(eval(as.call(list(as.symbol("s"), 
                                                  as.symbol(object$term[1]), bs = object$xt$bsmargin[1], 
                                                  k = object$bs.dim, m = object$p.order))), data = data, 
                                knots = object$knots[object$term[1]])
    smooth2 <- smooth.construct(eval(as.call(list(as.symbol("s"), 
                                                  as.symbol(object$term[2]), bs = object$xt$bsmargin[2], 
                                                  k = object$bs.dim, m = object$p.order))), data = data, 
                                knots = object$knots[object$term[2]])
    object$X <- tensor.prod.model.matrix(X = list(smooth1$X, smooth2$X))
    Sm[[1]] <- smooth1$S[[1]]
    Sm[[2]] <- smooth2$S[[1]]
    if (object$xt$kroneckersum) {
      S <- tensor.prod.penalties(list(Sm[[1]], Sm[[2]]))
      S <- S[[1]] + S[[2]]
    }
    else {
      S <- Sm[[1]] %x% Sm[[2]]
    }
    Z <- make_summation_matrix(F = object$bs.dim, 
                               skew = object$xt$skew,
                               cyclic.degree = object$xt$cyclic * (m[1]+1))
    object$margin < list()
    object$margin[[1]] <- smooth1
    object$margin[[2]] <- smooth2
    object$knots <- k1
    object$m <- m
    bs.dim <- ncol(Z)
  }
  
  # __ general _____________________________________________________
  # apply Z trafo and prepare and return object 
  
  object$X <- object$X %*% Z
  object$S <- list(crossprod(Z, S) %*% Z)
  object$Z <- Z
  # object$bs.dim <- bs.dim
  object$rank <- qr(object$S[[1]])$rank
  object$null.space.dim <- bs.dim - object$rank
  # no sum-to-zero constraint for skew-symm bases:
  if(object$xt$skew) object$C <- matrix(0, 0, bs.dim)
  class(object) <- "sympspline.smooth"
  object
}


# Predict.matrix function -------------------------------------------------


# adapted from sparseFLMM:::Predict.matrix.symm.smooth

#' Predict matrix method for (skew-)symmetric smooths
#' 
#' Method generating the predict matrix for a (skew-)symmetric (cyclic) univariate
#' or bivariate P-spline smooth. Adapted from 
#' \code{\link{sparseFLMM::Predict.matrix.symm.smooth}}.
#'
#' @param object a \code{sympspline.smooth} object created by 
#' \code{\link[smooth.construct.sps.smooth.spec]{smooth.construct.sps.smooth.spec}}.
#' See \code{\link[smooth.construct]{smooth.construct}}.
#' @param data a list containing the data. 
#' See \code{\link[smooth.construct]{smooth.construct}}.
#'
#' @return A matrix which will map the parameters associated with the smooth to 
#' the vector of values of the smooth evaluated at the covariate values given in object.
#' @importFrom mgcv Predict.matrix
#' @export
Predict.matrix.sympspline.smooth <- function (object, data) {
  
  # __ 1D case _______________________________________________
  # determine designmat X and apply Z transformation matrix 
  
  # almost identical to Predict.matrix.pspline.smooth
  # only with (skew)-symmetric basis in the end
  
  if(length(object$term) == 1) {
    m <- object$m[1] + 1
    ll <- object$knots[m + 1]
    ul <- object$knots[length(object$knots) - m]
    m <- m + 1
    x <- data[[object$term]]
    n <- length(x)
    ind <- x <= ul & x >= ll
    if (is.null(object$deriv)) 
      object$deriv <- 0
    if (sum(ind) == n) {
      X <- splines::spline.des(object$knots, x, m, rep(object$deriv, 
                                                       n))$design
    }
    else {
      D <- splines::spline.des(object$knots, c(ll, ll, ul, 
                                               ul), m, c(0, 1, 0, 1))$design
      X <- matrix(0, n, ncol(D))
      nin <- sum(ind)
      if (nin > 0) 
        X[ind, ] <- splines::spline.des(object$knots, x[ind], 
                                        m, rep(object$deriv, nin))$design
      if (object$deriv < 2) {
        ind <- x < ll
        if (sum(ind) > 0) 
          X[ind, ] <- if (object$deriv == 0) 
            cbind(1, x[ind] - ll) %*% D[1:2, ]
        else matrix(D[2, ], sum(ind), ncol(D), byrow = TRUE)
        ind <- x > ul
        if (sum(ind) > 0) 
          X[ind, ] <- if (object$deriv == 0) 
            cbind(1, x[ind] - ul) %*% D[3:4, ]
        else matrix(D[4, ], sum(ind), ncol(D), byrow = TRUE)
      }
    }
    # apply (skew-)symmetry constraint
    if(object$xt$skew) {
      bs.dim <- floor(object$bs.dim/2)
      X <- X[, 1:bs.dim] - 
        X[, ncol(X)+1 - (1:bs.dim)]
    } else {
      bs.dim <- ceiling(object$bs.dim/2)
      X <- X[, 1:bs.dim] + 
        X[, ncol(X)+1 - (ifelse(ncol(X)%%2, 2, 1):bs.dim)]
    }
    
    if (object$mono == 0) 
      return(X)
    else return(X %*% object$Bs)
  }
  
  # __ 2D case _______________________________________________
  # determine designmat X and apply Z transformation matrix 
  
  # almost identical to Predict.matrix.symm.smooth
  # only also allowing for the skew-symmetric option
  # in make_summation_matrix
  
  if(length(object$term) == 2) {
    m <- length(object$margin)
    X <- list()
    for (i in 1:m) {
      term <- object$margin[[i]]$term
      dat <- list()
      for (j in 1:length(term)) {
        dat[[term[j]]] <- data[[term[j]]]
      }
      X[[i]] <- PredictMat(object$margin[[i]], dat, n = length(dat[[1]]))
    }
    X <- tensor.prod.model.matrix(X)
    if(is.null(object$Z)) {
      Z <- make_summation_matrix(F = object$bs.dim, skew = object$xt$skew, 
                                 cyclic.degree = object$xt$cyclic * 
                                   (object$m[1]+1) )
    } else {
      Z <- object$Z
    }
      
    X %*% Z
  }
  
}

