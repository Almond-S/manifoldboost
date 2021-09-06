# Define orthogonalization operator on mboost base-learners, orthogonalizing the first with respect to the second
# Code is adapted from FDboost::Xc

#' Base Learner Basis Orthogonalization Operator
#' 
#' Operator orthogonalizing the design matrix of one `mboost` baselearner `bl1` 
#' with respect to the design matrix of another baselearner `bl2`, i.e., 
#' the design matrix of `bl1` is linearly transformed, such that all columns in
#' it are orthogonal to the columns of the design matrix of `bl2`.  
#'
#' @param bl1 a baselearner generator of class `blg` 
#' @param bl2 another baselearner generator of class `blg` which is typically 
#' contained in `bl1` in the sense that the columns of its design matrix are
#' in the span of the design matrix of `bl1`. 
#'
#' @return object of class `blg` with the transformed baselearner `bl1`.
#' @export
#' @import Matrix
#'
'%-%' <- function (bl1, bl2) 
{
  if (is.list(bl1) && !inherits(bl1, "blg")) 
    return(lapply(bl1, "%-%", bl2 = bl2))
  if (is.list(bl2) && !inherits(bl2, "blg")) 
    return(lapply(bl2, "%-%", bl1 = bl1))
  cll <- paste(bl1$get_call(), "%-%", bl2$get_call(), collapse = "")
  stopifnot(inherits(bl1, "blg"))
  stopifnot(inherits(bl2, "blg"))
  # used_bl <- c(deparse(match.call()$bl1[[1]]), deparse(match.call()$bl2[[1]]))
  
  # if (any(!used_bl %in% c("bols", "brandom", "bbs"))) {
    # warning("%-% is intended to operate on base-learners bols, brandom and bbs.")
  # }
  # stopifnot(!any(colnames(mboost_intern(bl1, fun = "model.frame.blg")) %in% 
  #                  colnames(mboost_intern(bl2, fun = "model.frame.blg"))))
  mf <- cbind(mboost_intern(bl1, fun = "model.frame.blg"), 
              mboost_intern(bl2, fun = "model.frame.blg"))
  index1 <- bl1$get_index()
  index2 <- bl2$get_index()
  if (is.null(index1)) 
    index1 <- 1:nrow(mf)
  if (is.null(index2)) 
    index2 <- 1:nrow(mf)
  mfindex <- cbind(index1, index2)
  index <- NULL
  CC <- all(mboost_intern(mf, fun = "Complete.cases"))
  if (!CC) 
    warning("base-learner contains missing values;\n", "missing values are excluded per base-learner, ", 
            "i.e., base-learners may depend on different", " numbers of observations.")
  DOINDEX <- (nrow(mf) > options("mboost_indexmin")[[1]])
  if (is.null(index)) {
    if (!CC || DOINDEX) {
      index <- mboost_intern(mfindex, fun = "get_index")
      mf <- mf[index[[1]], , drop = FALSE]
      index <- index[[2]]
    }
  }
  vary <- ""
  ret <- list(model.frame = function() if (is.null(index)) return(mf) else return(mf[index, 
                                                                                     , drop = FALSE]), get_call = function() {
                                                                                       cll <- deparse(cll, width.cutoff = 500L)
                                                                                       if (length(cll) > 1) cll <- paste(cll, collapse = "")
                                                                                       cll
                                                                                     }, get_data = function() mf, get_index = function() index, 
              get_vary = function() vary, get_names = function() colnames(mf), 
              set_names = function(value) attr(mf, "names") <<- value)
  class(ret) <- "blg"
  args <- environment(bl1$dpp)$args[c("lambda", "df")]
  l1 <- args$lambda
  Xfun <- function(mf, vary, args) {
    if (is.null(args$prediction)) 
      args$prediction <- FALSE
    newX1 <- environment(bl1$dpp)$newX
    newX2 <- environment(bl2$dpp)$newX
    X1 <- newX1(mf[, bl1$get_names(), drop = FALSE], prediction = args$prediction)
    K1 <- X1$K
    X1 <- X1$X
    if (!is.null(l1)) 
      K1 <- l1 * K1
    MATRIX <- options("mboost_useMatrix")$mboost_useMatrix
    if (MATRIX & !is(X1, "Matrix")) 
      X1 <- Matrix(X1)
    if (MATRIX & !is(K1, "Matrix")) 
      K1 <- Matrix(K1)
    X2 <- newX2(mf[, bl2$get_names(), drop = FALSE], prediction = args$prediction)
    X2 <- X2$X
    if (MATRIX & !is(X2, "Matrix")) 
      X2 <- Matrix(X2)
    if (is.null(args$Z)) {
      C <- t(X1) %*% X2
      qr_C <- qr(C)
      if (any(class(qr_C) == "sparseQR")) {
        rank_C <- qr_C@Dim[2]
      }
      else {
        rank_C <- qr_C$rank
      }
      Q <- qr.Q(qr_C, complete = TRUE)
      args$Z <- Q[, (rank_C + 1):ncol(Q)]
    }
    X <- X1 %*% args$Z
    K <- t(args$Z) %*% K1 %*% args$Z
    list(X = X, K = K, args = args)
  }
  temp <- Xfun(mf = mf, vary = vary, args = args)
  args$Z <- temp$args$Z
  rm(temp)
  ret$dpp <- mboost_intern(ret, Xfun = Xfun, args = args, fun = "bl_lin")
  return(ret)
}


# function directly applying basis transformation by matrix multip --------

#' Baselearner Basis Transformation Operator
#' 
#' A function transforming the basis of a baselearner linearly with a 
#' basis transformation matrix \eqn{Z}, i.e., if \eqn{B} is the design matrix 
#' of the baselearner, the design matrix of the transformed learner is 
#' \eqn{\tilde{B} = BZ}.
#'
#' @param bl1 a baselearner generator object of class `blg`.
#' @param Z a basis transformation matrix with its number of rows corresponding
#' to the number of basis functions of `bl1`.
#'
#' @return an object of class `blg` with the transformed baselearner.
#' @importFrom Matrix crossprod
#'
btrafo <- function (bl1, Z) 
{
  if (is.list(bl1) && !inherits(bl1, "blg")) 
    return(lapply(bl1, "%T%", Z = Z))
  cll <- paste(bl1$get_call(), "%T%", match.call()$Z, collapse = "")
  stopifnot(inherits(bl1, "blg"))
  stopifnot(inherits(Z, "matrix"))
  
  mf <- mboost_intern(bl1, fun = "model.frame.blg")
  index1 <- bl1$get_index()
  if (is.null(index1)) 
    index1 <- 1:nrow(mf)
  mfindex <- index1
  index <- NULL
  CC <- all(mboost_intern(mf, fun = "Complete.cases"))
  if (!CC) 
    warning("base-learner contains missing values;\n", "missing values are excluded per base-learner, ", 
            "i.e., base-learners may depend on different", " numbers of observations.")
  DOINDEX <- (nrow(mf) > options("mboost_indexmin")[[1]])
  if (is.null(index)) {
    if (!CC || DOINDEX) {
      index <- mboost_intern(mfindex, fun = "get_index")
      mf <- mf[index[[1]], , drop = FALSE]
      index <- index[[2]]
    }
  }
  vary <- ""
  ret <- list(model.frame = function() if (is.null(index)) return(mf) else return(mf[index, 
                                                                                     , drop = FALSE]), get_call = function() {
                                                                                       cll <- deparse(cll, width.cutoff = 500L)
                                                                                       if (length(cll) > 1) cll <- paste(cll, collapse = "")
                                                                                       cll
                                                                                     }, get_data = function() mf, get_index = function() index, 
              get_vary = function() vary, get_names = function() colnames(mf), 
              set_names = function(value) attr(mf, "names") <<- value)
  class(ret) <- "blg"
  args <- environment(bl1$dpp)$args
  l1 <- args$lambda
  # if (!is.null(l1) && !is.null(l2)) {
  #   args <- list(lambda = 1, df = NULL)
  # }
  # else {
  #   args <- list(lambda = NULL, df = ifelse(is.null(args$df), 
  #                                           1, args$df) * ifelse(is.null(args2$df), 1, args2$df))
  # }
  Xfun <- function(mf, vary, args) {
    if (is.null(args$prediction)) 
      args$prediction <- FALSE
    newX1 <- environment(bl1$dpp)$newX
    X1 <- newX1(mf[, bl1$get_names(), drop = FALSE], prediction = args$prediction)
    K1 <- X1$K
    X1 <- X1$X
    if (!is.null(l1)) 
      K1 <- l1 * K1
    MATRIX <- options("mboost_useMatrix")$mboost_useMatrix
    if (MATRIX & !is(X1, "Matrix")) 
      X1 <- Matrix(X1)
    if (MATRIX & !is(K1, "Matrix")) 
      K1 <- Matrix(K1)
    
    ### apply basis trafo
    if (is.null(args$Z)) {
      if(nrow(Z) != ncol(X1)) 
        stop(paste("The transformation matrix Z has to have", ncol(X1), "rows like the effect design matrix."))
      if(ncol(Z) > ncol(X1)) 
        stop(paste("The transformation matrix Z has to have at most", ncol(X1), "cols, 
                  the number of columns of the effect design matrix."))
      args$Z <- Z
    }
    X <- X1 %*% args$Z
    K <- crossprod(args$Z, K1) %*% args$Z
    list(X = X, K = K, args = args)
  }
  temp <- Xfun(mf = mf, vary = vary, args = args)
  args$Z <- temp$args$Z
  rm(temp)
  ret$dpp <- mboost_intern(ret, Xfun = Xfun, args = args, fun = "bl_lin")
  return(ret)
}
