#' Complex planar curve interpreter
#'
#' @param formula a formula in the form `shape ~ value^dim ~ arg | id`.
#' @param data data containing the shape in a format convertible to `shape_dat`. 
#' Shape should be contained as `tbl_cube` (`shape_cube`) if data is regular.
#' @param weights 
#' @param mf a mfGeometry object
#'
#' @return a list containing 
#' - mod0_data: a data.frame with the dim, id, and arg values needed
#' to set up a basis for the response object 
#' - weights: a vector/matrix (regular case) or a list of vectors/matrices of 
#' (integration) weights
#' - response: an array (regular case) or vector (irregular case) containing the 
#' response values
#'
# #' @examples
mfInterpret_cplx_curves <- function(
  formula,
  data,
  weights = NULL,
  arg.range = NULL
) {
  
  ### check formulas
  stopifnot(class(formula) == "formula")
  
  # remove everything irrelevant form data
  data <- data[which(names(data) %in% all.vars(formula))]
  
  # convert to shape_dat
  data <- as_shape_dat(data, formula)
  # 
  ### prepare data
  # -> separate values (response) from dim, arg and id (data)
  #  -> in the regular case, response: array
  #  -> in the irregular case, response: long
  
  # format correctly and get response(-values)
  if(inherits(data, "shape_cube")) {
    # variables in formula
    av <- all.vars(formula)
    names(av) <- c("shape", "value", "dim", "arg", "id")
    shp <- av["shape"]
    yname <- av["value"]
    # arrange and store response array
    ydims <- names(data[[shp]]$dims)
    stopifnot(all(ydims %in% av[c("arg", "dim", "id")]))
    ydimorder <- ordered(av[c("arg", "dim", "id")], levels = ydims)
    ydimorder <- as.numeric(ydimorder)
    response <- data[[shp]]$mets[[yname]]
    response <- aperm(response, ydimorder)
    mod0_data <- data[[shp]]$dims[ydimorder]
    # response in FDboost format
    mod0_data[[av["value"]]] <- t(matrix(response, ncol = dim(response)[3]))
  } else {
    # convert to shape_dat_frame_long and store response
    mod0_data <- as_shape_frame_long(data)
    # variables in formula
    av <- all.vars(formula(mod0_data))
    names(av) <- c("value", "dim", "arg", "id")
    shp <- shape_names(mod0_data)[["shape"]]
    yname <- av["value"]
    response <- mod0_data[[yname]]
  }
  
  ret <- list(
    variable_names = av,
    mod0_data = mod0_data,
    response = response
  )
  
  ### initialize base-learners
  if(is.array(response)) {
    if(length(unique(mod0_data[[av["dim"]]]))!=2) 
      stop("Sorry, I can only interpret planar 2-dimensional shapes.")
  } else 
    mod0_data <- as_shape_frame_complex(mod0_data)
  
  ### extract id
  nameid <- av["id"]
  .id <- mod0_data[[nameid]]
  
  ### check id/response
  if(!is.array(response)) {
    ## basic checks
    if(!is.null(dim(response))) 
      stop("If not supplied as an array, the response has to be a vector.")
    if(length(.id) != length(response)/2)
      stop("id must be of same length as response.")
  }
  
  ### prepare arg
  nameyarg <- av["arg"]
  .arg <- mod0_data[[nameyarg]]
  if(!is.array(response)) {
    stopifnot(length(.arg) == length(response)/2)
    id_ <- factor(.id, levels = unique(.id))
  }
  
  ### make numerical integration weights for arg
  if(is.character(weights)) {
    ### function setting up integration weights
    make_integration_weights <- function(weights, arg = NULL, range = NULL) {
      switch (weights,
              equal = NULL,
              trapezoidal = trapez_weights(arg, range = range),
              simpson = simpson_weights(arg, range = range)
      )
    } 
    if(is.array(response)) {
      A <- make_integration_weights(weights, .arg, arg.range)
    } else { 
        ## combine relevant info in one matrix
        arg_ <- split(.arg, id_)
        A <- lapply(arg_, make_integration_weights, 
                                   weights = weights, range = arg.range)
    } 
    
    # store integration weights
    ret$weights <- A
  } else {
    if(!is.null(weights)) stopifnot(length(weights) == length(.arg))
  }
  
  ret
} # mfInterpret_cplx_curves





#' Full Procrustres mean of planar shapes 
#'
#' @param formula a formula in the form `shape ~ value^dim ~ blearner(arg) | id`,
#' where `blearner` is typically `bbs` if `arg` is continuous and `bols` if discrete.
#' @param data data containing the shape in a format convertible to `shape_dat`. 
#' Shape should be contained as `tbl_cube` (`shape_cube`) if data is regular.
#' @param cyclic 
#' @param smoothed.cov 
#' @param cov.k 
#' @param arg.grid.len 
#' @param weights 
#' @param numInt
#' @param mf
#' @param mfboost.return 
#'
#' @return in the default setting, the estimated full procrustres mean in 
#' shape_frame_long format 
#' (data.frame with coordinates as real values stored in a single column).
#' @import mboost
#' @export
#'
# #' @examples
planarshape_full_proc <- function(
  formula,
  data,
  cyclic = FALSE,
  smoothed.cov = FALSE,
  cov.k = 10,
  arg.grid.len = 50,
  weights = NULL,
  arg.range = NULL,
  mf = NULL,
  mfboost.return = FALSE
) {
  
  basics <- mfInterpret_cplx_curves(formula = formula, 
                                    data = data, 
                                    weights = weights,
                                    arg.range = arg.range)
  mod0_data <- basics$mod0_data
  response <- basics$response
  av <- basics$variable_names
  
  if(is.array(response))
    mod0_data[[av["value"]]] <- response
  
  ### initialize base-learners
  
  if(is.array(response)) {
    arg_form <- if(formula[[3]][[1]] == "|") 
      formula[[3]][[2]] else 
        formula[[3]]
    val_form <- formula[[2]][[3]][[2]]
  } else {
    arg_form <- if(formula[[3]][[1]] == "|") 
      formula[[3]][[2]] else 
        formula[[3]]
    val_form <- formula(mod0_data)[[2]][[2]]
  }
  mod0_formula <- as.formula(
    paste0(deparse(val_form), "~", 
           deparse(arg_form)),
    env = environment(formula))
  
  g2 <- Gaussian()
  g2@check_y <- function(y) y # just to avoid any errors
  g2@offset <- function(y, weights) 0
  mod0 <- mboost(mod0_formula, mod0_data, 
                 family = g2, control = boost_control(mstop = 0))
  if(length(mod0$baselearner) != 1) 
    stop("Exactly one base-learner has to be specified for the pole.")
  bl <- mod0$baselearner[[1]]
  
  ### extract id
  nameid <- av["id"]
  .id <- mod0_data[[nameid]]
  
  ### extract dim
  ydimvar <- av["dim"]
  nameydim <- ydimvar
  .dim <- mod0_data[[ydimvar]]
  if(length(unique(.dim))!=2) 
    stop("Only planar 2-dimensional shapes provided, yet.")
  
  ### check id/response
  if(!is.array(response)) {
    ## basic checks
    if(!is.null(dim(response))) 
      stop("If not supplied as an array, the response has to be a vector.")
    if(length(.id) != length(response))
      stop("id must be of same length as response.")
  }
  
  ### prepare arg
  nameyarg <- av["arg"]
  .arg <- mod0_data[[nameyarg]]
  if(!is.array(response)) stopifnot(length(.arg) == length(response))
  if(is.null(arg.range)) arg.range <- environment(bl$dpp)$args$knots[[1]]$boundary.knots
  # make arg.grid for numerical integrations
  arg.grid <- if(is.array(response)) .arg else {
    if(is.numeric(.arg)) 
      head(seq(arg.range[1], arg.range[2], len = arg.grid.len), arg.grid.len - cyclic)  else
        stop("Irregular (long data format) version for non-numeric arg-values 
           not implemented, yet.") #unique(.arg)
  }
  
  ### make numerical integration weights for arg.grid
  ## TODO: improve implementation of weights using mfInterpret and
  # allowing for custom weights !!!
  if(is.null(weights)) weights <- "equal"
  ### function setting up integration weights
  make_integration_weights <- function(weights, arg = NULL, range = NULL) {
    switch (weights,
            equal = NULL,
            trapezoidal = trapez_weights(arg, range = range),
            simpson = simpson_weights(arg, range = range)
    )
  } 
  A <- make_integration_weights(weights, arg.grid, arg.range)
  
  B <- as.matrix(extract(bl, "design"))
  # Design matrix evaluated on grid
  B.grid <- if(is.array(response)) 
    B else {
      newdat <- data.frame(.arg = arg.grid)
      names(newdat) <- nameyarg
      # evaluate B on integration grid
      B.grid <- as.matrix(with(newdat, 
                               mboost::extract(eval(str2expression(bl$get_call())), "design")))
    }
  
  
  if(is.null(smoothed.cov)) {
    smoothed.cov <- !is.array(response) & is.numeric(.arg)
  }
  
  ### get estimate of ByyB, the complex covariance matrix of the basis coefs
  ## function returning ByyB
  compute_ByyB <- function(value, B, int.weights = NULL, yy) {
    if(missing(yy)) yy <- tcrossprod(value, Conj(value))
    ByyB <- if(is.null(int.weights))
      crossprod(B, yy) %*% B else {
        if(is.null(dim(int.weights)))
          crossprod(int.weights*B, yy) %*% (int.weights*B) else
            crossprod(int.weights %*% B, yy) %*% int.weights %*% B
      }
    ByyB
  }
  
  ### ... without covariance smoothing (regular / irregular case) or with cov smooth
  if(!smoothed.cov) {
    # could be extended to regularized FPCA
    # a la Silverman (1996): smoothed functional principal component analysis 
    # by choice of norm
    # but didn't work so well and is, thus, not implemented here
    # => if penalization necessary always use smoothed.cov
    
    
    if(is.array(response)) {
      as_complex <- function(x) {
        matrix(complex(real = x[, 1, ], imaginary = x[, 2, ]), nrow = dim(x)[1])
      }
      ByyB <- compute_ByyB(as_complex(response), B, A)
    } else {
      as_complex <- function(x) {
        complex(real = x[, 1], imaginary = x[, 2])
      }
      ## combine relevant info in one matrix
      B_ <- as.data.frame(cbind(response, .arg, B))
      B_ <- split(B_, .id)
      B_ <- lapply(B_, as.matrix)
      ByyB <- sapply(B_, function(B_id) {
        A <- make_integration_weights(weights, 
                                      arg = B_id[1:(nrow(B_id)/2), 2], arg.range)
        compute_ByyB(
          value = as_complex(matrix(B_id[,1], ncol = 2)),
          B = B_id[1:(nrow(B_id)/2), -(1:2)], A)
      })
      ByyB <- matrix( rowSums(ByyB), ncol = ncol(B) )
    }
  } else {
    # -> smoothed.cov == TRUE
    # ==> no distinction between regular and irregular design needed
    # use a PACE-type procedure generalizing the 
    # Cederbaum, Scheipl, Greven (2016): Fast symmetric additive covariance smoothing
    # to complex covariance matrices
    require(mgcv)
    
    if(!is.numeric(.arg)) 
      stop("If 'smoothed.cov = TRUE', the arg variable in the formula has to be numeric.")
    
    # we model the part of E[Conj(y(t1))y(t2)] where t1 < t2 
    # assuming that the time points are ordered
    cresponse <- as_shape_frame_complex(mod0_data)
    
    cov_dat <- lapply(split(cresponse, cresponse[[av["id"]]]), function(x) {
      combs <- combn(1:nrow(x), 2)
      data.frame(
        yy = x[[av["value"]]][combs[1,]] * Conj(x[[av["value"]]][combs[2,]]), 
        arg1 = x[[av["arg"]]][combs[1,]],
        arg2 = x[[av["arg"]]][combs[2,]]
      )
    })
    cov_dat <- do.call(rbind, cov_dat)
    
    cov_fit_re <- bam( Re(yy) ~ s(arg1, arg2, bs = "sps", k = cov.k, 
                                  xt = list(cyclic = cyclic)), 
                       data = cov_dat, knots = list(arg1 = arg.range) )
    cov_fit_im <- bam( Im(yy) ~ -1 + s(arg1, arg2, bs = "sps", k = cov.k,  
                                       xt = list(skew = TRUE, cyclic = cyclic)), 
                       data = cov_dat, knots = list(arg1 = arg.range ) )
    
    # predict smoothed covariance
    cov_dat <- expand.grid(arg1 = arg.grid, arg2 = arg.grid)
    yy <- matrix( complex(
      real = predict(cov_fit_re, newdata = cov_dat),
      imaginary = predict(cov_fit_im, newdata = cov_dat) ), 
      ncol = length(arg.grid))
    
    ByyB <- compute_ByyB(yy = yy, B = B.grid, int.weights = A)
    
  }
  
  ### compute matrix R = BB
  if(is.null(A)) 
    R <- crossprod(B.grid) else {
      if(is.null(dim(A))) 
        R <- crossprod(B.grid, A * B.grid) else {
          R <- crossprod(B.grid, A %*% B.grid)
        }
    }
  R <- as.matrix(R)
  
  ### Compute leading eigenvector
  ei <- eigen(as.matrix(solve(R)) %*% ByyB)
  pole_coefs = ei$vectors[,1]
  
  #### Prepare output
  
  if(is.array(response)) {
    pole_ <- B %*% pole_coefs
  } else {
    dim1 <- .dim == unique(.dim)[1]
    pole_ <- B[dim1, ] %*% pole_coefs
  }
  
  if(mfboost.return) {
    ### return pole_ (together with pole_ coefficients and tangent space)
    tangent_basis <- complex2realmat(Null(R %*% cbind(pole_coefs, 1)))
    tangent <- function( bl ) btrafo(bl, tangent_basis )
    # pole_X <- complex2realmat(pole_)
    # tangent <- function( bl ) bl %-% buser(pole_X)
    
    pole_ <- if(is.array(response)) 
      lapply(seq_len(dim(response)[3]), function(i) pole_) else {
        if(!exists("cresponse")) 
          cresponse <- as_shape_frame_complex(data)
        # keep id arrangement  
        ids <- cresponse[[nameid]]
        ids <- ordered(ids, levels = unique(ids))
        split(pole_, ids)
      }
    attr(pole_, "tangent") <- tangent
  } else {
    ### return pole_ (together with pole_ coefficients and tangent space)
    pole_ <- if(is.array(response)) cbind(Re(pole_), Im(pole_)) else {
      if(!exists("cresponse")) cresponse <- as_shape_frame_complex(mod0_data)
      cresponse[[av["value"]]] <- pole_
      as_shape_frame_long(cresponse)
    }
  }
  
  pole_
} # planarshape_full_proc