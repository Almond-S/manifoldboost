

o_control <- function(smoothed.cov = NULL, 
                      smoothed.cov_control = NULL, 
                      weights = NULL) {
  list(smoothed.cov = smoothed.cov, smoothed.cov_control = smoothed.cov_control,
       weights = weights)
}

mfboost <- function(formula,               # response ~ xvars 
                       argformula,           # formula for response argument
                       dim = NULL, # dim variable if resp is in long format
                       id = NULL,  # id variable if response is in long format
                       numInt = "equal",      # numerical integration method
                       data,                  # list of response (array), etc.
                       offset = NULL,         # so far only full procrustres
                       offset_control = o_control(), # further specifications. 
                       #to come
                       ...) {                 # further arguments passed to mboost
  ### 
  # Great parts of this function (basic logic, formula and argument handling, etc.) 
  # are directly adapted from the R function FDboost::FDboost.
  # Many thanks to David Ruegamer and in particular to Sarah Brockhaus!!!
  ###
  
  dots <- list(...)
  
  ret <- list() # list items returned in the end by the function
  
  ### save formula of shapeboost before it is changed
  formula_shapeboost <- formula
  
  tf <- terms.formula(formula, specials = c("c"))
  trmstrings <- attr(tf, "term.labels")
  equalBrackets <- NULL
  if(length(trmstrings) > 0){
    ## insert id at end of each base-learner
    trmstrings2 <- paste(substr(trmstrings, 1 , nchar(trmstrings)-1), ", index=", id[2],")", sep = "")
    ## check if number of opening brackets is equal to number of closing brackets
    equalBrackets <- sapply(1:length(trmstrings2), function(i)
    {
      sapply(regmatches(trmstrings2[i], gregexpr("\\(", trmstrings2[i])), length) ==
        sapply(regmatches(trmstrings2[i], gregexpr("\\)", trmstrings2[i])), length)
    })
  }
  
  ### check formulas
  if(missing(argformula) || class(try(argformula)) == "try-error") 
    stop("argformula must either be NULL or a formula object.")
  stopifnot(class(formula) == "formula")
  if(!is.null(argformula)) stopifnot(class(argformula) == "formula")
  if(class(try(id)) == "try-error") stop("id must either be NULL or a formula object.")
  if(class(try(dim)) == "try-error") stop("dim must either be NULL or a formula object.")
  
  ### split object
  
  
  
  
  ### insert the id variable into the formula, to treat it like the other variables
  if(!is.null(id)){
    stopifnot(class(id) == "formula")
    ##tf <- terms.formula(formula, specials = c("c"))
    ##trmstrings <- attr(tf, "term.labels")
    ##equalBrackets <- NULL
    if(length(trmstrings) > 0){
      ## insert index into the other base-learners of the tensor-product as well
      for(i in 1:length(trmstrings)){
        if(grepl( "%X", trmstrings2[i])){
          temp <- unlist(strsplit(trmstrings2[i], "%X"))
          temp1 <- temp[-length(temp)]
          ## http://stackoverflow.com/questions/2261079
          ## delete all trailing whitespace
          trim.trailing <- function (x) sub("\\s+$", "", x) 
          temp1 <- trim.trailing(temp1)
          temp1 <- paste(substr(temp1, 1 , nchar(temp1)-1), ", index=", id[2],")", sep = "")
          trmstrings2[i] <- paste0(paste0(temp1, collapse = " %X"), " %X", temp[length(temp)]) 
        } 
        ## do not add index to base-learners bhistx()
        if( grepl("bhistx", trmstrings[i]) ) trmstrings2[i] <- trmstrings[i] 
        ##  do not add an index if an index is already part of the formula
        if( grepl("index[[:blank:]]*=", trmstrings[i]) ) trmstrings2[i] <- trmstrings[i]
        ##  do not add an index if an index for %A%, %A0%, %O%
        if( grepl("%A%", trmstrings[i]) ) trmstrings2[i] <- trmstrings[i]
        if( grepl("%A0%", trmstrings[i]) ) trmstrings2[i] <- trmstrings[i]
        if( grepl("%O%", trmstrings[i]) ) trmstrings2[i] <- trmstrings[i]
        ##  do not add an index for base-learner that do not have brackets
        if( i %in% which(!equalBrackets) ) trmstrings2[i] <- trmstrings[i]
      }
      trmstrings <- trmstrings2
    } 
    xpart <- paste(as.vector(trmstrings), collapse = " + ")
    if(xpart != ""){
      if(any(substr(tf[[3]], 1, 1) == "1")) xpart <- paste0("1 + ", xpart)
    }else{
      xpart <- 1
    }
    formula <- as.formula(paste(tf[[2]], " ~ ", xpart))
    #print(formula)
    nameid <- paste(id[2])
    .id <- data[[nameid]]
  }else{
    nameid <- NULL
    .id <- NULL
  }
  
  ### extract time (argument) from argformula 
  yarg <- all.vars(argformula)[[1]]
  stopifnot(length(yarg) == 1)
  nameyarg <- yarg
  assign(yarg, data[[yarg]])
  .arg <- data[[yarg]]
  if(!is.numeric(.arg) & !is.factor(.arg)) 
    message("Non-numeric arg values are interpreted as factor
                               and lead to landmark shape models.")
  attr(.arg, "nameyarg") <- nameyarg
  
  ### extract response; a numeric array or vector
  yname <- all.vars(formula)[1]
  response <- data[[yname]]
  if(is.null(response)) stop("The response <", yname, 
                             "> is not contained in data.")
  ydim <- dim(response)
  
  if(is.array(response)){
    ### check dimensions
    if(length(ydim) != 3) 
      stop("If supplied as array, response has to be three dimensional.")
    if(ydim[2] != 2)
      stop("Only planar response curves provided, yet.")
    
    stopifnot(ydim[1] == length(.arg))
    
    ### extract dim if supplied
    if(!is.null(dim)) {
      stopifnot(class(dim) == "formula")
      ydimvar <- all.vars(dim)[[1]]
      nameydim <- ydimvar
      stopifnot(length(ydimvar) == 1)
      .dim <- data[[ydimvar]]
      if(length(unique(.dim))!=2) 
        stop("Only planar 2-dimensional shapes provided, yet.")
    } else {
      nameydim <- NULL
      if(is.array(response)) {
        .dim <- array( rep(c(1,2), each = length(response)/2), dim = ydim[c(1,3,2)])
        .dim <- c(aperm(.dim, c(1,3,2)))
      } else {
        .dim <- sapply(table(.id), function(x) rep(c(1,2), each = x/2))
      }
    }
    .dim <- factor(.dim)
    
    # column-wise stacking of response 
    dresponse <- as_shape_long(response, 
                               arg = .arg, 
                               dim = NULL, 
                               id = .id)$.value
    nobs <- ydim[3] # number of observed trajectories
    
    
    if(is.null(offset_control$smoothed.cov)) 
      offset_control$smoothed.cov = FALSE
    offset_control$arg.grid = unique(.arg)
  }else{
    stopifnot(is.null(dim(response))) 
    # check length of response and its argument and index
    stopifnot(length(response) == length(.arg) & length(response) == length(.id))
    
    if(any(is.na(response))) 
      warning("For non-grid observations the response should not contain missing values.")
    # if( !all(sort(unique(id)) == 1:length(unique(id))) ) stop("id has to be integers 1, 2, 3,..., N.")
    
    # response already in long format  
    dresponse <- response
    nobs <- length(unique(.id)) # number of observed trajectories
    
    if(is.null(offset_control$smoothed.cov)) 
      offset_control$smoothed.cov = TRUE
    offset_control$arg.grid = NULL
  }
  
  ## convert complex linear operator matrix to real matrix
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
  vec2op <- function(x) c(-tail(x, length(x)/2),
                                        head(x, length(x)/2))
  vec2cplx <- function(x) complex(
    real = tail(x, length(x)/2), im = head(x, length(x)/2))
  
  
  ##### handle offset
  if(!is.null(offset)) {
    if(is.array(response) & length(offset) != length(.arg)*2)
      stop("For regular response given as kx2xn array, 
           a user-specified offset has to be of length k*2.")
    if(is.array(response) & length(offset) != length(.arg))
      stop("For irregular response given in long format,
           a user-specified offset has to have the same length as the response.")
  }
  
  offsetformula <- update(argformula, as.formula(paste(yname, "~.")))
  if(is.null(offset)) {
    ### create offset from response and argformula
    # make offset formula
    shp_offset <- shapeboost_offset(offsetformula, dim, id, data, 
                                    cyclic, 
                                    smoothed.cov = offset_control$smoothed_cov,
                                    smoothed.cov_control = offset_control$smoothed.cov_control,
                                    weights = offset_control$weights)
    Toffset <- shp_offset$offset_tangent_coefs
    shp_offset <- shp_offset$offset
    
    Toffset <- cmat2vmat(Toffset)
    if(is.array(response)) offset <- c(Re(shp_offset), Im(shp_offset)) else
      offset <- tapply(shp_offset, id[table(id)/2], function(x)
        c(Re(x), Im(x)))
  } else { # for a user-specified offset
    
    ### initialize base-learner
    g2 <- Gaussian()
    g2@check_y <- function(y) y # just to avoid any errors
    g2@offset <- function(y, weights) 0
    mod0 <- mboost(offsetformula, data, 
                   family = g2, control = boost_control(mstop = 0))
    if(length(mod0$baselearner) != 1) 
      stop("Exactly one base-learner has to be specified for offset.")
    bl <- mod0$baselearner[[1]]
    
    ### find orthogonalization matrix with respect to offset
    X <- extract(bl, "design")
    shp_offset <- if(is.array(response))
      vec2cplx(offset) else 
        unlist(tapply(offset, .id, vec2cplx, simplify = FALSE))
    
    C <- crossprod(as.matrix(X), shp_offset)
    qr_C <- qr(C)
      rank_C <- qr_C$rank
    Q <- qr.Q(qr_C, complete = TRUE)
    Toffset <- Q[, (rank_C + 1):ncol(Q)]
    # convert to long format
    Toffset <- cmat2vmat(Toffset)
  }
  
  ### delete response variables from data
  data[[yname]] <- NULL
  data[[yarg]] <- NULL
  if(!is.array(response)) {
    data[[ydimvar]] <- NULL
    data[[nameid]] <- NULL
  }
  
  ### extract covariates
  allCovs <- unique(c(nameid, all.vars(formula)))
  if(length(allCovs) > 1){
    data <- data[allCovs[!allCovs %in% c(yname, nameyarg, nameydim)] ]
    if( any(is.na(names(data))) ) data <- data[ !is.na(names(data)) ]
  }else{
    data <- list(NULL)  # <SB> intercept-model without covariates
  } 
  
  ##### compose mboost formula
  
  ## get formula over time
  tfm <- paste(deparse(argformula), collapse = "")
  tfm <- strsplit(tfm, "~")[[1]]
  tfm <- strsplit(tfm[2], "\\+")[[1]]
  
  ## transform argformula to 2D and into tangent space of offset
  tfm <- paste("brandom(.dim, lambda = 0) %Xa0%", tfm, 
               "%T% Toffset")
  
  ## get formula in covariates 
  cfm <- paste(deparse(formula), collapse = "") 
  cfm <- strsplit(cfm, "~")[[1]]
  cfm0 <- cfm
  #xfm <- strsplit(cfm[2], "\\+")[[1]]
  xfm <- trmstrings
  
  # expand formula as Kronecker or tensor product 
  if(is.array(response)){
    tmp <- outer(xfm, tfm, function(x, y) paste(y, x, sep = " %O% "))
  }else{
    ## expand the bl according to id
    # do not expand for terms without brackets, which is equal to having an unequal number of brackets
    # in the generation part of trmstrings
    if(is.null(equalBrackets)){ # for intercept models y ~ 1
      which_equalBrackets <- 0
    } else{
      which_equalBrackets <- which(equalBrackets)
    }
    xfmTemp <- paste(substr(xfm[which_equalBrackets], 1 , 
                            nchar(xfm[which_equalBrackets]) - 1 ), ")", sep = "") # , index=id is done in the beginning
    xfm[which_equalBrackets] <- xfmTemp
    rm(xfmTemp)
    tmp <- outer(xfm, tfm, function(x, y) paste(y, x, sep = "%X%"))
  }
  
  ####### find the number of df for each base-learner 
  ## for a fair selection of bl the df must be equal in all bl
  get_df <- function(bl){
    split_bl <- unlist(strsplit(bl, split = "%.{1,3}%"))
    all_df <- c()
    for(i in 1:length(split_bl)){
      parti <- parse(text = split_bl[i])[[1]] 
      parti <- FDboost:::expand.call(definition = get(as.character(parti[[1]])), call = parti)
      dfi <- parti$df # df of part i in bl 
      if(is.symbol(dfi) || (!is.numeric(dfi) && is.numeric(eval(dfi)))) dfi <- eval(dfi) 
      lambdai <- parti$lambda # if lambda is present, df is ignored 
      if(is.symbol(lambdai)) lambdai <- eval(lambdai)
      if(!is.null(dfi)){
        all_df[i] <- dfi 
      }else{ ## for df = NULL, the value of lambda is used 
        if(lambdai == 0){
          all_df[i] <-  NCOL(extract(with(data, eval(parti)), "design"))
        }else{
          all_df[i] <- "" ## dont know df 
        }
        if(grepl("%X.{0,3}%", bl)){ ## special behaviour of %X%
          all_df[i] <- 1
        } 
      }
    }
    if(any(all_df == "")){
      ret <- NULL
    }else{
      ret <- prod(all_df) # global df for bl is product of all df 
      if( identical(ret, numeric(0)) ) ret <- NULL
    } 
    return(ret)
  }
  
  #### get the specified df for each base-learner
  ## does not take into account base-learners that do not have brackets
  if(length(tmp) == 0){
    bl_df <- NULL
  }else{
    bl_df <- vector("list", length(tmp))
    bl_df[equalBrackets] <- lapply(tmp[equalBrackets], function(x) try(get_df(x)))
    bl_df <- unlist(bl_df[equalBrackets & (!sapply(bl_df, class) %in% "try-error")])
    #print(bl_df)
    
    if( !is.null(bl_df) && any(abs(bl_df - bl_df[1]) > .Machine$double.eps * 10^10) ){
      warning("The base-learners differ in the degrees of freedom.")
    }
    
    if(!is.null(bl_df)){
      df_argformula <- get_df(tfm) 
      df_effects <- min(bl_df)
    }
  }
  
  ####### put together the model formulae 
  xpart <- paste(as.vector(tmp), collapse = " + ")
  fm <- as.formula(paste("dresponse ~ ", xpart))
  
  ## find variables that are defined in environment(formula) but not in environment(fm) or in data 
  fm_vars <- all.vars(fm) # all variables of fm 
  
  ## vars_envir_formula <- fm_vars[ ! fm_vars %in% c(names(data), "dresponse" , "ONEx", "ONEtime", yind) ]
  # variables that exist in environment(fm) 
  vars1 <- sapply(fm_vars, exists, envir = environment(fm), inherits = FALSE)   
  if(!is.null(names(data))){
    vars2 <- sapply(fm_vars, function(x){ x %in% names(data)} ) # variables that exist in data
  }else{
    vars2 <- FALSE 
  }
  
  # variables that exist neither in environment(fm) nor in data... 
  vars_envir_formula <- fm_vars[ !(vars1 | vars2) ]
  # ... take those from the environment of the formula with which FDboost was called 
  for(i in seq_along(vars_envir_formula)){
    if(! exists(vars_envir_formula[i], envir = environment(formulaFDboost)))
      stop("Variable <", vars_envir_formula[i], "> does not exist.")
    tmp <- get(vars_envir_formula[i], envir = environment(formulaFDboost))
    assign(x = vars_envir_formula[i], value = tmp,  envir = environment(fm))
  }
  rm(tmp)
  
  ## TODO: implement weights
  # ### expand weights for observations
  # if (is.null(weights)) weights <- rep(1, nr)
  # w <- weights
  # if(is.array(response)){
  #   if (length(w) == ydim[1]) w <- rep(w, prod(ydim[2:3])) # expand weights if they are only on the columns
  #   if (length(w) == prod(ydim[1:2])) w <- rep(w, ydim[3]) # expand weights if they are only on the columns
  #   # check dimensions of w
  #   if(length(w) != prod(ydim)) stop("Dimensions of weights do not match the dimensions of the response.")   
  # }
  # 
  # ### set weights of missing values to 0
  # if(sum(is.na(dresponse)) > 0){
  #   w[which(is.na(dresponse))] <- 0
  # }
  # 
  # if(all(w == 0)) stop("All weights are zero!")
  browser()
  ### Model fit
  if (length(data) > 0 && !(is.list(data) && length(data) == 1 && is.null(data[[1]]))) {
    ### mboost isn't happy with nrow(data) == 0 / list(NULL)
    ret <- mboost(fm, data = data, offset = offset, ...) 
  } else {
    ret <- mboost(fm, offset = offset, ...)
  }
  
  ### Return results
  
  ### assign new class (e.g., for specialized predictions)
  class(ret) <- c("FDboost", class(ret))
  if(!is.null(id)) class(ret) <- c("FDboostLong", class(ret))
  if(scalarResponse) class(ret) <- c("FDboostScalar", class(ret))
  ## generate an id-variable for a regular response
  if(is.null(id)){
    if(scalarResponse){
      id <- 1:NROW(response)
    }else{
      id <- rep(1:ydim[1], times = ydim[2])
    }
  } 
  
  ### reset weights for cvrisk etc., expanding works OK in bl_lin_matrix!
  # ret$"(weights)" <- weights
  # <SB> do not reset weights as than the integration weights are lost
  
  ret$yname <- yname
  ret$ydim <- ydim  
  ret$yind <- time
  ret$data <- data
  ret$id <- id
  attr(ret$id, "nameid") <- nameid
  
  # save the call
  ret$call <- match.call()
  
  # save the evaluated call
  ret$callEval <- ret$call
  ret$callEval[-1] <- lapply(ret$call[-1], function(x){  
    eval(x, parent.frame(3)) # use the environment of the call to FDboost()
  })
  ret$callEval$data <- NULL # do not save data and weights in callEval to save memory
  
  # save formulas as character strings to save memory
  ret$argformula <- paste(deparse(argformula), collapse = "")
  if(scalarNoFLAM) ret$argformula <- ""
  ret$formulaFDboost <- paste(deparse(formulaFDboost), collapse = "")
  ret$formulaMboost <- paste(deparse(fm), collapse = "")
  
  ret
}