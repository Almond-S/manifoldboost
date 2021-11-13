
setOldClass("RiemannL2sim")

#' R6 class for a simulation of RiemannL2 models
#' 
#' @description 
#' This class is created to set up a simulations scenario, simulate data from it,
#' fit models, do cross-validation, predict quantities of interest, and summarize
#' relevant statistics in an easy and transparent way.
#' 
#' @field model0 the underlying true model object.
#' @field newdata0_FDboost the data used for sampling in FDboost format.
#' @field new_mf0 a clone of the mfGeometry of model0 with newdata0.
#' @field new_pred0_ predictions of model0 on newdata0.
#' @field sessionInfo the result of sessionInfo() executed at initialization.
#' @field new_effect0 effects of model0 evaluated on newdata0.
#' @field new_fac0 tensor-product factorization of model0 on newdata0.
#' @field seed the random seed used for the simulation.
#' @field n a vecor of length three of the form c(n_sim = NA, n_obj = NA, n_grid = NA)
#' containing the specified (average) sample sizes.
#' @field simdata the simulated dataset.
#' @field cluster the cluster id used for parallel computing.
#' @field family the family object used.
#' @field mfboost_control the list specified for controlling mboost in the model fits.
#' @field models the fitted model objects.
#' @field runtime the runtime of the model fits.
#' @field warnings warnings produced during model fits.
#' @field cvrisk_control a list of argument supplied to control the cross-validation.
#' @field cv_seeds seeds used for cross-validation.
#' @field cvs cross-validation objects.
#' @field runtime_cv runtimes of the cross-validations.
#' @field warnings_cv warnings produced during cross-validations.
#' @field preds_ predictions of model fits in internal geometry format.
#' @field poles_ poles predicted in internal geometry format.
#' @field effects effects estimated on fitted models.
#' @field facs factorized model fits.
#' @field simsum first simulation summary.
#' @field bot TGbot object of R package telegram used for automized notifications.
#' @field name name of the simulation.
#' 
#' @details 
#' The methods depend on the magrittr \code{%>%} operator and \code{dplyr::bind_rows}. 
#' And, if messaging is desired, on the R-package \code{telegram}.
#' 
#' @import R6 
#' @import Matrix
#' @importFrom formula.tools lhs
RiemannL2sim <- R6Class("RiemannL2sim",
                  public = list(
                    
                    #' @description initialize simulation scenario.
                    #' @param model0 the model object with the underlying truth.
                    #' @param newdata0 new data to simulate from. The default 
                    #' \code{NULL} will take the data of \code{model0}.
                    #' @param name character string, a name for the simulation scenario.
                    initialize = function(model0, newdata0 = NULL, name = NULL) {
                      
                      # set simulation name
                      if(!is.null(name))
                        self$name <- name
                      
                      self$sessionInfo <- sessionInfo()
                      self$model0 <- model0
                      
                      # initialize geometry for newdata0
                      formula <- model0$formula
                      obj.formula <- model0$obj.formula
                      private$variablenames <- v <- mfInterpret_objformula(obj.formula)
                      self$newdata0_FDboost <- data_FDboost <- as_FD(newdata0, 
                                            formula = formula, 
                                            obj.formula = obj.formula)
                      
                      self$new_mf0 <- mf <- model0$family@mf$clone(deep = TRUE)
                      mf$initialize(data = data_FDboost, formula = obj.formula)
                      
                      # get predictions for newdata0
                      suppressWarnings(
                        pred0 <- predict(model0, 
                                       newdata = newdata0, 
                                       type = "response")
                        )
                      self$new_pred0_ <- pred0 %>% mf$structure() %>% mf$register()
                      
                      # ... and the pole of the regression model
                      suppressWarnings(
                      mf$pole_ <- predict(model0, newdata = newdata0, 
                                          type = "response", which = 0) %>%
                        mf$structure()
                      )
                      # compute residuals
                      fam <- RiemannL2(mf)
                      eps <- fam@ngradient(y = data_FDboost[[v$value]], f = pred0) %>%
                        split(data_FDboost[[v$id]])
                      
                      # now, get basis representations of residuals...
                      
                      ## factorize to get response baselearner
                      suppressWarnings(
                        self$new_fac0 <- factorize(model0, newdata = newdata0)
                        )
                      
                      ## get design matrix for new data
                      e <- environment(self$new_fac0$resp$baselearner[[1]]$dpp)
                      # attach new tangent space normals or similar
                      data_FDboost <- c(data_FDboost, 
                                        attr(fam@update_formula(~0, pole_ = mf$pole_),
                                             "update_formula_vars"))
                      B <- e$newX(newdata = data_FDboost, prediction = TRUE)
                      K <- B$K
                      B <- B$X
                      
                      ## break down to single design matrices per curve
                      ids <- split(1:length(data_FDboost[[v$id]]), 
                                   data_FDboost[[v$id]])
                      private$new_B <- B <- lapply(ids, function(i) B[i, , drop = FALSE])
                      
                      ## compute basis representations
                      ## (with slight penalty against singularity problems)
                      theta <- Map(function(b, y) 
                        solve(crossprod(b) + 1e-10*K, crossprod(b, y)), B, eps)
                      theta <- lapply(theta, as.vector)
                      
                      ## center residual coefs
                      mean_theta <- do.call(cbind, theta) %>% rowMeans() 
                      private$theta <- lapply(theta, `-`, mean_theta)
                      
                      # also evaluate effects on predictor level for later
                      suppressWarnings(
                        self$new_effect0 <- predict(
                                 object = self$model0,
                                 which = self$model0$which(),
                                 newdata = newdata0, 
                                 type = "link")
                      )
                      # done
                      invisible(self)
                    },
                    model0 = NULL,
                    newdata0_FDboost = NULL,
                    new_mf0 = NULL,
                    new_pred0_ = NULL,
                    sessionInfo = NULL,
                    new_effect0 = NULL, 
                    new_fac0 = NULL,
                    
                    #' @description function sampling from the evaluations of
                    #' a single response observation
                    #' @param x vector to be subsampled from
                    #' @param size average sample size
                    gridsampler = function(x, size) {
                      stopifnot(size >= 3)
                      idx <- seq_along(x)
                      fx <- sample(idx, 3)
                      if(size > 3) {
                        idx <- setdiff(idx, fx)
                        n <- length(idx)
                        fx <- c(fx, idx[runif(n-3) < (size-3)/(n-3)])
                      }
                      x[fx]
                    },
                    #' @description function for sampling observations
                    #' @param n_sim,n_obj,n_grid sample sizes
                    #' @param seed random seed
                    #' @param sample_obj function to sample from objects.
                    #' @param sample_grid function to sample within objects.
                    sample_topdown = function(n_sim = 1, n_obj = 1, n_grid = 40, 
                                      seed = NULL, 
                                      sample_obj = sample,
                                      sample_grid = self$gridsampler) {
                      if(!is.null(seed))
                        self$seed <- seed
                      if(is.null(self$seed))
                        self$seed <- sample(seq(1:99999), size = 1)
                      set.seed(self$seed)
                      
                      self$n <- c(n_sim = n_sim, n_obj = n_obj, n_grid = n_grid)
                      
                      sample_dataset <- function(i) {
                        mf <- self$new_mf0$clone(deep = TRUE)
                        v <- private$variablenames
                        data_FDboost <- self$newdata0_FDboost
                        
                        # randomly distribute residuals
                        eps_ <- Map(`%*%`, 
                                   private$new_B, 
                                   sample(private$theta, replace = TRUE)) %>% 
                          lapply(as.vector) %>%
                          unsplit(data_FDboost[[v$id]]) %>%
                          mf$structure() %>% mf$register_v(y0_ = mf$pole_)
                        
                        # store residual variances
                        residual_var <- mean(Re(unlist(mf$innerprod(eps_))))
                        
                        # transport smooth residuals to model predictions
                        mf$y_ <- eps_ %>% mf$transport(mf$pole_, self$new_pred0_) %>% 
                          mf$exp(y0_ = self$new_pred0_)
                        
                        # subsample observations (ids)
                        ids <- sort(sample_obj(seq_along(eps_), n_obj))
                        mf$slice(ids)
                        # subsample grid (inner ids)
                        mf$y_ <- lapply(mf$y_, function(y_) {
                          idx <- sample_grid(seq_along(y_), n_grid)
                          y_[-idx] <- NA
                          y_
                        })
                        
                        # build new dataset
                        lens <- lengths(data_FDboost)
                        leny <- length(data_FDboost[[v$value]])
                        data0 <- list(
                          resp = as.data.frame(
                            data_FDboost[lens == leny]),
                          cov = as.data.frame(
                            data_FDboost[lens != leny])
                        )
                        # remove unselected ids
                        data0$resp <- data0$resp[data_FDboost[[v$id]] %in% ids, ]
                        data0$resp[[v$value]] <- mf$unstructure(mf$y_)
                        # remove unselected grid points
                        data0$resp <- data0$resp[!is.na(data0$resp[[v$value]]), ]
                        # re-define ids for subdata
                        data0$resp[[v$id]] <- as.numeric(factor(data0$resp[[v$id]]))
                        #remove unselected ids from cov data
                        data0$cov <- data0$cov[unique(ids), ]
                        
                        ret <- unlist(setNames(data0, NULL), FALSE)
                        attr(ret, "residual_var") <- residual_var
                      }
                  
                      self$simdata <- lapply(seq_len(n_sim), sample_dataset)
                      
                      invisible(self)
                    },
                    #' @description function for sampling observations in a blockwise
                    #' strategy to preserve covariate composition.
                    #' @param n_sim,n_grid sample sizes.
                    #' @param n_blocks number of blocks in each sample.
                    #' @param seed random seed. For the defaul \code{NULL}, the seed
                    #' is randomly drawn.
                    #' @param random_trafo function taking \code{y_} an returning it
                    #' after random transformations.
                    #' @param sample_grid function to sample within objects.
                    sample_blockwise = function(n_sim = 1, n_blocks = 1, n_grid = 40, 
                                              seed = NULL,
                                              # function returning a randomized version of mf$y_
                                              # e.g. with respect to translation / scale / rotation
                                              random_trafo = function(y_) y_, 
                                              sample_grid = self$gridsampler) {
                      if(!is.null(seed))
                        self$seed <- seed
                      if(is.null(self$seed))
                        self$seed <- sample(seq(1:99999), size = 1)
                      set.seed(self$seed)
                      
                      n_obj <- max(self$model0$id) * n_blocks
                      
                      self$n <- c(n_sim = n_sim, n_obj = n_obj, n_grid = n_grid)
                      
                      suppressWarnings(fac0 <- factorize(self$model0))
                      
                      B <- extract(fac0$resp, "design")[[1]]
                      ## break down to single design matrices per curve
                      ids <- split(1:nrow(B), 
                                   self$model0$id)
                      B <- lapply(ids, function(i) B[i, , drop = FALSE])
                      
                      pred0_ <- predict(self$model0, type = "response") %>% 
                        self$model0$family@mf$structure()
                      
                      v <- private$variablenames
                      
                      sample_block <- function(i) {
                        mf <- self$model0$family@mf$clone(deep = TRUE)
                        
                        # randomly distribute residuals
                        eps_ <- Map(`%*%`, 
                                    B, 
                                    sample(private$theta, 
                                           size = n_obj / n_blocks, 
                                           replace = TRUE)) %>% 
                          lapply(as.vector) %>%
                          unsplit(self$model0$id) %>%
                          mf$structure() %>% mf$register_v(y0_ = mf$pole_)
                        
                        # store residual variances
                        residual_var <- mean(Re(unlist(mf$innerprod(eps_))))
                        
                        # transport smooth residuals to model predictions
                        mf$y_ <- eps_ %>% mf$transport(mf$pole_, pred0_) %>% 
                          mf$exp(y0_ = pred0_)
                        
                        # subsample grid (inner ids)
                        mf$y_ <- lapply(mf$y_, function(y_) {
                          idx <- sample_grid(seq_along(y_), n_grid)
                          y_[-idx] <- NA
                          y_
                        })
                        
                        # build new dataset
                        data0 <- c(self$model0$yind, self$model0$data)
                        lens <- lengths(data0)
                        leny <- length(self$model0$response)
                        data0 <- list(
                          resp = as.data.frame(
                            data0[lens == leny]),
                          cov = as.data.frame(
                            data0[lens != leny])
                        )
                        # randomize resp with respect to invariances
                        mf$y_ <- random_trafo(mf$y_)
                        
                        # add response
                        data0$resp[[v$value]] <- mf$unstructure(mf$y_)
                        # remove unselected grid points
                        data0$resp <- data0$resp[!is.na(data0$resp[[v$value]]), ]
                        
                        attr(data0, "residual_var") <- residual_var
                        data0
                        # unlist(setNames(data0, NULL), FALSE)
                      }
                      
                      # do sampling 
                      for(i in 1:n_sim) {
                        blocks <- lapply(seq_len(n_blocks), sample_block)
                        simdata <- list()
                        simdata$resp <- lapply(blocks, `[[`, "resp") %>% bind_rows(.id = "block_id")
                        simdata$cov <- lapply(blocks, `[[`, "cov") %>% bind_rows()
                        
                        simdata$resp[[v$id]] <- as.numeric(interaction(
                          simdata$resp[[v$id]], simdata$resp$block_id))
                        simdata$resp$block_id <- NULL
                        
                        simdata <- unlist(setNames(simdata, NULL), FALSE)
                        attr(simdata, "residual_var") <- mean(sapply(blocks, attr, "residual_var"))
                        self$simdata[[i]] <- simdata 
                      } 
                      
                      invisible(self)
                    },
                    seed = NULL,
                    n = c(n_sim = NA, n_obj = NA, n_grid = NA),
                    simdata = NULL,
                    
                    #' @description initialize parallel computing
                    #' @param nc number of cores.
                    parallel_setup = function(nc = detectCores()) {
                      self$cluster <- makePSOCKcluster(1:nc)
                      invisible(
                        clusterEvalQ(cl = self$cluster, library(shapeboost))
                        )
                    },
                    cluster = NULL, 
                    
                    #' @description end parallel computing
                    parallel_stop = function() {
                      stopCluster(self$cluster)
                      self$cluster <- NULL
                    },
                    
                    #' @description fit simulated datasets
                    #' @param family mfFamily object used for fitting.
                    #' @param ... arguments supplied to mboost control argument.
                    #' @param verbose logical, should live info be printed.
                    fit = function(family = self$model0$family,
                      ...,
                      verbose = FALSE) {
                      
                      self$mfboost_control <- list(...)
                      self$family <- family
                      
                      self$models <- list()
                      
                      # modify formula for direct FDboost data format
                      new_formula <- self$model0$formula
                      lhs(new_formula) <- NULL
                      
                      # prepare counter
                      i <- 1
                      if(!is.null(self$cluster))
                        clusterExport(self$cluster, "i", environment())
                      
                      fit_dat <- function(dat, family, ...) {
                        
                        e <- environment(family@ngradient)
                        fam <- RiemannL2(
                          mf = family@mf$clone(deep = TRUE),
                          pole.type = e$pole.type,
                          pole.control = e$pole.control)
                        
                        warn <- list()
                        if(verbose) cat("Fit model", i, "\n")
                        withCallingHandlers(
                          runtime <- system.time(
                              mod <- try(mfboost( 
                              data = dat, family = fam, ...))
                            ),
                          warning = function(w) {
                            warn[[length(warn)+1]] <<- conditionMessage(w)
                          }
                          )
                        if(verbose) print(runtime)
                        i <<- i+1
                        
                        attr(mod, "warnings") <- warn
                        attr(mod, "system.time") <- runtime
                        mod
                      }
                      
                      self$models <- private$my_lapply(self$simdata, fit_dat, 
                                               formula = new_formula,
                                               obj.formula = self$model0$obj.formula,
                                               family = family, ...)
                      
                      self$runtime <- self$warnings <- list()
                      for(j in seq_along(self$models)) {
                        self$runtime[[j]] <- attr(self$models[[j]], "system.time")
                        self$warnings[[j]] <- attr(self$models[[j]], "warnings")
                        attr(self$models[[j]], "system.time") <- NULL
                        attr(self$models[[j]], "warnings") <- NULL
                      }
                      
                      invisible(self)
                    },
                    family = NULL,
                    mfboost_control = list(),
                    models = list(),
                    runtime = list(),
                    warnings = list(),
                    
                    #' @description conduct cross-validations of fitted models.
                    #' @param type type of resampling (as for mboost).
                    #' @param B number of resampling folds.
                    #' @param ... other arguments supplied to cvLong.
                    #' @param verbose logical, should live info be printed?
                    #' @param seeds vector of random seeds for cross-validation. For the 
                    #' default \code{NULL}, seeds are randomly drawn.
                    crossvalidate = function(type = "kfold", B = 5, ..., verbose = FALSE, seeds = NULL) {
                      
                      self$cvrisk_control <- list(type = type, B = B, ...)
                      # specify seeds
                      if(!is.null(seeds)) {
                        if(length(seeds) == 1)
                          seeds <- rep(seeds, length(self$models))
                        self$cv_seeds <- seeds
                      }
                      if(is.null(self$cv_seeds))
                        self$cv_seeds <- sample(seq(1:99999), size = length(self$models))
                      stopifnot(length(self$cv_seeds) == length(self$models))
                      
                      # prepare counter
                      i <- 1
                      if(!is.null(self$cluster)) {
                        clusterExport(self$cluster, list("i", "type", "B"), environment())
                      }
                      
                      self$cvs <- list()
                      
                      cv_mod <- function(mod, seed, ...) {
                        set.seed(seed)
                        warn <- list()
                        if(verbose) cat("CV model", i, "\n")
                        withCallingHandlers(
                          runtime <- system.time(
                          cv <- try(FDboost:::cvrisk.FDboost(mod, 
                                                             folds = cvLong(
                                                                            id = mod$id, 
                                                                            weights = mod$`(weights)`, 
                                                                            type = type, 
                                                                            B = B), 
                                                             grid = 0:mstop(mod),...)
                            )),
                          warning = function(w) {
                            warn[[length(warn)+1]] <<- conditionMessage(w)
                          }
                          )
                        
                        if(verbose) print(runtime)
                        # increase counter
                        i <<- i+1
                        
                        attr(cv, "warnings") <- warn
                        attr(cv, "system.time") <- runtime
                        cv
                      }
                      
                      self$cvs <- private$my_Map(cv_mod, 
                                                 mod = self$models, 
                                                 seed = self$cv_seeds, 
                                                 MoreArgs = list(...))
                      
                      # apply cv
                      self$runtime_cv <- self$warnings_cv <- list()
                      for(j in seq_along(self$cvs)) {
                        try({
                          self$runtime_cv[[j]] <- attr(self$cvs[[j]], "system.time")
                          self$warnings_cv[[j]] <- attr(self$cvs[[j]], "warnings")
                          attr(self$cvs[[j]], "system.time") <- NULL
                          attr(self$cvs[[j]], "warnings") <- NULL
                          
                          mstop(self$models[[j]]) <- mstop(self$cvs[[j]])
                        })
                      }
                      
                      invisible(self)
                    },
                    cvrisk_control = list(),
                    cv_seeds = NULL,
                    cvs = list(),
                    runtime_cv = list(),
                    warnings_cv = list(),
                    
                    #' @description predict model0 and models and extract estimated effects.
                    predict = function() {
                      ## end-points of interest:
                      ## (all evaluated on self$data0_FDboost)
                      
                      mf <- self$model0$family@mf
                      newdata <- c(self$model0$yind, self$model0$data)
                      newdata[[self$model0$yname]] <- self$model0$response
                      
                      suppressWarnings({
                      
                      # response-level predictions
                        self$preds_ <- self$models %>% 
                          lapply(predict, 
                                 newdata = newdata, 
                                 type = "response") %>% 
                          lapply(mf$structure) %>% 
                          lapply(mf$register)
                      
                      # pole predictions
                        self$poles_ <- self$models %>% 
                          lapply(predict,
                                 which = 0,
                                 newdata = newdata,
                                 type = "response") %>% 
                          lapply(mf$structure) %>% 
                          lapply(mf$register)
                      
                      # single predictor-level effects
                        self$effects <- self$models %>% lapply(predict, 
                                   which = self$model0$which(),
                                   newdata = newdata,
                                   type = "link")
                      
                      # response directions
                      self$facs <- lapply(self$models, factorize, newdata = newdata)
                      
                      })
                      
                      invisible(self)
                    },
                    preds_ = NULL, # compare to pred0_ 
                    poles_ = NULL, # compare mf0$pole_
                    effects = NULL,
                    facs = NULL,
                    
                    #' @description generate simulation summary
                    summarize = function() {
                      
                      mf <- self$model0$family@mf
                      pred0_ <- predict(self$model0, type = "response") %>%
                        mf$structure()
                      suppressWarnings({
                        effect0 <- predict(
                          object = self$model0,
                          which = self$model0$which(), 
                          type = "link")
                        fac0 <- factorize(self$model0)
                      })
                      
                      # store sample sizes
                      self$simsum$n <- self$n
                      self$simsum$n_grids <- lapply(self$models, function(m) {
                        table(m$id) %>% as.numeric() / 2 
                      })
                      
                      # store selected mstops
                      self$simsum$mstops <- sapply(self$models, mstop)
                      
                      # compute predictor variance
                      self$simsum$var0 <- mean(
                        Re(mf$distance(
                          mf$pole_, pred0_))^2
                      )
                      
                      # store residual variance
                      self$simsum$residual_var <- sapply(self$simdata, attr, "residual_var")
                      
                      # evaluate predictions 
                      self$simsum$mse_predictions <- sapply(self$preds_, function(p_) 
                        mean(
                          Re(mf$distance(
                            p_, pred0_))^2
                      ))
                      
                      # evaluate poles
                      self$simsum$mse_poles <- sapply(self$poles_, function(p_) 
                        mean(
                          Re(mf$distance(
                            p_, mf$pole_))^2
                      ))
                      
                      # evaluate effects
                      self$simsum$mse_effects <- list()
                      
                      for(bl in seq_along(self$model0$baselearner)) {
                        eff0 <- effect0[, bl] %>% mf$structure()
                        # transport effects to common tangent space
                        effs <- Map(function(x, p_) {
                          x[, bl] %>% mf$structure() %>% 
                            mf$register_v(p_) %>% 
                            mf$transport(p_, mf$pole_)
                        }, self$effects, self$poles_)
                        # compute mse
                        mse <- sapply(effs, function(ef) mean(
                          Mod(unlist(mf$innerprod(ef))) - 
                            2*Mod(unlist(mf$innerprod(ef, eff0)))))
                        mse <- mse + mean(Re(unlist(mf$innerprod(eff0))))
                        self$simsum$mse_effects[[bl]] <- mse
                      }
                      names(self$simsum$mse_effects) <- extract(self$model0, "bnames", 
                                                                which = self$model0$which())
                      
                      # compare effect directions
                      self$simsum$angle_dirs <- list()
                      dir0 <- predict(fac0$resp, 
                                      which = seq_along(fac0$resp$baselearner), 
                                      type = "link")
                      dirs <- lapply(self$facs, function(fc) {
                        dir <- predict(fc$resp, 
                                       which = seq_along(fac0$resp$baselearner), 
                                       type = "link")
                      } )
                      
                      for(bl in seq_along(fac0$resp$baselearner)) {
                        this_dir0 <- dir0[, bl] %>% mf$structure()
                        # transport effects to common tangent space
                        this_dirs <- Map(function(x, p_) {
                          x[, bl] %>% mf$structure() %>% 
                            mf$register_v(p_) %>% 
                            mf$transport(p_, mf$pole_)
                        }, dirs, self$poles_)
                        # compute mse
                        ip <- sapply(this_dirs, function(ef) sum(
                            Mod(unlist(mf$innerprod(ef, this_dir0)))
                            ) / sqrt(sum(Mod(unlist(mf$innerprod(ef)))))
                        )
                        ip <- ip / sqrt(sum(Re(unlist(mf$innerprod(this_dir0)))))
                        self$simsum$angle_dirs[[bl]] <- acos(ip) / pi * 180
                      }
                      names(self$simsum$angle_dirs) <- extract(fac0$resp, "bnames", 
                                                               which = fac0$resp$which())
                      
                      # list run times
                      self$simsum$runtime <- list(
                        fit = sapply(self$runtime, `[`, 3),
                        cv = sapply(self$runtime_cv, `[`, 3)
                        )
                      
                      # extract variable importance measures
                      self$simsum$varimp <- c(
                        list("0" = varimp(self$model0)),
                        setNames(lapply(self$models, varimp), 
                                 seq_along(self$models))) %>% 
                        lapply(as.data.frame) %>% bind_rows(.id = "run")
                      
                      invisible(self)
                    },
                    simsum = NULL,
                    
                    #' @description send automated notifications via telegram using
                    #' the R package \link{telegram}. For getting started see their
                    #' package vignette.
                    #' @param token telegram bot token, see \code{?telegram::bot_token}.
                    #' @param id telegram user id, see \code{?telegram::user_id}.
                    #' @param name simulation name attached to the notification.
                    notify = function(token = bot_token('RBot'), id = user_id("me"), name = NULL) {
                      require(telegram)
                      ## set up telegram notifications
                      if(is.null(self$bot)) {
                        self$bot <- TGBot$new(token)
                        self$bot$set_default_chat_id(id)
                      }
                      if(!is.null(name)) 
                        self$name <- name
                      
                      self$bot$sendMessage(paste0("RiemannL2Sim\n", 
                                                 self$name, "\n(", 
                                                 paste(rbind(paste0(names(self$n),":"), self$n), 
                                                       collapse = " "),
                                                 ")"))
                      },
                    bot = NULL,
                    name = NULL,
                    #' @description make first summary for first inspection of results.
                    plot = function() {
                      require(ggplot2)
                      require(gridExtra)
                      require(ggExtra)
                      
                      pl <- list()
                      
                      s <- self$simsum
                      n_sim <- length(s$mse_predictions)
                      
                      pl$resid <- ggplot(data = NULL, aes(x = seq(-.1, .1, len = n_sim), y = s$residual_var / s$var0)) + 
                          geom_boxplot() + geom_text(aes(label = 1:n_sim)) + 
                          ylab("Residual variance / predictor variance") + 
                          xlim(-.3, .3) +
                          theme(axis.title.x = element_blank(), axis.text.x = element_blank())
                      
                      pl$pred <- ggplot(data = NULL, aes(x = seq(-.1, .1, len = n_sim), y = s$mse_predictions / s$var0)) + 
                          geom_boxplot() + geom_text(aes(label = 1:n_sim)) + 
                          ylab("MSE / predictor Variance") + xlim(-.3, .3) +
                          theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
                          ggtitle("Predicted vs. Original Mean")
                      pl$pole <- ggplot(data = NULL, aes(x = seq(-.1, .1, len = n_sim), y = s$mse_poles / s$var0)) + 
                          geom_boxplot() + geom_text(aes(label = 1:n_sim)) + 
                          ylab("MSE / predictor Variance") + xlim(-.3, .3) +
                          theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
                          ggtitle("Predicted vs. Original Pole")
                      
                      effmse <- as.data.frame(s$mse_effect)
                      effmse$run <- factor(1:n_sim)
                      effmse <- pivot_longer(effmse, cols = -which(names(effmse)=="run"), 
                                             names_to = "Baselearner", values_to = "MSE")
                      effmse$Baselearner <- factor(effmse$Baselearner, levels = unique(effmse$Baselearner))
                      levels(effmse$Baselearner) <- paste(seq_along(levels(effmse$Baselearner)), 
                                                          str_sub(levels(effmse$Baselearner), 1, 20), sep = ": ")
                      pl$effects <- ggplot(data = effmse, aes(x = scale(as.numeric(run)) / 8, 
                                                               y = MSE / s$var0,
                                                               col = Baselearner)) + 
                          geom_boxplot() + geom_text(aes(label = run)) + 
                          ylab("MSE / predictor Variance") + xlim(-.3, .3) +
                          theme(axis.title.x = element_blank(), 
                                axis.text.x = element_blank(),
                                # legend.position = "bottom", legend.direction = "vertical"
                          ) +
                          ggtitle("Predicted vs. Original Effects")
                      
                      dirangle <- as.data.frame(s$angle_dirs)
                      dirangle$run <- factor(1:n_sim)
                      dirangle <- pivot_longer(dirangle, cols = -which(names(dirangle)=="run"), 
                                               names_to = "Direction", values_to = "angle")
                      dirangle$Direction <- factor(dirangle$Direction, levels = unique(dirangle$Direction))
                      levels(dirangle$Direction) <- names(s$angle_dirs)
                      levels(dirangle$Direction) <- paste(seq_along(levels(dirangle$Direction)), 
                                                          str_extract(levels(dirangle$Direction), "\\[(.*?)\\]"), sep = ": ")
                      pl$dirs <- ggplot(data = dirangle, aes(x = scale(as.numeric(run)) / 8, 
                                                              y = angle,
                                                              col = Direction)) + 
                          geom_boxplot() + geom_text(aes(label = run)) + 
                          ylab("Angle in degree") + xlim(-.3, .3) +
                          theme(axis.title.x = element_blank(), 
                                axis.text.x = element_blank(),
                                # legend.position = "bottom", legend.direction = "vertical"
                          ) +
                          ggtitle("Angle to original response directions")
                      
                      rt <- as.data.frame(s$runtime)
                      rt$run <- seq_len(n_sim)
                      pl$runtime <- ggplot(data = rt, aes(x = fit, y = cv)) + 
                        geom_point() + ylab("Cross-validation time [sec]") + xlab("Initial fitting time [sec]") +
                        geom_text(aes(label = run), hjust = "right", vjust = "right")
                      pl$runtime <- ggMarginal(pl$runtime, type = "boxplot")
                      
                      sum(unlist(s$runtime))
                      # arrange common plot
                      grid.arrange(arrangeGrob(grobs = pl, 
                                               layout_matrix = matrix(seq_along(pl), ncol = 2), 
                                               widths = c(1/3, 2/3), 
                                               top = paste0("RiemannL2Sim ", 
                                                             self$name, " (", 
                                                             paste(rbind(paste0(names(self$n),":"), self$n), 
                                                                   collapse = " "),
                                                            ")")))
                    }
                    ), private = list(
                      theta = NULL,
                      new_B = NULL,
                      variablenames = NULL,
                      my_lapply = function(X, FUN, ...) {
                        if(is.null(self$cluster)) 
                            suppressWarnings(
                              suppressMessages(
                                lapply(X, FUN, ...)
                              )
                            ) else 
                              parLapply(cl = self$cluster, X, FUN, ...)
                      },
                      my_Map = function(FUN, ..., MoreArgs = NULL) {
                        if(is.null(self$cluster)) 
                          suppressWarnings(
                            suppressMessages(
                              mapply(FUN, ..., MoreArgs = MoreArgs, SIMPLIFY = FALSE)
                            )
                          ) else 
                            clusterMap(cl = self$cluster, FUN, ..., MoreArgs = MoreArgs)
                      }
                    ))
