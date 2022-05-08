setOldClass("RiemannL2sim_elastic")

RiemannL2sim_elastic <- R6Class("RiemannL2sim_elastic", inherit = RiemannL2sim,
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
                            self$newdata0_FDboost <- as_FD(newdata0, 
                                                           formula = formula, 
                                                           obj.formula = obj.formula)
                            
                            # done
                            invisible(self)
                          },
                          sample_bootstrap = function(n_sim = 1, n_obj = 10, per = NULL,
                                                    n_grid = 40, seed = NULL,
                                                      # function returning a randomized version of mf$y_
                                                      # e.g. with respect to translation / scale / rotation
                                                      random_trafo = function(y_) y_, 
                                                      sample_grid = self$gridsampler) {
                            if(!is.null(seed))
                              self$seed <- seed
                            if(is.null(self$seed))
                              self$seed <- sample(seq(1:99999), size = 1)
                            set.seed(self$seed)
                            
                            self$n <- c(n_sim = n_sim, n_obj = n_obj, n_grid = n_grid)
                            v <- private$variablenames
                            
                            # prepare (stratified) subsampling of observations
                            if(is.null(per)) {
                              ids_selected <- lapply(1:n_sim, function(i) 
                                sample(seq_along(pred0_), size = n_obj))
                            } else {
                              ids_selected <- lapply(1:n_sim, function(i) {
                                ids <- split(seq_along(pred0_), self$newdata0_FDboost[[per]])
                                sapply(ids, sample, size = n_obj)
                              })
                            }
                            
                            B <- private$new_B
                            
                            # sample observations
                            sample_obs <- function(ids) {
                              mf <- self$new_mf0$clone(deep = TRUE)
                              
                              # randomly distribute residuals
                              eps_ <- Map(`%*%`, 
                                          B, 
                                          sample(private$theta, # bootstrap errors
                                                 size = length(B), 
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
                              data0 <- self$newdata0_FDboost
                              lens <- lengths(data0)
                              data0 <- as.data.frame(data0[lens == length(mf$y_)])
                              data0 <- as.list(data0[ids, ])
                              
                              # randomize resp with respect to invariances
                              mf$y_ <- random_trafo(mf$y_)
                              
                              # add response
                              data0[[v$value]] <- t(sapply(mf$y_[ids], function(x) 
                                c(Re(x[!is.na(x)]), Im(x[!is.na(x)])) ))
                              data0[[v$arg]] <- rep(seq(
                                min(self$newdata0_FDboost[[v$arg]]), 
                                max(self$newdata0_FDboost[[v$arg]]), len = n_grid), 2)
                              
                              mf$slice(ids)
                              data0$resp[[v$value]] <- mf$unstructure(mf$y_)
                              # remove unselected grid points
                              data0$resp <- data0$resp[!is.na(data0$resp[[v$value]]), ]
                              
                              attr(data0, "residual_var") <- residual_var
                              data0
                              # unlist(setNames(data0, NULL), FALSE)
                            }
                            
                            # do sampling 
                            self$simdata <- lapply(ids_selected, sample_obs)
                            
                            invisible(self)
                          }
                        ))
