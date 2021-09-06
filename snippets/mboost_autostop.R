# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# function for stopping mboost fit automated based on threshold criterion  #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# function arguments
## object: a fitted mboost object
## score.tol: numeric, threshold for the 99% quantile of the absolute value 
#             of the score functions (ngradient * design) of the baselearners.
#             To ignore this, set to NULL.
## rel.d.risk.tol: numeric, threshold for the absolute relative change in risk 
#                  from one iteration to another. NULL to ignore.
## d.risk.tol: numeric, threshold for the total (signed) risk change in risk
#              from one iteration to another. NULL to ignore.
## max.iter:   integer, maximal number of boosting iterations. NULL to ignore.
## contidion:  function, taking TRUE (below threshold) / FALSE (above threshold)
#              for score.tol, rel.d.risk.tol, d.risk.tol, max.iter as input to
#              return TRUE as long as a combined threshold is not reached.
#              Defaults, to all() stopping as soon as any of the thresholds is met.

autostop <- function( object, score.tol = 1e-6, rel.d.risk.tol = score.tol / 10, d.risk.tol = NULL, max.iter = 200, condition = all, verbose = TRUE ) {
  if(is.null(d.risk.tol)) d.risk.tol <- -Inf
  if(is.null(rel.d.risk.tol)) rel.d.risk.tol <- -Inf
  if(is.null(score.tol)) mean.grad.tol <- -Inf
  if(is.null(max.iter)) max.iter <- Inf
  
  designs <- lapply(object$baselearner, extract, "design")
  correct_format <- sapply(designs, function(X) 
    nrow(X) == length(object$response))
  if(!any(correct_format)) {
    warning("No baselearner has design matrix in appropriate format and, thus,
            score threshold is ignored.")
    score.tol <- -Inf } else {
      if(!all(correct_format)) message(
        paste( "Design matrix of", names(designs)[!correct_format], 
               "has wrong format and is ignored in score computation.\n", 
               collapse = ", ")
      )
    }
  correct_format <- which(correct_format)
  
  # initialize
  r <- tail(risk(object), 2)
  d.risk <- - diff(r)
  if(length(d.risk) == 0) d.risk <- Inf
  score <- if(length(correct_format)>0) 
    quantile( abs(unlist(
    lapply(designs[correct_format], function(X) as.matrix(crossprod(X, object$resid())))
    )), .99) else Inf
  mstp <- mstop(object)
  if(verbose)
    writeLines("\nBoosting iteration:")

  # threshold
  while( condition( 
          abs(score) > score.tol,
          abs(d.risk) > abs(r[1]) * rel.d.risk.tol, 
          d.risk > d.risk.tol,
          mstp <= max.iter ) ) {
    mstp <- mstp + 1
    r <- tail( risk(object[mstp]), 2 )
    d.risk <- - diff(r)
    score <- if(length(correct_format)>0) 
      quantile( abs(unlist(
        lapply(designs[correct_format], function(X) as.matrix(crossprod(X, object$resid())))
      )), .99) else Inf
    if(verbose) 
      cat(mstp, "|")
  }
  
  if(verbose) {
    reached <- function(x) paste(">", x, "reached")
    cat(
      if(d.risk <= d.risk.tol) reached("d.risk.tol"),
      if(d.risk <= abs(r[1]) * rel.d.risk.tol) reached("rel.d.risk.tol"),
      if(abs(score) <= score.tol ) reached("mean.grad.tol"),
      if(mstp >= max.iter) reached("max.iter")
    )
  }
  
  if(verbose) cat("\nIteration: ", mstp,
                           "\nRisk: ", r[2],
                           "\nRisk diff.: ", d.risk,
                           "\n99% score: ", score)
}
