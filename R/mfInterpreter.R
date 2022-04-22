#' Extract structural info from obj.formula
#' 
#' Extract info from obj.formula relevant for the basic structure of the (response) object.
#'
#' @param x a formula.
#' @param varnames logcial, should character variable names be returned (default).
#' Alternatively formula parts are directly passed.
#' 
#' @return a named list containing the `value`, `dim`, `arg` and `id` 
#' extracted from the formula.
#' @export
#' @importFrom formula.tools lhs rhs
#'
mfInterpret_objformula <- function(x, varnames = TRUE, ...) {
  UseMethod("mfInterpret_objformula")
}

#' @importFrom formula.tools lhs rhs
mfInterpret_objformula <- function(x, varnames = TRUE, ...) {
  lf <- lhs(x)
  rf <- rhs(x)
  # list with formula parts
  v <- list()
  
  # interpret lhs
  if(length(lf) == 1) {
    v$value <- lf
    # v$dim <- ".dim_"
  } else if(length(lf) == 3) {
    if(lf[[1]] == "^") {
      v$value <- lf[[2]]
      v$dim <- lf[[3]]
    } 
  } else 
    stop("Left-hand side has to be a single variable or of the form 'value ^ dim'.")
  
  # interpret rhs
  if(rf[[1]] == "|") {
    v$args <- rf[[2]]
    v$id <- rf[[3]]
  } else {
    v$args <- rf
    v$id <- as.symbol(".id_")
  }
  
  if(varnames)
    v <- lapply(v, all.vars)
  v
}

#' Convert obj.formula and friends to FDboost formulae
#'
#' @param x an obj.formula needed for interpretation.
#' @param ... other arguments needed to interpret `x`.
#'
#' @return A list containing the `formula`, `timeformula` and `id`-formula passed
#' to FDboost 
#' 
#' @export
#'
mfFormulae_for_FDboost <- function(x, ...) {
  UseMethod("mfFormulae_for_FDboost")
}

#' Interpret obj.formula as FDboost timeformula
#' 
#' @param x an obj.formula of the form `value^dim ~ arg | id` reflecting the structure 
#' of the response. `^dim` and/or `|id` can be missing, when no dimensional/id
#' information has to be passed.
#' @param formula a formula describing the model on a meta level with for the response
#' object
#' @param ... additional arguments to be passed to or from methods.
#'
#' @importFrom formula.tools lhs rhs
mfFormulae_for_FDboost.formula <- function(x, formula, ...) {
  v <- mfInterpret_objformula(x, varnames = FALSE)
  
  # build timeformula from template
  timeformula <- ~ arg
  if(is.null(v$dim))
    timeformula[[2]] <- v$args else {
      timeformula <- update(timeformula, ~ brandom(dim, df = nlevels(dim)) %Xa0% (arg))
      timeformula[[2]][[2]][[2]] <- timeformula[[2]][[2]][[3]][[2]] <- v$dim
      timeformula[[2]][[3]] <- v$args
    }
  environment(timeformula) <- environment(x)
  
  yvalname <- v$value
  yidname <- v$id
  
  id <- ~ id
  id[[2]] <- v$id
  
  # change response name in formula to response value name
  lhs(formula) <- v$value
  
  list(
    formula = formula,
    timeformula = timeformula,
    id = id
  )
}
