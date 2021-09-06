
# FD-class ---------------------------------------------------------------

#' @importFrom methods setOldClass
#' @exportClass FD
#' @exportClass FD_regular
#' @exportClass FD_irregular

setOldClass("FD")
setOldClass("FD_regular")
setOldClass("FD_irregular")

#' An `FD` S3 class for FDboost data formats
#'
#' @description The `FD` class is an abstract class for lists
#' containing a sample of 'functional data' objects together 
#' with other variables carrying additional information. 
#' Objects of these class belong either to the subclass `FD_regular` or 
#' `FD_irregular` reflecting either of the formats desctibed in `?FDboost`.
#'
#' @name FD-class
#' @aliases FD_regular-class FD_irregular-class
#' @seealso [FDboost()]
#'    
NULL


# Convert to FD -----------------------------------------------------------

#' Convert other data to the FDboost format
#'
#' @param x the data object.
#' @param ... additional arguments passed to methods.
#'
#' @return an object of class `FD` with subclass `FD_regular` or `FD_irregular`.
#' #' @seealso [FDboost(), FD-class]
#' @export
#'
#' @examples
as_FD <- function(x, ...) {
  UseMethod("as_FD")
}


#' Convert data.frame to FDboost data format
#'
#' @param x a data.frame including the response values, functional response arguments,
#' dimension indicator for multidimensional response functions and an id column 
#' indicating the curve id
#' @param obj.formula a formula of the form `value^dim ~ arg | id` indicating the respective 
#' columns of `x`. Missing `^dim` or `|id` are interpreted as only one dimension/id
#' occuring in the data.
#' @param ... additional arguments to be passed to or from methods.
#' 
#' @details A data.frame will basically directly passed to FDboost without much further 
#' transformation. In principle, only the id-column is replaced by a numeric vector 
#' ranging from 1 to length(unique(id)), and the dim-colum is converted to factor.
#'
#' @return a data.frame in 'irregular' FDboost format `FD_irregular`.
#' @seealso [FDboost(), FD-class]
#' @export
#'
#' @examples
as_FD.data.frame <- function(x, obj.formula, ...) {
  v <- mfInterpret_objformula(obj.formula)
  
  if(any(names(x) == v$dim))
    x[[v$dim]] <- factor(x[[v$dim]])
  if(any(names(x) == v$id))
    x[[v$id]] <- as.numeric(factor(x[[v$id]], levels = unique(x[[v$id]])))
  class(x) <- c("FD_irregular", "FD", class(x))
  x
}

#' Interpret response given as tbl_cube for FDboost
#'
#' @param x a 'tbl_cube' containing the response values as measures and other variables
#' as dimensions.
#' @param obj.formula a formula of the form `value^dim ~ arg | id` reflecting the meaning
#' of the variables. Missing `^dim` or `|id` are interpreted as only one dimension/id
#' occuring in the data.
#' @param ... additional arguments to be passed to or from methods.
#'
#' @return a list with a matrix containing the values and vectors of length nrow(matrix)
#' for the other variables corresponding to the 'regular' FDboost data format
#' @seealso [FDboost(), tbl_cube, FD-class]
#' @export
#'
#' @examples
as_FD.tbl_cube <- function(x, obj.formula, ...) {
  v <- mfInterpret_objformula(obj.formula)
  # arrange array dimensions such that id is first and dim is last 
  ret <- x$dims
  pos <- c("dim" = which(names(ret)==v$dim), 
           "id" = which(names(ret)==v$id))
  pos <- c(pos["id"], setdiff(seq_along(ret), pos), pos["dim"])
  ret[pos[-1]] <- expand.grid(ret[pos[-1]]) # expand without id
  ret[[v$value]] <- matrix( aperm(x$mets[[v$value]], pos), nrow = length(ret[[v$id]]))
  # ensure correct format of dim and id
  if(any(names(ret) == v$dim))
    ret[[v$dim]] <- factor(ret[[v$dim]])
  if(any(names(ret) == v$id))
    ret[[v$id]] <- factor(ret[[v$id]], levels = ret[[v$id]])
  class(ret) <- c("FD_regular", "FD", class(ret))
  ret
}

#' Interpret response given as list for FDboost
#'
#' @param x either a list with an entry per response unit containing the 
#' response data or a list with one element containing the response data indicated 
#' by the left-hand side of the `formula`.
#' @param formula a formula indicating the element inside of `x` to be converted. 
#' If its left-hand side is NULL, `x` will be considered in FDboost format. If 
#' the response name is not contained in the list,
#' the list will be considered the response object itself. 
#' 
#' @param ... additional arguments to be passed to or from methods.
#' 
#' @details This interpreter applies mfInterpret_for_FDboost to all entries of the response 
#' data list, stacks them to a large data.frame and adds an id variable called 
#' '_id_'.
#'
#' @return a data.frame in 'irregular' FDboost format
#' @seealso FDboost
#' @export
#' @importFrom formula.tools lhs.vars
#'
#' @examples
as_FD.list <- function(x, formula, ...) {
  
  if(is.null(formula)) {
    # just assume x is already in FDboost format
    return(x)
  }
  
  # identify response name
  obj.name <- lhs.vars(formula)
  
  if(is.null(obj.name)) {
    # just assume x is already in FDboost format
    return(x)
  }
  
  if(obj.name %in% names(x)) {
    # if obj contained in list, go into object and apply method again
    xo <- x[names(x) != obj.name]
    x <- eval(formula[[2]], x)
    x <- as.list(as_FD(x, formula = formula, ...))
    x[names(xo)] <- xo
    return(x)
  } else {
    # hover over list entries
    x <- lapply(x, as_FD, formula = formula, ...)
    if(!all(sapply(x, inherits, "FD_irregular"))) {
      stop("FD conversion for list responses in regular data format
           not implemented, yet. Sorry.")
    }
    names(x) <- NULL
    x <- as.list(bind_rows(x, .id = ".id_"))
    x[[".id_"]] <- as.numeric(x[[".id_"]])
    class(x) <- c("FD_irregular", "FD", class(x))
    return(x)
  }
}



