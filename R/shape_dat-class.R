
# helper functions ------------------------------------------------------


#' Match variables from formula to template
#' 
#' Internal functions matching the variables in `formula` to those provided in a `template`.
#' If formulas do not have the same form, the functions returns an error. 
#'
#' @param formula,template formulas of the same form, i.e., equal up to their 
#' variable names. 
#'
#' @return A named character vector with the variable names in formula as values
#' and the corresponding variable names in the template as names.
#'
match_template <- function(formula, template) {
  fvars <- all.vars(formula)
  tvars <- all.vars(template)
  fall <- all.names(formula)
  # remove anything which is not ^, |, cbind, or ~
  fall <- fall[fall %in% c(fvars, "^", "|", "cbind", "~")]
  tall <- all.names(template)
  if(length(fall) != length(tall))
    stop("Number of all.names() in formula differs from template.")
  names(fall) <- names(tall) <- tall
  tall <- replace(tall, tvars, fvars)
  eq <- all.equal(fall, tall)
  if(!eq)
    stop(paste("all.names(formula) don't differ from template in the sense that:", eq))
  fall[tvars]
}


#' Extract data with structure formula
#' 
#' Creates a object data data.frame from a formula stating the structure of the 
#' object according to a formula template. Used in the `shape_*` functions.
#'
#' @param data an optional `data.frame` or `list` containing the object structure variables.
#' @param formula a formula stating the object structure.
#'
#' @return `data.frame` or `list` with a `formula` attribute. 
#'
# #' @examples
get_dat <- function(data, formula) {
  
  .dimvars <- all.vars(formula)
  
  if(missing(data)) data <- NULL
  ext_vars <- setdiff(.dimvars, names(data))
  if(length(ext_vars) != 0) {
    ext_vars <- mget(ext_vars, environment(formula))
    data <- c(data, ext_vars)
  }
  
  attr(data, "formula") <- formula
  data
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Extract `shape_dat` formula
#' 
#' Extracts the formula stating the structure of the shape object in the data.
#'
#' @param x an object of class `shape_dat`.
#'
#' @return a formula
#' @export
#'
formula.shape_dat <- function(x) {
  attr(x, "formula")
} 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Extract `shape_dat` stored variable names
#' 
#' Extracts the variable names of the shape object stored in the data.
#'
#' @param x an object of class `shape_dat`.
#'
#' @return a named character vector.
#' @export
#'
shape_names <- function(x) {
  .shape_names <- attr(x, "shape_names")
  if(is.null(.shape_names)) 
    .shape_names <- list()
  .shape_names
} 

store_names <- function(.shape_names, ...) {
  nl <- list(...)
  nl.names <- names(nl)
  nl <- Filter(Negate(is.null), nl)
  nl_0 <- setdiff(nl.names, names(nl))
  if(length(nl)!=0)
    if(any(sapply(nl, Negate(is.character))))
      stop("Shape names have to be given as character.")
  .shape_names <- replace(.shape_names, names(nl), nl)
  if(length(nl_0) != 0) {
    nl_0 <- structure( as.list(nl_0), names = nl_0 )
    .shape_names <- replace(nl_0, names(.shape_names), .shape_names)
  }
  .shape_names
}

# shape_dat-class ---------------------------------------------------------------

#' @importFrom methods setOldClass
#' @exportClass shape_dat
#' @exportClass shape_frame
#' @exportClass shape_frame_default
#' @exportClass shape_frame_long
#' @exportClass shape_cube
#' @exportClass shape_cube_default
#' @exportClass shape_col
#' @exportClass shape_col_default
#' @exportClass shape_col_long

setOldClass("shape_dat")
setOldClass("shape_frame")
setOldClass("shape_frame_default")
setOldClass("shape_frame_long")
setOldClass("shape_cube")
setOldClass("shape_cube_default")
setOldClass("shape_col")
setOldClass("shape_col_default")
setOldClass("shape_col_long")

#' `shape_dat` S3 class and subclasses
#'
#' @description The `shape_dat` class is an abstract class for lists/data frames
#' containing a sample of 'shape' objects together with other variables carrying 
#' additional information. A formula attribute stores the information on the 
#' variables the shape is composed of. Its subclasses `shape_frame`, `shape_cube`,
#' `shape_col` identify storage schemes where the different variables belonging
#' to a shape are stored seperately in a `data.frame`, in a `tbl_cube` or as `list`
#' column in a data frame. Each of them are again parents to two subclasses 
#' reflecting the internal format of an shape object, which are `shape_frame_default`,
#' `shape_frame_long`, etc.
#'
#' @name shape_dat-class
#' @aliases shape_frame-class shape_frame_default-class shape_frame_long-class
#' @aliases shape_cube-class shape_cube_default-class 
#' @aliases shape_col-class shape_col_default-class shape_col_long-class
#' @seealso [shape_dat()]
#'  
#' @example tests/shape_dat-class_example.R
# #' @exampleDontrun R/shape_dat-class_donotrun.R  
NULL





# shape_dat constructors -------------------------------------------------------------


#' Build a `shape_dat` object
#' 
#' @description 
#' 
#' `shape_dat()` constructs a data frame of shapes of class `shape_dat`.
#' 
#' @param formula formula stating identifying '`shape`' variable and its internal 
#' strucure in the format `value^dim ~ arg | id`.
#' @param ... optional variables passed to `data.frame()` or `list()`. 
#'
#' @name shape_dat
#' @rdname shape_dat  
#' @example tests/shape_dat-class_example.R  
#'
new_shape_dat <- function(formula, ...) {
  if(missing(formula))
    stop("formula has to be provided.")
  template <- c(
    shape_frame_default = cbind(value1, value2) ~ arg | id,
    shape_frame_long = value^dim ~ arg | id,
    shape_cube_default = shape ~ value^dim ~ arg | id,
    shape_col_default = shape ~ cbind(value1, value2) ~ arg,
    shape_col_long = shape ~ value^dim ~ arg
  )
  
  which_format <- sapply(template, function(tp) {
    !inherits(try(match_template(formula, tp), silent = TRUE), "try-error")
  })
  
  if(sum(which_format) == 0)
    stop(paste("\nformula does not match any of the templates. Accepted formats are", print(template), collapse = "\n"))
  if(sum(which_format) != 1)
    stop("formula matches multiple templates. Sorry, this should not happen!")
  
  form <- names(which_format)[which_format]
  form <- paste0("new_", form)
  
  eval(parse(text = form))(formula, ...)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Build a `shape_frame_long` object
#' 
#' @description 
#' 
#' `shape_frame_long()` constructs a data frame of shapes in long format 
#' of class `shape_frame_long`.
#' 
#' @param formula formula stating identifying '`shape`' variable and its internal 
#' strucure in the format `value^dim ~ arg | id`. If variables 
#' in `formula` are provided in `...`, they are taken from there.
#' @param ... optional variables passed to `data.frame()`. 
#'  
#' @example tests/shape_dat-class_example.R 

new_shape_frame_long <- function(formula, ...) {
  
  template <- value^dim ~ arg | id
  .dimvars <- match_template(formula, template)
  x <- get_dat(data.frame(...), formula)
  
  class(x) <- unique(
    c("shape_frame_long", "shape_frame", "shape_dat", class(x)))
  
  if(validate) validate_shape_frame_long(x)
  return(x)
}

validate_shape_frame_long <- function(x) {
  .formula <- attr(x, "formula")
  .dimvars <- all.vars(.formula)
  assign(.dimvars, x[.dimvars]) 
  
  stopifnot(is.numeric(value))
  n <- length(value)
  m <- length(unique(dim))
  stopifnot(length(arg) == n & length(dim == n) )
  stopifnot(length(id) == n | length(id) == 1)
  if(m!=2) stop("So far, only planar shapes with two dimensions are provided.")
  
  if(all(table(paste(arg, dim, id)))!=1) 
    stop("Each arg-dim-id combination may only occur once.")
  if(any(tapply(dim, paste(arg, id), 
                function(x) length(x) != m | length(unique(x)) != m)))
    stop("For each point, there must be values provided in each dimension.")
  x
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Build a `shape_frame_default` object
#' 
#' @description 
#' 
#' `shape_frame_default()` constructs a data frame of shapes in long format 
#' of class `shape_frame_default`.
#' 
#' @param formula formula stating identifying '`shape`' variable and its internal 
#' strucure in the format `cbind(value1, value2) ~ arg | id`.
#' @param ... optional variables passed to `data.frame()`.
#'  
#' @example tests/shape_dat-class_example.R 
#' 

new_shape_frame_default <- function(formula, ...) {
  
  template <- cbind(value1, value2) ~ arg | id
  .dimvars <- match_template(formula, template)
  x <- get_dat(data.frame(...), formula)
  
  class(x) <- unique(
    c("shape_frame_default", "shape_frame", "shape_dat", class(x)))
  if(validate) validate_shape_frame_default(x)
  return(x)
}

validate_shape_frame_default <- function(x) {
  .formula <- attr(x, "formula")
  .dimvars <- all.vars(.formula)
  assign(.dimvars, x[.dimvars]) 
  
  stopifnot(is.numeric(value1) & is.numeric(value2))
  n <- length(value1)
  stopifnot(length(arg) == n & length(value2) == n )
  stopifnot(length(id) == n | length(id) == 1)
  if(all(table(paste(arg, id)))!=1) 
    stop("Each arg-id combination may only occur once.")
  x
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Build a `shape_cube_default` object
#' 
#' @description 
#' 
#' `shape_cube_default()` constructs a `tbl_cube` of shapes of class `shape_cube_default` with 
#' possibly an additional slot `info`.
#' 
#' @param formula formula stating identifying '`shape`' variable and its internal 
#' strucure in the format `shape ~ value^dim ~ arg | id`. These variables are used to construct
#' a `tbl_cube(measures = list(value), dimensions = list(arg, dim, id)` named `shape`.
#' @param ... optional variables passed to `tbl_cube()` as described above or to `list()`.
#'  
#' @example tests/shape_dat-class_example.R  
#' 
#' @import dplyr
new_shape_cube_default <- function(formula, ...) {
  
  template <- shape ~ value^dim ~ arg | id
  .shape <- match_template(formula, template)[1]
  
  inner_formula <- update(formula, paste(
    deparse(formula[[2]][[3]]), ~.)
  )
  
  .dimvars <- all.vars(inner_formula)
  .value <- .dimvars[1]
  .dimvars <- .dimvars[-1]
  
  x <- get_dat(list(...), inner_formula)
  x[[.shape]] <- tbl_cube(
    measures = x[.value],
    dimensions = x[.dimvars[c(2,1,3)]]
    )
  
  inner_vars <- which(names(x) %in% all.vars(inner_formula))
  x <- x[-inner_vars]
  
  class(x) <- c("shape_cube", "shape_cube_default", class(x))
  return(x)
}

validate_shape_cube_default <- function(x) {
  .formula <- formula(x)
  innervars <- all.vars(x)[-1]
  stopifnot(length(setdiff(
    innervars,
    names(x[[all.vars(x)[1]]]$dimensions))) == 0)
  x
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Build a `shape_col_long` object
#' 
#' @description 
#' 
#' `shape_col_long()` constructs a `data.frame` of class `shape_dat` with one 
#' list column containing shapes as a list of `data.frame`s 
#' (or objects convertible to such).
#' 
#' @param formula formula stating identifying '`shape`' variable and its internal 
#' strucure in the format `shape ~ value^dim ~ arg`.
#' @param ...
#'    
#' @example tests/shape_dat-class_example.R  
#'
new_shape_col_long <- function(formula, ...) {
  
  template <- shape ~ value^dim ~ arg
  .shape <- match_template(formula, template)[1]
  
  inner_formula <- update(formula, paste(
    deparse(formula[[2]][[3]]), ~.)
  )
  
  .dimvars <- all.vars(inner_formula)
  .value <- .dimvars[1]
  .dim <- .dimvars[2]
  .arg <- .dimvars[3]
  
  x <- get_dat(list(...), inner_formula)
  x[[.shpvar]] <- Map(
    function(value, arg) {
      structure( 
        data.frame(.value = value, .arg = arg, .dim = x[[.dim]]),
        names = c(.value, .arg, .dim))
    },
    x[[.value]],
    x[[.arg]]
  )
  
  inner_vars <- which(names(x) %in% all.vars(inner_formula))
  x <- do.call(as.data.frame, x[-inner_vars])
  
  class(x) <- unique(
    c("shape_col_long", "shape_col", class(x)))
  x
}

validate_shape_col_long <- function(x) {
  .formula <- formula(x)
  innervars <- all.vars(x)[-1]
  shp_names <- lapply(x[[all.vars(x)[1]]], names)
  stdf <- lapply(shp_names, setdiff, y = innervars)
  stopifnot(all(sapply(stdf, length) == 0))
  x
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Build a `shape_col_default` object
#' 
#' @description 
#' 
#' `shape_col_default()` constructs a `data.frame` of class `shape_dat` with one list column
#' containing shapes as a list of `data.frame`s (or objects convertible to such).
#' 
#' @param formula formula stating identifying '`shape`' variable and its internal 
#' strucure in the format `shape ~ cbind(value1, value2) ~ arg`.
#' @param ... 
#'  
#' @example tests/shape_dat-class_example.R  
#' 
new_shape_col_default <- function(x, formula, validate = TRUE) {
  
  template <- shape ~ cbind(value1, value2) ~ arg
  .shape <- match_template(formula, template)[1]
  
  inner_formula <- update(formula, paste(
    deparse(formula[[2]][[3]]), ~.)
  )
  
  .values <- all.vars(inner_formula[[2]])
  .arg <- all.vars(inner_formula[[3]])
  
  x <- get_dat(list(...), inner_formula)
  x[[.shpvar]] <- Map(
    function(arg, ...) {
      structure( 
        data.frame(..., arg, .dim = x[[.dim]]),
        names = c(.values, .arg, .dim))
    },
    x[[.arg]],
    x[.values]
  )
  
  inner_vars <- which(names(x) %in% all.vars(inner_formula))
  x <- do.call(as.data.frame, x[-inner_vars])
  
  class(x) <- unique(
    c("shape_col_default", "shape_col", class(x)))
  x
}


validate_shape_col_default <- function(x) {
  validate_shape_col_long(x)
}

# converters -------------------------------------------------------------------

#' Convert objects to `shape_dat` subclasses
#' 
#' @description 
#' `as_shape_dat` and friends are S3 generics converting other data formats to
#' `shape_dat` subclasses.
#' 
#' @param x an object that should be coerced.
#' @param ... other arguments passed to the individual methods.
#' 
#' @export
#' @name as_shape_dat
#' @rdname as_shape_dat
#' @example tests/shape_dat-class_example.R
#'
# #' @examples
as_shape_dat <- function(x, ...) {
  UseMethod("as_shape_dat")
} 

# default converters for a data.frame or list ----------------------------------

#' Convert a `data.frame` or `list` to a `shape_dat` object
#' 
#' @description 
#' 
#' `as_shape_dat.default()` converts a `data.frame` or `list` to an object of 
#' class `shape_dat`. The assigned subclass depends on the structure of the `formula` 
#' and `x` provided.  
#' 
#' @param x optional `data.frame` or `list` with relevant data.
#' @param formula formula defining the structure of the 'shape'. 
#' @param validate logical, should an additional validation of the constructed `shape_dat`
#' be performed (if any)? Defaults to TRUE. FALSE might speed up computation time.
#' 
#' @return usually just `x` with an extra attribute `'formula'`, plus eventually 
#' extra variables provided in the formula, but not in `x`.
#' @seealso `as_shape_frame_default`, `as_shape_cube_default`, `as_shape_col_default`, ...
#'
#' @name shape_dat
#' @rdname shape_dat  
#' @example tests/shape_dat-class_example.R  
#' @export
#'

as_shape_dat.default <- function(x, formula, validate = TRUE) {
  if(missing(formula))
    stop("formula has to be provided.")
  template <- c(
    shape_frame_default = cbind(value1, value2) ~ arg | id,
    shape_frame_long = value^dim ~ arg | id,
    shape_cube_default = shape ~ value^dim ~ arg | id,
    shape_col_default = shape ~ cbind(value1, value2) ~ arg,
    shape_col_long = shape ~ value^dim ~ arg
  )
  
  which_format <- sapply(template, function(tp) {
    !inherits(try(match_template(formula, tp), silent = TRUE), "try-error")
  })
  
  if(sum(which_format) == 0)
    stop(paste("\nformula does not match any of the templates. Accepted formats are", print(template), collapse = "\n"))
  if(sum(which_format) != 1)
    stop("formula matches multiple templates. Sorry, this should not happen!")
  
  form <- names(which_format)[which_format]
  form <- paste0("as_", form, ".default")
  
  eval(parse(text = form))(x, formula, validate)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Convert a `data.frame` or `list` to a `shape_frame_long`
#' 
#' @description 
#' 
#' `shape_frame_long()` constructs a data frame of shapes in long format 
#' of class `shape_frame_long`.
#' 
#' @param x optional data frame or list with relevant data.
#' @param formula formula stating identifying '`shape`' variable and its internal 
#' strucure in the format `value^dim ~ arg | id`.
#' @param validate logical, should an additional validation of the constructed `shape_dat`
#' be performed (if any)? Defaults to TRUE. FALSE might speed up computation time.
#'  
#' @export  
#' @example tests/shape_dat-class_example.R

as_shape_frame_long <- function(x, ...) {
  UseMethod("as_shape_frame_long")
} 

#' @export
#' @rdname as_shape_frame_long
as_shape_frame_long.default <- function(x, formula, validate = TRUE) {
  
  template <- value^dim ~ arg | id
  .allvars <- match_template(formula, template)
  x <- get_dat(x, formula)
  
  class(x) <- unique(
    c("shape_frame_long", "shape_frame", "shape_dat", class(x)))
  
  if(!validate) return(x)
  validate_shape_frame_long(x)
}

validate_shape_frame_long <- function(x) {
  .formula <- attr(x, "formula")
  template <- value^dim ~ arg | id
  .dimvars <- match_template(.formula, template)
  e <- environment()
  mapply(assign, names(.dimvars), x[.dimvars], MoreArgs = list(pos = e)) 
  
  stopifnot(is.numeric(value))
  n <- length(value)
  m <- length(unique(dim))
  stopifnot(length(arg) == n & length(dim == n) )
  stopifnot(length(id) == n | length(id) == 1)
  if(m!=2) stop("So far, only planar shapes with two dimensions are provided.")
  
  if(all(table(paste(arg, dim, id)))!=1) 
    stop("Each arg-dim-id combination may only occur once.")
  if(any(tapply(dim, paste(arg, id), 
                function(x) length(x) != m | length(unique(x)) != m)))
    stop("For each point, there must be values provided in each dimension.")
  x
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Convert a `data.frame` or `list` to a `shape_frame_default`
#' 
#' @description 
#' 
#' `shape_frame_default()` constructs a data frame of shapes in long format 
#' of class `shape_frame_default`.
#' 
#' @param x optional data frame or list with relevant data.
#' @param formula formula stating identifying '`shape`' variable and its internal 
#' strucure in the format `cbind(value1, value2) ~ arg | id`.
#' @param validate logical, should an additional validation of the constructed `shape_dat`
#' be performed (if any)? Defaults to TRUE. FALSE might speed up computation time.
#'  
#' @export
#' @name as_shape_frame_default
#' @example tests/shape_dat-class_example.R 
#' 
as_shape_frame_default <- function(x, ...) {
  UseMethod("as_shape_frame_default")
}

#' @export
#' @rdname as_shape_frame_default
as_shape_frame_default.default <- function(x, formula, validate = TRUE) {
  
  template <- cbind(value1, value2) ~ arg | id
  .allvars <- match_template(formula, template)
  x <- get_dat(x, formula)
  
  class(x) <- unique(
    c("shape_frame_default", "shape_frame", "shape_dat", class(x)))
  if(!validate) return(x)
  validate_shape_frame_default(x)
}

validate_shape_frame_default <- function(x) {
  .formula <- attr(x, "formula")
  template <- cbind(value1, value2) ~ arg | id
  .dimvars <- match_template(.formula, template)
  e <- environment()
  mapply(assign, names(.dimvars), x[.dimvars], MoreArgs = list(pos = e)) 
  
  stopifnot(is.numeric(value1) & is.numeric(value2))
  n <- length(value1)
  stopifnot(length(arg) == n & length(value2) == n )
  stopifnot(length(id) == n | length(id) == 1)
  if(all(table(paste(arg, id)))!=1) 
    stop("Each arg-id combination may only occur once.")
  x
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Convert a `data.frame` or `list` to a `shape_cube_default`
#' 
#' @description 
#' 
#' `shape_cube_default()` constructs a `tbl_cube` of shapes of class `shape_cube_default` with 
#' possibly an additional slot `info`.
#' 
#' @param x optional data frame or list with relevant data.
#' @param formula formula stating identifying '`shape`' variable and its internal 
#' strucure in the format `shape ~ value^dim ~ arg | id`.
#' @param validate logical, should an additional validation of the constructed `shape_dat`
#' be performed (if any)? Defaults to TRUE. FALSE might speed up computation time.
#'  
#' @export
#' @name as_shape_cube_default
#' @example tests/shape_dat-class_example.R  
#' 
#' @import dplyr

as_shape_cube_default <- function(x, ...) {
  UseMethod("as_shape_cube_default")
}

#' @export
#' @rdname as_shape_cube_default
as_shape_cube <- function(x, ...) {
  UseMethod("as_shape_cube_default")
}


#' @export
#' @rdname as_shape_cube_default
as_shape_cube_default.default <- function(x, formula, validate = TRUE) {
  
  template <- shape ~ value^dim ~ arg | id
  .allvars <- match_template(formula, template)
  shp <- .allvars[1]
  shpformula <- as.formula(paste(shp, "~0"), env = environment(formula))
  x <- get_dat(x, shpformula)
  stopifnot(inherits(x[[shp]], "tbl_cube"))
  
  class(x) <- c("shape_cube", "shape_cube_default", "shape_dat", class(x))
  if(!validate) return(x)
  validate_shape_cube_default(x)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Convert a `data.frame` or `list` to a `shape_col_long`
#' 
#' @description 
#' 
#' `shape_col_long()` constructs a `data.frame` of class `shape_dat` with one 
#' list column containing shapes as a list of `data.frame`s 
#' (or objects convertible to such).
#' 
#' @param x optional data frame or list with relevant data.
#' @param formula formula stating identifying '`shape`' variable and its internal 
#' strucure in the format `shape ~ value^dim ~ arg`.
#' @param validate logical, should an additional validation of the constructed `shape_dat`
#' be performed (if any)? Defaults to TRUE. FALSE might speed up computation time.
#'  
#' @export
#' @example tests/shape_dat-class_example.R  
#'
as_shape_col_long <- function(x, ...) {
  UseMethod("as_shape_col_long")
}

#' @export
#' @rdname as_shape_col_long
as_shape_col_long.default <- function(x, formula, validate = TRUE) {
  
  template <- shape ~ value^dim ~ arg
  .allvars <- match_template(formula, template)
  shp <- .allvars[1]
  shpformula <- as.formula(paste(shp, "~0"), env = environment(formula))
  x <- get_dat(x, shpformula)
  stopifnot(is.list(x[[shp]]))
  
  attr(x, "formula") <- formula
  class(x) <- unique(
    c("shape_col_long", "shape_col", "shape_dat", class(x)))
  x
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Convert a `data.frame` or `list` to a `shape_col_default`
#' 
#' @description 
#' 
#' `shape_col_default()` constructs a `data.frame` of class `shape_dat` with one list column
#' containing shapes as a list of `data.frame`s (or objects convertible to such).
#' 
#' @param x optional data frame or list with relevant data.
#' @param formula formula stating identifying '`shape`' variable and its internal 
#' strucure in the format `shape ~ cbind(value1, value2) ~ arg`.
#' @param validate logical, should an additional validation of the constructed `shape_dat`
#' be performed (if any)? Defaults to TRUE. FALSE might speed up computation time.
#'  
#' @export
#' @example tests/shape_dat-class_example.R  
#' 

as_shape_col_default <- function(x, ...) {
  UseMethod("as_shape_col_default")
}

#' @export
#' @rdname as_shape_col_default
as_shape_col_default.default <- function(x, formula, validate = TRUE) {
  
  template <- shape ~ cbind(value1, value2) ~ arg
  .allvars <- match_template(formula, template)
  shp <- .allvars[1]
  shpformula <- as.formula(paste(shp, "~0"), env = environment(formula))
  x <- get_dat(x, shpformula)
  stopifnot(is.list(x[[shp]]))
  
  attr(x, "formula") <- formula
  class(x) <- unique(
    c("shape_col_default", "shape_col", "shape_dat", class(x)))
  x
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# mutual shape_dat conversions -------------------------------------------------

#' @export
#' @rdname as_shape_col_default
as_shape_col <- function(x, ...) {
  UseMethod("as_shape_col")
}

#' @export
#' @rdname as_shape_col_default
as_shape_col.shape_frame <- function(x, shape.name = NULL, ...) {
  .shape_names <- store_names( shape_names(x), shape = shape.name )
  shape <- .shape_names$shape
  .class <- class(x)
  .formula <- formula(x)
  .allvars <- all.vars(.formula)
  shp <- x[.allvars[-length(.allvars)]]
  x <- x[setdiff(names(x), names(shp))]
  idvar <- .allvars[length(.allvars)]
  # make sure order is preserved
  x[[idvar]] <- ordered(x[[idvar]], levels = unique(x[[idvar]]))
  shp <- split(shp, x[[idvar]])
  x <- x %>% group_by_at(idvar) %>% summarise_all( ~ head(., 1) ) %>% ungroup
  x[[idvar]] <- NULL
  x[[shape]] <- lapply(shp, function(x) 
    structure(x, class = class(x)[-grep("shape_", class(x))]) )
  .formula[[3]] <- .formula[[3]][[2]]
  attr(x, "shape_names") <- store_names( .shape_names, id = idvar )
  .formula[[2]] <- str2lang(paste(shape, "~", deparse(.formula[[2]])))
  attr(x, "formula") <- .formula
  structure(x, 
            class = sub("shape_frame", "shape_col", .class))
}

#' @export
#' @rdname as_shape_frame_default
as_shape_frame <- function(x, ...) {
  UseMethod("as_shape_frame")
}

#' @export
#' @rdname as_shape_frame_default
as_shape_frame.shape_col <- function(x, id.name = NULL, ...) {
  .shape_names <- store_names( shape_names(x), id = id.name )
  id <- .shape_names$id
  .class <- class(x)
  .class <- .class[grep("^shape_", .class)]
  .class <- sub("shape_col", "shape_frame", .class)
  .formula <- formula(x)
  .allvars <- all.vars(.formula)
  .shp <- .allvars[1]
  shp <- x[[.shp]]
  if(is.null(names(shp))) names(shp) <- seq_along(shp)
  shp_lens <- sapply(shp, nrow)
  shp <- bind_rows(shp, .id = id)
  x <- as.data.frame( x[-which(names(x) == .shp)] )
  x <- x[ rep(seq_len(nrow(x)), shp_lens), , drop = FALSE]
  x <- cbind(shp, x)
  .formula[[2]] <- .formula[[2]][[3]]
  attr(x, "shape_names") <- store_names( .shape_names, shape = .shp )
  .formula[[3]] <- str2lang(paste(deparse(.formula[[3]]), "|", id))
  attr(x, "formula") <- .formula
  structure(x, 
            class = c(.class, class(x)))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @export
#' @rdname as_shape_frame_default
#' @importFrom tidyr pivot_wider
as_shape_frame_default.shape_frame_long <- function(x, ...) {
  .formula <- formula(x)
  .valuevar <- all.vars(.formula)[1]
  .dimvar <- all.vars(.formula)[2]
  value.vars <- unique(x[[.dimvar]])
  numdim <- !all( is.na( suppressWarnings( as.numeric( as.character( value.vars) ) ) ) )
  if(numdim) {
    x[[.dimvar]] <- paste0(.valuevar, x[[.dimvar]])
    value.vars <- paste0(.valuevar, value.vars)
  }
  .shape_names <- store_names( shape_names(x), 
                               value = .valuevar,
                               dim = .dimvar )
  # .othervars <- setdiff(names(x), all.vars(.formula))
  attr(x, "formula") <- attr(x, "shape_names") <- NULL
  x <- pivot_wider(x, 
               # id_cols = c(all.vars(.formula[[3]]), .othervars),
               values_from = .valuevar, 
               names_from = .dimvar)
  .formula[[2]] <- str2lang(
    paste("cbind(", paste(value.vars, collapse = ","), ")"))
  x <- as_shape_dat(x, formula = .formula)
  attr(x, "shape_names") <- .shape_names
  x
}

#' @export
#' @rdname as_shape_frame_long
as_shape_frame_long.shape_frame_default <- function(x, 
                                                    value.name = NULL, 
                                                    dim.name = NULL, ...) {
  .shape_names <- store_names( shape_names(x), 
                                         value = value.name,
                                         dim = dim.name )
  value <- .shape_names$value
  dim <- .shape_names$dim
  .formula <- formula(x)
  .othervars <- setdiff(names(x), all.vars(.formula))
  x <- reshape(x, 
               varying = all.vars(.formula[[2]]), 
               idvar = all.vars(.formula[[3]]), 
               v.names = value, 
               direction = "long", 
               timevar = dim)
  .formula[[2]] <- str2lang(paste(value, dim, sep = "^"))
  attr(x, "reshapeLong") <- NULL
  attr(x, "shape_names") <- .shape_names
  attr(x, "formula") <- .formula
  structure(x, class = sub("shape_frame_default", "shape_frame_long", class(x)))
}
#' @export
#' @rdname as_shape_frame_default
as_shape_frame_default.shape_dat <- function(x, ...) {
  as_shape_frame_default(as_shape_frame(x), ...)
}
#' @export
#' @rdname as_shape_frame_long
as_shape_frame_long.shape_dat <- function(x, ...) {
  as_shape_frame_long(as_shape_frame(x), ...)
}

#' @export
#' @rdname as_shape_col_default
as_shape_col_default.shape_dat <- function(x, ...) {
  as_shape_col(as_shape_frame_default(x), ...)
}
#' @export
#' @rdname as_shape_col_long
as_shape_col_long.shape_dat <- function(x, ...) {
  as_shape_col(as_shape_frame_long(x, ...), ...)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @export
#' @rdname as_shape_frame_long
as_shape_frame_long.shape_cube <- function(x, ...) {
  .formula <- formula(x)
  .allvars <- all.vars(.formula)
  .shp <- .allvars[1]
  .id <- .allvars[length(.allvars)]
  shp <- as.data.frame(x[[.shp]])
  # order preserving ids
  idvar <- as.numeric(ordered(shp[[.id]], levels = unique(shp[[.id]])))
  x <- x[ idvar, -which(names(x) == .shp)]
  x <- cbind(shp, x)
  .formula[[2]] <- .formula[[2]][[3]]
  attr(x, "formula") <- .formula
  attr(x, "shape_names") <- store_names( shape_names(x), shape = .shp )
  structure(x, 
            class = c("shape_frame_long", "shape_frame", 
                      setdiff(class(x), c("shape_cube_default", "shape_cube"))))
}

#' @export
#' @rdname as_shape_cube_default
as_shape_cube_default.shape_frame_long <- function(x, shape.name = NULL, ...) {
  .shape_names <- store_names( shape_names(x), shape = shape.name )
  shape <- .shape_names$shape
  .formula <- formula(x)
  .allvars <- all.vars(.formula)
  shp <- x[.allvars]
  shp <- as.tbl_cube(shp, 
                     dim_names = .allvars[c(3,2,4)],
                     met_name = .allvars[1])
  x <- x[setdiff(names(x), .allvars[-length(.allvars)])]
  idvar <- .allvars[length(.allvars)]
  x <- x %>% group_by_at(idvar) %>% summarise_all( ~ head(., 1) ) %>% ungroup
  x <- as.list(x)
  x[[shape]] <- shp
  x[[idvar]] <- NULL
  .formula[[2]] <- str2lang(paste(shape, "~", deparse(.formula[[2]])))
  attr(x, "formula") <- .formula
  structure(x, 
            class = c("shape_cube_default", "shape_cube", 
                      setdiff(class(x), c("shape_frame_long", "shape_frame"))))
}

#' @export
#' @rdname as_shape_cube_default
as_shape_cube_default.shape_dat <- 
  as_shape_cube_default.shape_dat <- 
    function(x, ...) {
      x <- as_shape_cube(as_shape_frame_long(x), ...)
    }

# non-exportet formats and converters ------------------------------------------

setOldClass("shape_frame_complex")

#' `shape_frame_complex` S3 class
#'
#' @description Internal `shape_dat` subclass for planar shapes corresponding to 
#' `shape_frame_default`, only that instead of having two colums `cbind(value1, value2)`
#' with numeric values, there is one complex column `value`. The formula is, 
#' accordingly, `value ~ arg | id`.
#'
#' @name shape_frame_complex-class
#'  
NULL

#' Convert objects to `shape_frame_complex`
#' 
#' @description 
#' Internal function converting `shape_dat` objects to `shape_frame_complex`.
#' As I didn't need it, yet, there's no constructers, validators, or converters
#' from other formats or form `shape_frame_complex` to something else implemented
#' so far.
#' 
#' @param x an object that should be coerced.
#' @param ... other arguments passed to the individual methods.
#' 
#' @seealso `shape_frame_complex-class`
#' 
#' @name as_shape_frame_complex
as_shape_frame_complex <- function(x, ...) 
  UseMethod("as_shape_frame_complex")

as_shape_frame_complex.shape_frame_default <- function(x, value.name = NULL, ...) {
  .formula <- formula(x)
  .shape_names <- store_names( shape_names(x), 
                              value = value.name )
  .class <- class(x)
  value <- .shape_names$value
  value_dims <- all.vars(.formula[[2]])
  x[[value]] <- complex(real = x[[value_dims[1]]], imag = x[[value_dims[2]]])
  value_dims <- setdiff(value_dims, value)
  x <- x[- which(names(x) %in% value_dims)]
  .formula[[2]] <- as.name(value)
  attr(x, "formula") <- .formula
  attr(x, "shape_names") <- .shape_names
  class(x) <- sub("shape_frame_default", "shape_frame_complex", .class)
  x
}


as_shape_frame_complex.shape_dat <- function(x, ...) {
  as_shape_frame_complex.shape_frame_default( as_shape_frame_default(x, ...) )
}

as_shape_frame_long.shape_frame_complex <- function(x, dim.name = NULL, ...) {
  .formula <- formula(x)
  .shape_names <- store_names(shape_names(x), dim = dim.name)
  .value <- all.vars(.formula)[1]
  .dim <- .shape_names[["dim"]]
  .class <- class(x)
  attr(x, "formula") <- attr(x, "shape_names") <- NULL
  x <- list(real = x, imaginary = x)
  x[[1]][[.value]] <- Re(x[[1]][[.value]])
  x[[2]][[.value]] <- Im(x[[2]][[.value]])
  x <- bind_rows(x, .id = .dim)
  attr(x, "formula") <- as.formula(
    paste0(deparse(.formula[[2]]), "^", .dim, "~", deparse(.formula[[3]])),
    env = environment(.formula))
  attr(x, "shape_names") <- .shape_names
  class(x) <- sub("shape_frame_complex", "shape_frame_long", .class)
  x
}

# trivial converters -----------------------------------------------------------

#' @export
as_shape_dat.shape_dat <-
 as_shape_col.shape_col <- 
  as_shape_col_default.shape_col_default <- 
  as_shape_col_long.shape_col_long <- 
  as_shape_frame.shape_frame <- 
  as_shape_frame_default.shape_frame_default <- 
  as_shape_frame_long.shape_frame_long <- 
  as_shape_cube.shape_cube <- 
  as_shape_cube_default.shape_cube_default <- function(x, ...) x

# `[` operator -----------------------------------------------------------------

#' @export
`[.shape_dat` <- function(x, ...) {
  .formula <- formula(x)
  .shape_names <- shape_names(x)
  x <- NextMethod()
  attr(x, "formula") <- .formula
  attr(x, "shape_names") <- .shape_names
  x
}



