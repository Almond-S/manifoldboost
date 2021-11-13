#' Cell outlines simulated from CPM 
#'
#' 2D cell outlines extracted from cells simulated by 
#' the Cellular Potts Model (CPM) proposed by Thueroff et al. 2019, a stochastic 
#' biophysical model, for calibration to real-life cells (Schaffer 2021).
#' 
#' @details The datasets contain the following covariate vectors:
#' \itemize{
#' \item \code{a} substrate adhesion (CPM parameter, numeric)
#' \item \code{r} signaling radius (CPM parameter, numeric)
#' \item \code{b} bulk stiffness (CPM parameter, numeric)
#' \item \code{m} membrane stiffness (CPM parameter, numeric)
#' }
#' For the irregular dataset \code{cells}:
#' \item \code{response} A tibble with the cell outline data:
#'  \itemize{
#'   \item \code{id} cell ID (factor, ordered according to covariate vectors)
#'   \item \code{arg} parameterization of the curves (numeric, between 0 and 1)
#'   \item \code{dim} dimension indicator (two-leveled factor)
#'   \item \code{value} x- and y-values of observed points of the outline (numeric)
#'  }
#' }
#' For the regular dataset \code{cellr}:
#' \item \code{response} A tibble with the cell outline data:
#'  \itemize{
#'   \item \code{id} cell ID (factor, ordered according to covariate vectors)
#'   \item \code{arg} parameterization of the curves (numeric, between 0 and 1)
#'   \item \code{dim} dimension indicator (two-leveled factor)
#'   \item \code{value} x- and y-values of observed points of the outline (numeric)
#'  }
#' }
#' 
#' @docType data
#' @keywords datasets
#' @usage data(cells, package = "manifoldboost")
#' @source Thankfully provided by Sophia Schaffer.
#' @references
#' Schaffer, S. A. (2021).
#' Cytoskeletal  dynamics  in  confined  cell  migration:  experiment  andmodelling.
#' PhD thesis, LMU Munich.  DOI: 10.5282/edoc.28480.
#'
#' Thueroff, F., A. Goychuk, M. Reiter, and E. Frey (2019).
#' Bridging the gap betweensingle-cell migration and collective dynamics.
#' eLife  8, e46842.
#' 
"cells"

#' @docType data
#' @usage data(cellr, package = "manifoldboost")
#' @rdname cells
"cellr"