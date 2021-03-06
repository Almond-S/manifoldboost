% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shape_dat-class.R
\name{shape_frame_complex-class}
\alias{shape_frame_complex-class}
\title{\code{shape_frame_complex} S3 class}
\description{
Internal \code{shape_dat} subclass for planar shapes corresponding to
\code{shape_frame_default}, only that instead of having two colums \code{cbind(value1, value2)}
with numeric values, there is one complex column \code{value}. The formula is,
accordingly, \code{value ~ arg | id}.
}
