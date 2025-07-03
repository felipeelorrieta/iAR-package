#' Plot Method for unidata Objects
#'
#' This method allows visualizing `unidata` and `multidata` objects using the `plot.zoo` 
#' function from the `zoo` package. It converts the `unidata` (`multidata`) object into 
#' a `zoo` object and then applies the plotting function to the time series 
#' contained in the object.
#' @name plot
#' 
#' @param x An object of class `unidata` or `multidata`. This object must contain the time 
#'   series in the `series` slot and the associated time points in the 
#'   `times` slot.
#' @param ... Additional arguments passed to the `plot.zoo` function for 
#'   customizing the plot (e.g., titles, colors, etc.).
#'
#' @importFrom zoo plot.zoo zoo
#' @export
S7::method(plot, unidata) <- function(x, ...) {
  zoo::plot.zoo(zoo::zoo(x@series, x@times), ...)
}
S7::method(plot, multidata) <- function(x, ...) {
  zoo::plot.zoo(zoo::zoo(x@series[,1], x@times), ...)
  readline("Press [Enter] to show the next plot...")
  zoo::plot.zoo(zoo::zoo(x@series[,2], x@times), ...)
}