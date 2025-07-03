#' Plot Fit Method for unidata Objects
#'
#' This method visualizes the `unidata` and `multidata` objects by plotting both the 
#' original time series and the fitted values. The `unidata` (`multidata`) object is 
#' converted into a `zoo` object for plotting, and the fitted values are 
#' added as a secondary line with a different style and color.
#'
#' @usage
#' plot_fit(x,...)
#' 
#' @param x An object of class `unidata` or `multidata`. This object must contain the 
#'   time series in the `series` slot, the associated time points in the 
#'   `times` slot, and the fitted values in the `fitted_values` slot.
#' @param ... Additional arguments passed to the `plot.zoo` function for 
#'   customizing the plot (e.g., titles, colors, etc.).
#'
#' @importFrom zoo plot.zoo zoo
#' @export
plot_fit<- S7::new_generic("plot_fit", "x")
S7::method(plot_fit, unidata) <- function(x, ...) {
  zoo::plot.zoo(zoo::zoo(x@series, x@times), ...)
  lines(zoo::zoo(x@fitted_values, x@times), col = 2, lty = 2)
}
S7::method(plot_fit, multidata) <- function(x, ...) {
  zoo::plot.zoo(zoo::zoo(x@series[,1], x@times), ...)
  lines(zoo::zoo(x@fitted_values[,1], x@times), col = 2, lty = 2)
  readline("Press [Enter] to show the next plot...")
  zoo::plot.zoo(zoo::zoo(x@series[,2], x@times), ...)
  lines(zoo::zoo(x@fitted_values[,2], x@times), col = 2, lty = 2)
}