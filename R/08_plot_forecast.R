#' Plot Forecast Method for unidata Objects
#'
#' This method visualizes the `unidata` and `multidata` objects by plotting both the 
#' original time series and the forecasted values. The `unidata` (`multidata`) object 
#' is converted into a `zoo` object for plotting, and the forecasted 
#' values are added as points with a different color.
#'
#' @usage
#' plot_forecast(x,...)
#' 
#' @param x An object of class `unidata` or `multidata`. This object must contain the 
#'   time series in the `series` slot, the associated time points in the 
#'   `times` slot, the forecasted values in the `forecast` slot, and the 
#'   forecast horizon in the `tAhead` slot.
#' @param ... Additional arguments passed to the `plot.zoo` function for 
#'   customizing the plot (e.g., titles, colors, etc.).
#'
#' @importFrom zoo plot.zoo zoo
#' @export
plot_forecast <- S7::new_generic("plot_forecast", "x")
S7::method(plot_forecast, unidata) <- function(x, ...) {
  zoo::plot.zoo(zoo::zoo(x@series,x@times), xlim=c(min(x@times),max(x@times)+x@tAhead), ...)
  points(zoo::zoo(x@forecast,max(x@times)+x@tAhead), col = 2)
}
S7::method(plot_forecast, multidata) <- function(x, ...) {
  zoo::plot.zoo(zoo::zoo(x@series[,1],x@times), xlim=c(min(x@times),max(x@times)+x@tAhead), ...)
  points(zoo::zoo(x@forecast[,1],max(x@times)+x@tAhead), col = 2)
  readline("Press [Enter] to show the next plot...")
  zoo::plot.zoo(zoo::zoo(x@series[,1],x@times), xlim=c(min(x@times),max(x@times)+x@tAhead), ...)
  points(zoo::zoo(x@forecast[,1],max(x@times)+x@tAhead), col = 2)
}

