#' Multidata Class
#'
#' The `multidata` class is an S7 class designed to represent multidimensional time series models, including the main time series and additional series (e.g., error standard deviations or related variables).
#'
#' @param times A numeric vector representing the time points of the time series.
#' @param series A numeric vector or matrix representing the main time series.
#' @param series_esd A numeric vector or matrix representing the additional series, such as error standard deviations or other related data.
#' @param series_names An optional character vector representing the name of the series. 
#'
#' @section Validation Rules:
#' - `@times` must be a numeric vector and strictly increasing.
#' - `@series` must be a numeric matrix or empty. Its number of rows must match the length of `@times`.
#' - `@series_esd` must be a numeric matrix or empty. If provided and not empty, its dimensions must match those of `@series`.
#' - If `@series_names` is provided, it must be a character vector of the same length as the number of columns of `@series`, and contain unique values.
#'
#' @examples
#' # Create a multidata object
#' multidata_instance <- multidata(
#'   times = c(1, 2, 3, 4),
#'   series = matrix(c(10, 20, 15, 25,
#'                     12, 18, 17, 23), 
#'                   nrow = 4, ncol = 2),
#'   series_esd = matrix(c(1, 1.5, 1.2, 1.8,
#'                         0.9, 1.3, 1.1, 1.7),
#'                       nrow = 4, ncol = 2)
#' )
#' @seealso [BiAR]
#' @export
multidata <- S7::new_class(
  "multidata",
  package = "iAR",
  properties = list(
    times = S7::class_numeric,
    series = S7::class_numeric,
    series_esd = S7::class_numeric,
    series_names = S7::new_property(S7::class_character, default = character(0)),
    series_esd_names = S7::new_property(
      class = S7::class_character,
      getter = function(self) {
        if (length(self@series_names) > 0) {
          paste0(self@series_names, "_esd")
        } else {
          character(0)
        }
      }
    )
  ), validator = function(self) {
    if (!(is.atomic(self@times) && is.null(dim(self@times)))) {
      return("@times must be a numeric vector without dimensions")
    }
    if (!(is.matrix(self@series) && is.numeric(self@series)) && length(self@series) != 0) {
      return("@series must be a numeric matrix or be empty")
    }
    if (!(is.matrix(self@series_esd) && is.numeric(self@series_esd)) && length(self@series_esd) != 0) {
      return("@series_esd must be a numeric matrix or empty")
    }
    if (length(self@series) > 0) {
    if (nrow(self@series) != length(self@times)) {
      return("Number of rows in @series must match length of @times")
    }
    }
    if (length(self@series_esd) > 2) {
      if (!all(dim(self@series) == dim(self@series_esd))) {
        return("Dimensions of @series_esd must match those of @series")
      }
    }
    # Validar series_names
    if (ncol(self@series) > 0 && length(self@series_names) > 0) {
      if (length(self@series_names) != ncol(self@series)) {
        return("Length of @series_names must match number of columns in @series")
      }
      if (anyDuplicated(self@series_names)) {
        return("@series_names must be unique")
      }
    }
    if (is.unsorted(self@times, strictly = TRUE)) {
      return("@times must be strictly increasing")
    }
    NULL
  }
)