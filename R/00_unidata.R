#' Unidata Class
#'
#' The `unidata` class is an S7 class designed to represent univariate irregularly observed time series models
#' with associated times, values, and optional error standard deviations.
#'
#' @param times A numeric vector representing the time points.
#' @param series A numeric vector representing the values of the time series.
#' @param series_esd A numeric vector representing the error standard deviations of the time series.
#' @param series_names An optional character vector of length 1 representing the name of the series. 
#'
#' @section Validation Rules:
#' - `@times`, `@series`, and `@series_esd` must be numeric vectors.
#' - `@times` must not contain `NA` values and must be strictly increasing.
#' - The length of `@series` must match the length of `@times`.
#' - The length of `@series_esd` must be 0, 1, or equal to the length of `@series`.
#' - `NA` values in `@series` must correspond exactly (positionally) to `NA` values in `@series_esd`.
#' - `@series_names`, if provided, must be a character vector of length 1.
#'
#' @examples
#' # Create a unidata object
#' unidata_instance <- unidata(
#'   times = c(1, 2, 3, 4),
#'   series = c(10, 20, 15, 25),
#'   series_esd = c(1, 1.5, 1.2, 1.8),
#'   series_names = "my_series")
#'
#' @seealso [iAR], [CiAR]
#' @export
unidata <- S7::new_class(
    "unidata",
    package = "iAR",
    properties = list(
      times = S7::class_numeric,
      series = S7::class_numeric,
      series_esd = S7::class_numeric,
      series_names = S7::new_property(S7::class_character, default = character(0)),
      series_esd_names = S7::new_property(
        class = S7::class_character,
        getter = function(self) {
          if (length(self@series_names) == 1) {
            return(paste0(self@series_names, "_esd"))
          } else {
            return(character(0))
          }
        }
      )
    ),
    validator = function(self) {
      if (!(is.atomic(self@times) && is.null(dim(self@times))) ||
          !is.vector(self@series) || !is.vector(self@series_esd)) {
        return("@times, @series, and @series_esd must be numeric vectors without dimensions")
      }
      if (anyNA(self@times)) {
        return("NA values are not allowed in @times")
      }
      if (length(self@series_esd) == length(self@series) & sum(na.omit(self@series_esd))>0) {
        mismatch <- is.na(self@series) != is.na(self@series_esd)
        if (any(mismatch)) {
          return("NA values in @series must correspond to NA values in @series_esd")
        }
      }
      if (length(self@series) > 0 && length(self@times) != length(self@series)) {
        return("The length of @times must match the length of @series")
      }
      if (!(length(self@series_esd) %in% c(0, 1, length(self@series)))) {
        return("The length of @series_esd must be 0, 1, or match @series")
      }
      if (is.unsorted(self@times, strictly = TRUE)) {
        return("@times must be strictly increasing")
      }
      # Validar series_names
      if (length(self@series_names) > 0) {
        if (length(self@series_names) != 1) return("Length of @series_names must be 1")
        # if (anyDuplicated(self@series_names)) return("@series_names must be unique")
      }
      NULL
    }
)