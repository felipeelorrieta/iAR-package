.onLoad <- function(...) {
  S7::methods_register()
}

# Format one coefficient table value; returns "  NA" when SE is unavailable.
.fmt_se <- function(x) if (is.null(x) || !is.finite(x)) "  NA" else format(round(x, 2), nsmall = 2)

# Clamp an iAR estimated coefficient to (0, 1) and issue warnings if it
# hits the boundary or is suspiciously close to it.
.clamp_iar_coef <- function(phi, eps = 1e-4, near_boundary = 0.995) {
  if (!is.finite(phi) || phi <= 0 || phi >= 1) {
    phi_clamped <- max(eps, min(1 - eps, phi))
    warning(paste0(
      "Estimated coefficient (", round(phi, 6), ") is at or beyond the boundary of (0, 1). ",
      "It has been clamped to ", phi_clamped, ". ",
      "Results may be unreliable. If negative autocorrelation is expected, consider using CiAR."
    ), call. = FALSE)
    return(phi_clamped)
  }
  if (phi > near_boundary) {
    warning(paste0(
      "Estimated coefficient (", round(phi, 4), ") is near the upper boundary of (0, 1). ",
      "The model may be near non-stationary; results should be interpreted with caution."
    ), call. = FALSE)
  } else if (phi < (1 - near_boundary)) {
    warning(paste0(
      "Estimated coefficient (", round(phi, 4), ") is near zero. ",
      "This may indicate weak autocorrelation or convergence issues."
    ), call. = FALSE)
  }
  phi
}
