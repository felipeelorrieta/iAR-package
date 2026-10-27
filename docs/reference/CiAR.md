# \`CiAR\` Class

Represents a complex irregular autoregressive (CiAR) time series model.
This class extends the \`unidata\` class and provides additional
properties for modeling, forecasting, and interpolating irregularly
observed time series data with both negative and positive
autocorrelation.

## Usage

``` r
CiAR(
  times = integer(0),
  series = integer(0),
  series_esd = integer(0),
  series_names = character(0),
  fitted_values = integer(0),
  kalmanlik = integer(0),
  coef = c(0.9, 0),
  tAhead = 1,
  forecast = integer(0),
  interpolated_values = integer(0),
  interpolated_times = integer(0),
  interpolated_series = integer(0),
  zero_mean = TRUE,
  standardized = TRUE,
  hessian = FALSE,
  summary = list()
)
```

## Arguments

- times:

  A numeric vector representing the time points.

- series:

  A complex vector representing the values of the time series.

- series_esd:

  A numeric vector representing the error standard deviations of the
  time series.

- series_names:

  An optional character vector of length 1 representing the name of the
  series.

- fitted_values:

  A numeric vector containing the fitted values from the model.

- kalmanlik:

  A numeric value representing the Kalman likelihood of the model.

- coef:

  A numeric vector of length 2 containing the real (phiR) and imaginary
  (phiI) parts of the autoregressive coefficient. Each value must be in
  \[-1, 1\] and the modulus \`sqrt(phiR^2 + phiI^2)\` must be strictly
  less than 1 (stationarity condition). Defaults to \`c(0.9, 0)\`.
  Unlike \[iAR\], which uses a scalar, both \`CiAR\` and \[BiAR\] use a
  length-2 vector.

- tAhead:

  A numeric value specifying the forecast horizon (default: 1).

- forecast:

  A numeric vector containing the forecasted values.

- interpolated_values:

  A numeric vector containing the interpolated values.

- interpolated_times:

  A numeric vector containing the times of the interpolated data points.

- interpolated_series:

  A numeric vector containing the interpolated series.

- zero_mean:

  A logical value indicating if the model assumes a zero-mean process
  (default: TRUE).

- standardized:

  A logical value indicating if the model assumes a standardized process
  (default: TRUE).

- hessian:

  A logical value indicating whether the Hessian matrix is computed
  during estimation (default: FALSE).

- summary:

  A list containing the summary of the model fit, including diagnostics
  and statistical results.

## Details

The \`CiAR\` class is designed to handle irregularly observed time
series data with either negative or positive autocorrelation using an
autoregressive approach. It extends the \`unidata\` class to include
functionalities specific to the \`CiAR\` model.

Key features of the \`CiAR\` class include: - Support for irregularly
observed time series data with negative or positive autocorrelation. -
Forecasting and interpolation functionalities for irregular time
points. - Configurable assumptions of zero-mean and standardized
processes.

## Validation

\- Inherits all validation rules from the \`unidata\` class: -
\`@times\`, \`@series\`, and \`@series_esd\` must be numeric vectors. -
\`@times\` must not contain \`NA\` values and must be strictly
increasing. - The length of \`@series\` must match the length of
\`@times\`. - The length of \`@series_esd\` must be 0, 1, or equal to
the length of \`@series\`. - \`NA\` values in \`@series\` must
correspond exactly (positionally) to \`NA\` values in \`@series_esd\`. -
\`@series_names\`, if provided, must be a character vector of length 1.

\- \`@coef\` must be a numeric vector of length 2 with no dimensions. -
Each value in \`@coef\` must be in the interval \[-1, 1\]. - The modulus
\`sqrt(phiR^2 + phiI^2)\` must be strictly less than 1 (stationarity
condition). - \`@tAhead\` must be a strictly positive numeric scalar.

## References

Elorrieta, F, Eyheramendy, S, Palma, W (2019). “Discrete-time
autoregressive model for unequally spaced time-series observations.”
*A&A*, **627**, A120.
[doi:10.1051/0004-6361/201935560](https://doi.org/10.1051/0004-6361/201935560)
.

## Examples

``` r
times=gentime(n=200, distribution = "expmixture",
             lambda1 = 130, lambda2 = 6.5, p1 = 0.15, p2 = 0.85)@times
my_CiAR <- CiAR(times = times,coef = c(0.9, 0))

# Access properties
my_CiAR@coef
#> [1] 0.9 0.0
```
