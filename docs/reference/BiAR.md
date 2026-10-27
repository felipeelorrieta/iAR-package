# \`BiAR\` Class

Represents a bivariate irregular autoregressive (BiAR) time series
model. This class extends the \`multidata\` class and provides
additional properties for modeling, forecasting, and interpolation of
bivariate time series data.

## Usage

``` r
BiAR(
  times = integer(0),
  series = integer(0),
  series_esd = integer(0),
  series_names = character(0),
  fitted_values = integer(0),
  loglik = integer(0),
  kalmanlik = integer(0),
  coef = c(0.8, 0),
  tAhead = 1,
  rho = 0,
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

  A numeric matrix or vector representing the values of the time series.

- series_esd:

  A numeric matrix or vector representing the error standard deviations
  of the time series.

- series_names:

  An optional character vector representing the name of the series.

- fitted_values:

  A numeric vector containing the fitted values from the model.

- loglik:

  A numeric value representing the log-likelihood of the model.

- kalmanlik:

  A numeric value representing the Kalman likelihood of the model.

- coef:

  A numeric vector of length 2 containing the autoregressive (phiR) and
  cross-correlation (phiI) coefficients of the model. Each value must be
  in \[-1, 1\] and the modulus \`sqrt(phiR^2 + phiI^2)\` must be less
  than 1. Defaults to \`c(0.8, 0)\`. Unlike \[iAR\], which uses a
  scalar, both \[CiAR\] and \`BiAR\` use a length-2 vector.

- tAhead:

  A numeric value specifying the forecast horizon (default: 1).

- rho:

  A numeric vector containing the estimated coefficients of the model.

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

The \`BiAR\` class is designed to handle bivariate irregularly observed
time series data using an autoregressive approach. It extends the
\`multidata\` class to include additional properties for modeling
bivariate time series.

Key features of the \`BiAR\` class include: - Support for bivariate time
series data. - Forecasting and interpolation functionalities for
irregular time points. - Assumptions of zero-mean and standardized
processes, configurable by the user. - Estimation of model parameters
and likelihoods, including Kalman likelihood.

## Validation Rules

\- \`@times\` must be a numeric vector without dimensions and strictly
increasing. - \`@series\` must be a numeric matrix with two columns
(bivariate) or be empty. - The number of rows in \`@series\` must match
the length of \`@times\`. - \`@series_esd\`, if provided, must be a
numeric matrix. Its dimensions must match those of \`@series\`, or it
must have one row and the same number of columns. - If \`@series_esd\`
contains NA values, they must correspond positionally to NA values in
\`@series\`. - \`@series_names\`, if provided, must be a character
vector with length equal to the number of columns in \`@series\`, and
all names must be unique. - \`@coef\` must be a numeric vector of length
2, with each element strictly between -1 and 1. - The modulus
\`sqrt(phiR^2 + phiI^2)\` must be strictly less than 1 (stationarity
condition). - \`@tAhead\` must be a strictly positive numeric scalar.

## References

Elorrieta F, Eyheramendy S, Palma W, Ojeda C (2021). “A novel bivariate
autoregressive model for predicting and forecasting irregularly observed
time series.” *Monthly Notices of the Royal Astronomical Society*,
**505**(1), 1105-1116. ISSN 0035-8711.
[doi:10.1093/mnras/stab1216](https://doi.org/10.1093/mnras/stab1216) .
https://academic.oup.com/mnras/article-pdf/505/1/1105/38391762/stab1216.pdf.

## Examples

``` r
times=gentime(n=200, distribution = "expmixture",
             lambda1 = 130, lambda2 = 6.5, p1 = 0.15, p2 = 0.85)@times
my_BiAR <- BiAR(times = times,coef = c(0.9, 0.3), rho = 0.9)

# Access properties
my_BiAR@coef
#> [1] 0.9 0.3
```
