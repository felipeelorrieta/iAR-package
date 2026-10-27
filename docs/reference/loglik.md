# Maximum Likelihood Estimation for iAR Models

Maximum Likelihood Estimation for irregular autoregressive (iAR) models,
supporting different distribution families: normal (\`iAR\`), t
(\`iAR-T\`), and gamma (\`iAR-Gamma\`).

## Usage

``` r
loglik(x, ...)
```

## Arguments

- x:

  An object of class `iAR`, containing the model specification,
  parameters, and the time series to be evaluated:

  - `series`: The observed time series.

  - `times`: A numeric vector specifying the time points of the series.

  - `series_esd`: (Optional) A standardized version of the series.

  - `zero_mean`: Logical. Indicates whether the series should be
    mean-centered.

  - `standardized`: Logical. Indicates whether the series is
    standardized.

  - `hessian`: Logical. If `TRUE`, the function computes the Hessian
    matrix for parameter estimation.

  - `init`: Initial values for parameter estimation. For
    `family = "norm"`, a single numeric value (default 0.1). For
    `family = "t"`, a numeric vector of length two. For
    `family = "gamma"`, a numeric vector of length three.

  - `family`: The distribution family of the iAR model (one of "norm",
    "t", or "gamma").

  - `df`: Degrees of freedom for the t-distribution (only for
    `family = "t"`).

  - `sigma`: The scale parameter for the t-distribution (only for
    `family = "t"`).

  - `mean`: The mean parameter for the gamma distribution (only for
    `family = "gamma"`).

  - `variance`: The variance parameter for the gamma distribution (only
    for `family = "gamma"`).

- ...:

  Additional arguments (unused).

## Value

An updated `iAR` object with the following additional attributes:

- `coef`: Estimated model coefficients.

- `loglik`: Log-likelihood value of the model.

- `summary`: A summary table containing parameter estimates, standard
  errors, and p-values.

- `sigma`: For t and gamma families, the estimated scale parameter.

- `mean`: For the gamma family, the estimated mean parameter.

- `variance`: For the gamma family, the estimated variance parameter.

## Details

This method estimates the parameters of an iAR model using the Maximum
Likelihood Estimation (MLE) approach. Depending on the chosen
distribution family, the corresponding likelihood function is maximized:

- "norm" maximizes the likelihood for a normally-distributed series.

- "t" maximizes the likelihood for a t-distributed series.

- "gamma" maximizes the likelihood for a gamma-distributed series.

The function updates the `iAR` object with the estimated parameters, the
log-likelihood value, and a summary table that includes standard errors
and p-values.

## References

Eyheramendy S, Elorrieta F, Palma W (2018). “An irregular discrete time
series model to identify residuals with autocorrelation in astronomical
light curves.” *Monthly Notices of the Royal Astronomical Society*,
**481**(4), 4311-4322. ISSN 0035-8711.
[doi:10.1093/mnras/sty2487](https://doi.org/10.1093/mnras/sty2487) .
https://academic.oup.com/mnras/article-pdf/481/4/4311/25906473/sty2487.pdf.

## Examples

``` r
# Example: Estimating parameters for a normal iAR model
library(iAR)
times <- 1:100
model <- iAR(family = "norm", times = times, coef = 0.9, hessian = TRUE)
model <- sim(model)  # Simulate the series
model <- loglik(model)  # Estimate parameters using MLE
#> Warning: Optimization may not have converged (code = 52).
print(model@coef)  # Access the estimated coefficients
#> [1] 0.8957576
print(model@loglik)  # Access the computed log-likelihood
#> [1] 60.81373
```
