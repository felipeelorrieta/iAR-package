# Simulate Time Series for Multiple iAR Models

Simulates a time series for irregular autoregressive (iAR) models,
including: 1. Normal iAR model (\`iAR\`) 2. T-distribution iAR model
(\`iAR-T\`) 3. Gamma-distribution iAR model (\`iAR-Gamma\`)

## Usage

``` r
sim(x, ...)
```

## Arguments

- x:

  An object of class `iAR`, `CiAR`, or `BiAR`, containing the model
  specification and parameters:

  - For `iAR` (irregular AR models), the model family could be "norm",
    "t", or "gamma", where:

    - `family`: The distribution family of the iAR model (one of "norm",
      "t", or "gamma").

    - `coef`: The coefficient(s) of the iAR model.

    - `times`: A numeric vector specifying the time points of the
      series.

    - `df`: Degrees of freedom for the t-distribution (only for
      `family = "t"`).

    - `sigma`: The scale parameter for the t-distribution (only for
      `family = "t"`).

    - `mean`: The mean parameter for the gamma distribution (only for
      `family = "gamma"`).

    - `variance`: The variance parameter for the gamma distribution
      (only for `family = "gamma"`).

  - For `CiAR` (complex irregular autoregressive models):

    - `coef`: The real and imaginary parts of the CiAR model's
      coefficients.

    - `times`: A numeric vector specifying the time points of the
      series.

    - `rho`: The correlation parameter for the CiAR model.

    - `c`: The scale parameter for the CiAR model.

  - For `BiAR` (BiAR models):

    - `coef`: The coefficients of the BiAR model (real and imaginary).

    - `times`: A numeric vector specifying the time points of the
      series.

    - `rho`: The correlation parameter for the BiAR model.

    - `series_esd`: The series for the error structure (optional, used
      internally).

- ...:

  Additional arguments for generating irregular times. This is used only
  if no time points have been provided in the `times` argument:

  n

  :   A positive integer. Length of observation times. Default is 100.
      This is used only if no time points are provided in `times`.

  distribution

  :   A character string specifying the distribution of the observation
      times. Default is \`"expmixture"\`. Available options are:

      expmixture

      :   A mixture of two exponential distributions.

      uniform

      :   A uniform distribution.

      exponential

      :   A single exponential distribution.

      gamma

      :   A gamma distribution.

  lambda1

  :   Mean (1/rate) of the exponential distribution or the first
      exponential distribution in a mixture of exponential
      distributions. Default is \`130\`.

  lambda2

  :   Mean (1/rate) of the second exponential distribution in a mixture
      of exponential distributions. Default is \`6.5\`.

  p1

  :   Weight of the first exponential distribution in a mixture of
      exponential distributions. Default is \`0.15\`.

  p2

  :   Weight of the second exponential distribution in a mixture of
      exponential distributions. Default is \`0.85\`.

  a

  :   Shape parameter of a gamma distribution or lower limit of the
      uniform distribution. Default is \`0\`.

  b

  :   Scale parameter of a gamma distribution or upper limit of the
      uniform distribution. Default is \`1\`.

## Value

An updated object of class `iAR`, `CiAR`, or `BiAR`, where the `series`
property contains the simulated time series.

## Details

This function simulates time series based on the specified model and its
parameters. Depending on the class of the input object:

- For `iAR` models, it supports three distribution families:

- "norm" for normal distribution.

- "t" for t-distribution.

- "gamma" for gamma distribution.

- For `CiAR` models, it uses complex autoregressive processes to
  generate the time series.

- For `BiAR` models, it simulates a BiAR process using specified
  coefficients and correlation.

The coefficients and any family-specific parameters must be set before
calling this function.

## References

Eyheramendy S, Elorrieta F, Palma W (2018). “An irregular discrete time
series model to identify residuals with autocorrelation in astronomical
light curves.” *Monthly Notices of the Royal Astronomical Society*,
**481**(4), 4311-4322. ISSN 0035-8711.
[doi:10.1093/mnras/sty2487](https://doi.org/10.1093/mnras/sty2487) .
https://academic.oup.com/mnras/article-pdf/481/4/4311/25906473/sty2487.pdf.
,Elorrieta, F, Eyheramendy, S, Palma, W (2019). “Discrete-time
autoregressive model for unequally spaced time-series observations.”
*A&A*, **627**, A120.
[doi:10.1051/0004-6361/201935560](https://doi.org/10.1051/0004-6361/201935560)
. ,Elorrieta F, Eyheramendy S, Palma W, Ojeda C (2021). “A novel
bivariate autoregressive model for predicting and forecasting
irregularly observed time series.” *Monthly Notices of the Royal
Astronomical Society*, **505**(1), 1105-1116. ISSN 0035-8711.
[doi:10.1093/mnras/stab1216](https://doi.org/10.1093/mnras/stab1216) .
https://academic.oup.com/mnras/article-pdf/505/1/1105/38391762/stab1216.pdf.

## Examples

``` r
# Example 1: Simulating a normal iAR model
library(iAR)
n=100
set.seed(6714)
times=gentime(n=n)@times
model_norm <- iAR(family = "norm", times = times, coef = 0.9,hessian=TRUE)
model_norm <- sim(model_norm)
plot(model_norm, type = "l", main = "Simulated iAR-Norm Series")


# Example 2: Simulating a CiAR model
set.seed(6714)
model_CiAR <- CiAR(times = times,coef = c(0.9, 0))
model_CiAR <- sim(model_CiAR)
plot(model_CiAR , type = "l", main = "Simulated CiAR Series")


# Example 3: Simulating a BiAR model
set.seed(6714)
model_BiAR <- BiAR(times = times,coef = c(0.9, 0.3), rho = 0.9)
model_BiAR <- sim(model_BiAR)
plot(times, model_BiAR@series[,1], type = "l", main = "Simulated BiAR Series")

```
