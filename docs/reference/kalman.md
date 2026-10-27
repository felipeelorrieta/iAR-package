# Maximum Likelihood Estimation of Parameters for iAR, CiAR, and BiAR Models using the Kalman Filter

Performs Maximum Likelihood Estimation (MLE) of the parameters of the
iAR, CiAR, and BiAR models by maximizing the likelihood function using
the Kalman filter. This method applies the Kalman filter to compute the
likelihood and estimate the model parameters that maximize the
likelihood for each model type.

## Usage

``` r
kalman(x, ...)
```

## Arguments

- x:

  An object of class `iAR`, `CiAR`, or `BiAR` containing the model
  specification and parameters:

  - For `iAR` models, the object should contain:

    - `series`: A numeric vector of the time series data.

    - `times`: A numeric vector specifying the time points of the
      series.

    - `series_esd`: A numeric vector specifying the error structure.

    - `zero_mean`: A logical value indicating if the series should be
      zero-centered.

    - `standardized`: A logical value indicating if the series should be
      standardized.

    - `init`: A numeric value specifying the initial value for parameter
      estimation. Default is 0.1.

  - For `CiAR` models, the object should contain:

    - `series`: A numeric vector of the time series data.

    - `times`: A numeric vector specifying the time points of the
      series.

    - `series_esd`: A numeric vector specifying the error structure.

    - `zero_mean`: A logical value indicating if the series should be
      zero-centered.

    - `standardized`: A logical value indicating if the series should be
      standardized.

    - `c`: A numeric value specifying the scale parameter for the CiAR
      model.

    - `init`: A numeric vector of length two specifying the initial
      values for parameter estimation. When provided, these values are
      included as an additional starting point.

  - For `BiAR` models, the object should contain:

    - `series`: A numeric matrix containing two columns for bivariate
      time series data.

    - `times`: A numeric vector specifying the time points of the
      series.

    - `series_esd`: A numeric matrix specifying the error structure for
      both series.

    - `zero_mean`: A logical value indicating if the series should be
      zero-centered.

    - `init`: A numeric vector of length two specifying the initial
      values for parameter estimation. When provided, these values are
      included as an additional starting point.

- ...:

  Additional arguments (unused).

## Value

An updated object of class `iAR`, `CiAR`, or `BiAR`, where the `coef`
property is updated with the estimated model parameters (using MLE) and
the `kalmanlik` property contains the log-likelihood value of the model.

## Details

This function applies the Kalman filter to perform Maximum Likelihood
Estimation (MLE) of the parameters of the autoregressive models (iAR,
CiAR, BiAR). The Kalman filter is used to maximize the likelihood
function based on the given time series data, and the parameters that
maximize the likelihood are estimated.

\- For `iAR`, the Kalman filter is applied to estimate the model
parameters by maximizing the likelihood. - For `CiAR`, the Kalman filter
is applied to estimate the parameters of the complex autoregressive
model by maximizing the likelihood. - For `BiAR`, the Kalman filter is
applied to estimate the parameters of the bivariate autoregressive model
by maximizing the likelihood.

The method returns the updated model object, including the estimated
parameters and the log-likelihood value.

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
# Example 1: Applying Kalman filter for MLE of iAR model parameters
library(iAR)
n=100
set.seed(6714)
times=gentime(n=n)@times
model_norm <- iAR(family = "norm", times = times, coef = 0.9,hessian=TRUE)
model_norm <- sim(model_norm)
model_norm <- kalman(model_norm)
print(model_norm@coef)  # Access the estimated coefficients
#> [1] 0.9088113
print(model_norm@kalmanlik)  # Access the Kalman likelihood value
#> [1] 0.2431375

# Example 2: Applying Kalman filter for MLE of CiAR model parameters
set.seed(6714)
model_CiAR <- CiAR(times = times,coef = c(0.9, 0))
model_CiAR <- sim(model_CiAR)
y=model_CiAR@series
y1=y/sd(y)
model_CiAR@series=y1
model_CiAR@series_esd=rep(0,n)
model_CiAR <- kalman(model_CiAR)
print(model_CiAR@coef)
#> [1] 0.9199923 0.1072453
print(model_CiAR@kalmanlik)
#> [1] 0.2272269

# Example 3: Applying Kalman filter for MLE of BiAR model parameters
set.seed(6714)
model_BiAR <- BiAR(times = times,coef = c(0.9, 0.3), rho = 0.9)
model_BiAR <- sim(model_BiAR)
y=model_BiAR@series
y1=y/apply(y,2,sd)
model_BiAR@series=y1
model_BiAR@series_esd=matrix(0,n,2)
model_BiAR <- kalman(model_BiAR)
print(model_BiAR@coef)
#> [1] 0.8997529 0.3154864
print(model_BiAR@kalmanlik)
#> [1] -27.07637
```
