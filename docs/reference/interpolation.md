# Interpolation for iAR, CiAR, and BiAR Classes

This method performs imputation of missing values in a time series using
an autoregressive model. The imputation is done iteratively for each
missing value, utilizing available data and model coefficients.
Depending on the model family, the interpolation is performed
differently: - For `norm`: A standard autoregressive model for normally
distributed data. - For `t`: A model for time series with t-distributed
errors. - For `gamma`: A model for time series with gamma-distributed
errors. - For `CiAR`: A complex irregular autoregressive model. - For
`BiAR`: A bivariate autoregressive model.

## Usage

``` r
interpolation(x, ...)
```

## Arguments

- x:

  An object of class `iAR`, `CiAR`, or `BiAR`, containing the model
  specification and parameters:

  - For `iAR`:

    - `family`: The distribution family of the iAR model (one of "norm",
      "t", or "gamma").

    - `series`: A numeric vector representing the time series to
      interpolate.

    - `coef`: The coefficients of the iAR model.

    - `zero_mean`: Logical, whether the model assumes a zero-mean
      series.

    - `standardized`: Logical, whether the model uses standardized data
      (only for "norm" family).

    - `mean`: The mean parameter (only for "gamma" family).

  - For `CiAR`:

    - `coef`: The coefficients of the CiAR model.

    - `series_esd`: The series of error standard deviations for the CiAR
      model.

    - `zero_mean`: Logical, whether the model assumes a zero-mean
      series.

    - `standardized`: Logical, whether the model uses standardized data.

    - `seed`: Optional seed for random number generation.

  - For `BiAR`:

    - `coef`: The coefficients of the BiAR model.

    - `series_esd`: The series of error standard deviations for the BiAR
      model.

    - `zero_mean`: Logical, whether the model assumes a zero-mean
      series.

    - `yini1`: Initial value for the first time series (for BiAR
      models).

    - `yini2`: Initial value for the second time series (for BiAR
      models).

    - `seed`: Optional seed for random number generation.

- ...:

  Additional arguments (unused).

## Value

An object of the same class as `x` with interpolated time series.

## Details

Performs interpolation on time series with missing values. This method
is implemented for: 1. Irregular Autoregressive models (iAR) 2. Complex
Irregular Autoregressive models (CiAR) 3. Bivariate Autoregressive
models (BiAR)

The method handles missing values (NA) in the time series by imputing
them iteratively. For each missing value, the available data is used to
fit the autoregressive model, and the missing value is imputed based on
the model's output. For the `CiAR` and `BiAR` models, the error standard
deviations and initial values are also considered during imputation.

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

## See also

[`forecast`](https://felipeelorrieta.github.io/iAR-package/reference/forecast.md)
for forecasting methods for these models.

## Examples

``` r
# Interpolation for iAR model
library(iAR)
n=100
set.seed(6714)
times=gentime(n=n)@times
model_norm <- iAR(family = "norm", times = times, coef = 0.9)  
model_norm <- sim(model_norm)
y=model_norm@series
y1=y/sd(y)
model_norm@series=y1
model_norm@series_esd=rep(0,n)
model_norm <- kalman(model_norm) 
print(model_norm@coef)
#> [1] 0.8850837
napos=10
model_norm@series[napos]=NA
model_norm <- interpolation(model_norm)
interpolation=na.omit(model_norm@interpolated_values)
mse=as.numeric(y1[napos]-interpolation)^2
print(mse)
#> [1] 2.019834
plot(times,y,type='l',xlim=c(times[napos-5],times[napos+5]))
points(times,y,pch=20)
points(times[napos],interpolation*sd(y),col="red",pch=20)


# Interpolation for CiAR model
model_CiAR <- CiAR(times = times,coef = c(0.9, 0))
model_CiAR <- sim(model_CiAR)
y=model_CiAR@series
y1=y/sd(y)
model_CiAR@series=y1
model_CiAR@series_esd=rep(0,n)
model_CiAR <- kalman(model_CiAR)
print(model_CiAR@coef)
#> [1]  9.182499e-01 -6.952136e-10
napos=10
model_CiAR@series[napos]=NA
model_CiAR <- interpolation(model_CiAR)
interpolation=na.omit(model_CiAR@interpolated_values)
mse=as.numeric(y1[napos]-interpolation)^2
print(mse)
#> [1] 2.333803
plot(times,y,type='l',xlim=c(times[napos-5],times[napos+5]))
points(times,y,pch=20)
points(times[napos],interpolation*sd(y),col="red",pch=20)

# Interpolation for BiAR model
model_BiAR <- BiAR(times = times,coef = c(0.9, 0.3), rho = 0.9)
model_BiAR <- sim(model_BiAR)
y=model_BiAR@series
y1=y/apply(y,2,sd)
model_BiAR@series=y1
model_BiAR@series_esd=matrix(0,n,2)
model_BiAR <- kalman(model_BiAR)
print(model_BiAR@coef) 
#> [1] 0.9031073 0.3014462
napos=10
model_BiAR@series[napos,1]=NA
model_BiAR@series_esd[napos,1]=NA
model_BiAR <- interpolation(model_BiAR)
interpolation=na.omit(model_BiAR@interpolated_values[,1])
mse=as.numeric(y1[napos,1]-interpolation)^2
print(mse)
#> [1] 0.3141242
par(mfrow=c(2,1))
plot(times,y[,1],type='l',xlim=c(times[napos-5],times[napos+5]))
points(times,y[,1],pch=20)
points(times[napos],interpolation*apply(y,1,sd)[1],col="red",pch=20)
plot(times,y[,2],type='l',xlim=c(times[napos-5],times[napos+5]))
points(times,y[,2],pch=20)
```
