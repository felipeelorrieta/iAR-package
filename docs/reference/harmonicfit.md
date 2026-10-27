# Harmonic Fit to Time Series

This function fit an k-harmonic function to time series data.

## Usage

``` r
harmonicfit(x, ...)
```

## Arguments

- x:

  An object of class \`utilities\`, or \`NULL\` (the default). When
  \`x\` is \`NULL\` or omitted, the \`utilities\` object is created
  internally, so \`harmonicfit(data = clcep, f1 = f1)\` is the
  recommended usage.

- ...:

  Additional arguments for pairing time series:

  data

  :   A data frame with three columns corresponding to the time, values,
      and standard errors of the irregularly observed time series.

  f1

  :   frequency (1 / period) of the time series.

  nham

  :   Number of harmonic components in the model.

  weighted

  :   logical; if true, a weighted least squares (WLS) estimation is
      performed using weights based on the standard deviations of the
      errors. Default is 'FALSE'.

  remove_trend

  :   logical; if true, the linear trend of time series will be removed
      before the the harmonic model is fitted.

## Value

An object of class \`utilities\` with the slots:

- fitted_values:

  Fitted values from the harmonic model.

- residuals:

  Residuals from the harmonic model.

- coef:

  Estimated coefficients of the harmonic model.

- summary:

  A summary object containing detailed model information.

## Details

The function fits a harmonic regression model with 'nham' components to
the input time series, optionally removing a linear trend and allowing
for weighted estimation when standard errors are available.

## References

Elorrieta F, Eyheramendy S, Palma W, Ojeda C (2021). “A novel bivariate
autoregressive model for predicting and forecasting irregularly observed
time series.” *Monthly Notices of the Royal Astronomical Society*,
**505**(1), 1105-1116. ISSN 0035-8711.
[doi:10.1093/mnras/stab1216](https://doi.org/10.1093/mnras/stab1216) .
https://academic.oup.com/mnras/article-pdf/505/1/1105/38391762/stab1216.pdf.

## Examples

``` r
data(clcep)
f1=0.060033386
o1=harmonicfit(data=clcep,f1=f1)
head(o1@fitted_values)
#>           1           2           3           4           5           6 
#> -0.01307640  0.16116566 -0.08189031 -0.08047315 -0.07301615  0.26554040 
head(o1@residuals)
#>            1            2            3            4            5            6 
#> -0.020520255 -0.022710532 -0.016366868  0.028516019  0.003159220 -0.002247028 
```
