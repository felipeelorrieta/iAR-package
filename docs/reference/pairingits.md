# Pairing two irregularly observed time series

This method pairs the observational times of two irregularly observed
time series.

## Usage

``` r
pairingits(x, ...)
```

## Arguments

- x:

  An object of class \`utilities\`, or \`NULL\` (the default). When
  \`x\` is \`NULL\` or omitted, the \`utilities\` object is created
  internally, so \`pairingits(data1 = ts1, data2 = ts2)\` is the
  recommended usage.

- ...:

  Additional arguments for pairing time series:

  data1

  :   A data frame with three columns corresponding to the first
      irregularly observed time series.

  data2

  :   A data frame with three columns corresponding to the second
      irregularly observed time series.

  tol

  :   A numeric value indicating the tolerance parameter.

## Value

An object of class \`utilities\` with the slots:

- series:

  A matrix containing the paired time series, where unmatched
  measurements are filled with \`NA\`.

- series_esd:

  A matrix containing the paired error standard deviations of the time
  series, where unmatched measurements are filled with \`NA\`.

- times:

  A numeric vector with the paired observational times.

## Details

The method checks the observational times in both input time series and
pairs the measurements if they fall within the specified tolerance
(\`tol\`). If a measurement in one series cannot be paired, it is filled
with \`NA\` values for the corresponding columns of the other series.

## References

Elorrieta F, Eyheramendy S, Palma W, Ojeda C (2021). “A novel bivariate
autoregressive model for predicting and forecasting irregularly observed
time series.” *Monthly Notices of the Royal Astronomical Society*,
**505**(1), 1105-1116. ISSN 0035-8711.
[doi:10.1093/mnras/stab1216](https://doi.org/10.1093/mnras/stab1216) .
https://academic.oup.com/mnras/article-pdf/505/1/1105/38391762/stab1216.pdf.

## Examples

``` r
data(cvnovag)
data(cvnovar)
o1=pairingits(data1=cvnovag, data2=cvnovar, tol=0.1)
pargr1=na.omit(o1@paired)
st=apply(pargr1[,c(1,4)],1,mean)
model_BiAR <- BiAR(times = st,series=pargr1[,c(2,5)],series_esd=pargr1[,c(3,6)])
model_BiAR <- kalman(model_BiAR)
```
