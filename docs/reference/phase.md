# Computing phased time series

This function computes a time series folded on its period.

## Usage

``` r
phase(x, ...)
```

## Arguments

- x:

  An object of class \`utilities\`, or \`NULL\` (the default). When
  \`x\` is \`NULL\` or omitted, the \`utilities\` object is created
  internally, so \`phase(data = clcep, f1 = f1)\` is the recommended
  usage.

- ...:

  Additional arguments for pairing time series:

  data

  :   A data frame with three columns corresponding to the time, values,
      and standard errors of the irregularly observed time series.

  f1

  :   frequency (1 / period) of the time series.

  twop

  :   logical; if TRUE, the phased series will be duplicated over two
      cycles (0–2).

## Value

An object of class \`utilities\` with the slots:

- series_phased:

  A numeric vector containing the time series values ordered by phase.

- series_esd_phased:

  A numeric vector containing the error standard deviations of the time
  series ordered by phase.

- times_phased:

  A numeric vector of phased times (values between 0 and 1, or 0 and 2
  if 'two.cycles = TRUE').

## Details

The phase \\\phi\\ of an observation is computed as \$\$\phi = \frac{t -
t_0}{p} - \mathrm{E}(t),\$\$ where \\t_0\\ is the reference time (by
default the first observation), \\p = 1/f_1\\ is the period, and
\\\mathrm{E}(t)\\ is the integer part of \\(t - t_0)/p\\.

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
o1=phase(data=clcep,f1=f1,twop=TRUE)
plot(o1@times_phased,o1@series_phased,pch=20)
```
