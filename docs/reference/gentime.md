# Generating Irregularly spaced times

A method for generating time points based on a statistical distribution.
The results are stored in the \`times\` slot of the \`utilities\`
object.

## Usage

``` r
gentime(x = NULL, ...)
```

## Arguments

- x:

  An object of class \`utilities\`, or \`NULL\` (the default). When
  \`x\` is \`NULL\` or omitted, the \`utilities\` object is created
  internally, so \`gentime(n = 200)\` is equivalent to the three-step
  pattern \`o \<- utilities(); o \<- gentime(o, n = 200); o\`.

- ...:

  Additional arguments for generating time points:

  n

  :   A positive integer. Length of observation times.

  distribution

  :   A character string specifying the distribution of the observation
      times. Default is \`"expmixture"\`. Available options are: -
      \`"expmixture"\`: A mixture of two exponential distributions. -
      \`"uniform"\`: A uniform distribution. - \`"exponential"\`: A
      single exponential distribution. - \`"gamma"\`: A gamma
      distribution.

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

An updated \`utilities\` object with the generated observation times
stored in the \`times\` slot.

## References

Eyheramendy S, Elorrieta F, Palma W (2018). “An irregular discrete time
series model to identify residuals with autocorrelation in astronomical
light curves.” *Monthly Notices of the Royal Astronomical Society*,
**481**(4), 4311-4322. ISSN 0035-8711.
[doi:10.1093/mnras/sty2487](https://doi.org/10.1093/mnras/sty2487) .
https://academic.oup.com/mnras/article-pdf/481/4/4311/25906473/sty2487.pdf.

## Examples

``` r
set.seed(12917)
st=gentime(n=200, distribution = "expmixture", lambda1 = 130,
lambda2 = 6.5, p1 = 0.15, p2 = 0.85)@times
mean(diff(st))
#> [1] 28.4244

st=gentime(n=200, distribution = "expmixture", lambda1 = 15,
lambda2 = 2.5, p1 = 0.15, p2 = 0.85)@times
mean(diff(st))
#> [1] 4.592272
```
