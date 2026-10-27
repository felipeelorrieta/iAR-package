# Active Galactic Nuclei

Time series of the AGN MCG-6-30-15 measured in the K-band between 2006
August and 2011 July with the ANDICAM camera mounted on the 1.3 m
telescope at Cerro Tololo Inter-American Observatory (CTIO)

## Usage

``` r
agn
```

## Format

A data frame with 237 observations on the following 3 variables:

- t:

  heliocentric Julian Day - 2450000

- m:

  Flux \$(10^(-15) ergs/s/cm^2 /A)\$

- merr:

  measurement error standard deviations.

## References

Lira P, Arévalo P, Uttley P, McHardy IMM, Videla L (2015). “Long-term
monitoring of the archetype Seyfert galaxy MCG-6-30-15: X-ray, optical
and near-IR variability of the corona, disc and torus.” *Monthly Notices
of the Royal Astronomical Society*, **454**(1), 368-379. ISSN 0035-8711.
[doi:10.1093/mnras/stv1945](https://doi.org/10.1093/mnras/stv1945) .

## Examples

``` r
data(agn)
plot(agn$t,agn$m,type="l",ylab="",xlab="")
```
