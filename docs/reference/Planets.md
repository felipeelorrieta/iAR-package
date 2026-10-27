# Transit of an extrasolar planet

Time series corresponding to the residuals of the parametric model
fitted by Jordan et al (2013) for a transit of an extrasolar planet.

## Usage

``` r
Planets
```

## Format

A data frame with 91 observations on the following 2 variables:

- t:

  Time from mid-transit (hours).

- r:

  Residuals of the parametric model fitted by Jordan et al (2013).

## References

Jordán A, Espinoza N, Rabus M, Eyheramendy S, Sing D~K, Désert J, Bakos
G~Á, Fortney J~J, López-Morales M, Maxted P~F~L, Triaud A~H~M~J,
Szentgyorgyi A (2013). “A Ground-based Optical Transmission Spectrum of
WASP-6b.” *The Astrophysical Journal*, **778**, 184.
[doi:10.1088/0004-637X/778/2/184](https://doi.org/10.1088/0004-637X/778/2/184)
. 1310.6048.

## Examples

``` r
data(Planets)
#plot(Planets[,1],Planets[,2],xlab='Time from mid-transit (hours)',ylab='Noise',pch=20)
```
