# Classical Cepheid

Time series of a classical cepheid variable star obtained from
HIPPARCOS.

## Usage

``` r
clcep
```

## Format

A data frame with 109 observations on the following 3 variables:

- t:

  heliocentric Julian Day

- m:

  magnitude

- merr:

  measurement error of the magnitude (in mag).

## Details

The frequency computed by GLS for this light curve is 0.060033386.
Catalogs and designations of this star: HD 1989: HD 305996 TYCHO-2
2000:TYC 8958-2333-1 USNO-A2.0:USNO-A2 0225-10347916 HIP: HIP-54101

## Examples

``` r
data(clcep)
f1=0.060033386
o1=phase(data=clcep,f1=f1,twop=TRUE)
plot(o1@times_phased,o1@series_phased,pch=20)
```
