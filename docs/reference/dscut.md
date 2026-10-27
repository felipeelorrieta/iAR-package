# Delta Scuti

Time series of a Delta Scuti variable star obtained from HIPPARCOS.

## Usage

``` r
dscut
```

## Format

A data frame with 116 observations on the following 3 variables:

- t:

  heliocentric Julian Day

- m:

  magnitude

- merr:

  measurement error of the magnitude (in mag).

## Details

The frequency computed by GLS for this light curve is 14.88558646.
Catalogs and designations of this star: HD 1989: HD 199757 TYCHO-2 2000:
TYC 7973-401-1 USNO-A2.0: USNO-A2 0450-39390397 HIP: HIP 103684

## Examples

``` r
data(dscut)
f1=14.88558646
o1=phase(data=dscut,f1=f1,twop=TRUE)
plot(o1@times_phased,o1@series_phased,pch=20)
```
