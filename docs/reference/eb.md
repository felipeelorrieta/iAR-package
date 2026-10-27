# Eclipsing Binaries (Beta Lyrae)

Time series of a Beta Lyrae variable star obtained from OGLE.

## Usage

``` r
eb
```

## Format

A data frame with 470 observations on the following 3 variables:

- t:

  heliocentric Julian Day

- m:

  magnitude

- merr:

  measurement error of the magnitude (in mag).

## Details

The frequency computed by GLS for this light curve is 1.510571586.
Catalogs and designations of this star:OGLE051951.22-694002.7

## Examples

``` r
data(eb)
f1=1.510571586
o1=phase(data=eb,f1=f1,twop=TRUE)
plot(o1@times_phased,o1@series_phased,pch=20)
```
