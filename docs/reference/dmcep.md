# Double Mode Cepheid.

Time series of a double mode cepheid variable star obtained from OGLE.

## Usage

``` r
dmcep
```

## Format

A data frame with 191 observations on the following 3 variables:

- t:

  heliocentric Julian Day

- m:

  magnitude

- merr:

  measurement error of the magnitude (in mag).

## Details

The dominant frequency computed by GLS for this light curve is
0.7410152. The second frequency computed by GLS for this light curve is
0.5433353. OGLE-ID:175210

## Examples

``` r
data(dmcep)
f1=0.7410152
o1=phase(data=dmcep,f1=f1,twop=TRUE)
plot(o1@times_phased,o1@series_phased,pch=20)

#fit=harmonicfit(dmcep,f1)
#f2=0.5433353
#foldlc(cbind(dmcep$t,fit$res,dmcep$merr),f2)
```
