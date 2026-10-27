# ZTF g-band Cataclysmic Variable/Nova

Time series of a cataclysmic variable/nova object observed in the g-band
of the ZTF survey and processed by the ALeRCE broker.ZTF Object code:
ZTF18aayzpbr

## Usage

``` r
cvnovag
```

## Format

A data frame with 67 observations on the following 3 variables:

- t:

  heliocentric Julian Day - 2400000

- m:

  magnitude

- merr:

  measurement error standard deviations.

## References

Förster F, Cabrera-Vives G, Castillo-Navarrete E, Estévez PA,
Sánchez-Sáez P, Arredondo J, Bauer FE, Carrasco-Davis R, Catelan M,
Elorrieta F, Eyheramendy S, Huijse P, Pignata G, Reyes E, Reyes I,
Rodríguez-Mancini D, Ruz-Mieres D, Valenzuela C, Álvarez-Maldonado I,
Astorga N, Borissova J, Clocchiatti A, Cicco DD, Donoso-Oliva C,
Hernández-García L, Graham MJ, Jordán A, Kurtev R, Mahabal A, Maureira
JC, Muñoz-Arancibia A, Molina-Ferreiro R, Moya A, Palma W,
Pérez-Carrasco M, Protopapas P, Romero M, Sabatini-Gacitua L, Sánchez A,
Martín JS, Sepúlveda-Cobo C, Vera E, Vergara JR (2021). “The Automatic
Learning for the Rapid Classification of Events (ALeRCE) Alert Broker.”
*The Astronomical Journal*, **161**(5), 242.
[doi:10.3847/1538-3881/abe9bc](https://doi.org/10.3847/1538-3881/abe9bc)
.

## Examples

``` r
data(cvnovag)
plot(cvnovag$t,cvnovag$m,type="l",ylab="",xlab="",col="green")
```
