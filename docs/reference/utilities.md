# Utilities Class

The \`utilities\` class is an S7 class designed to group utility
functions that are not directly tied to other specific objects in the
package. These include the functions \`gentime\`, \`paringits\`,
\`harmonicfit\`, and \`foldlc\`.

## Usage

``` r
utilities(
  times = integer(0),
  series = integer(0),
  series_esd = integer(0),
  times_phased = integer(0),
  series_phased = integer(0),
  series_esd_phased = integer(0),
  fitted_values = integer(0),
  residuals = integer(0),
  coef = integer(0),
  summary = list(),
  paired = integer(0)
)
```

## Arguments

- times:

  A numeric vector storing the time points.

- series:

  A numeric vector representing the values of the time series.

- series_esd:

  A numeric vector representing the error standard deviations of the
  time series.

- times_phased:

  A numeric vector of phased times.

- series_phased:

  A numeric vector containing the time series values ordered by phase.

- series_esd_phased:

  A numeric vector containing the error standard deviations of the time
  series ordered by phase.

- fitted_values:

  A numeric vector containing the fitted values from the harmonic model.

- residuals:

  A numeric vector containing the residuals from the harmonic model.

- coef:

  A numeric vector containing the estimated coefficients of the harmonic
  model.

- summary:

  A summary object containing the harmonic model information.

- paired:

  Data Frame with the paired datasets.

## Description

This class acts as a container for standalone methods that perform
independent operations within the package. By grouping them under a
single class, the package achieves better modularity and organization,
facilitating maintenance and extensibility.

## Available Methods

\- \`gentime\`: Generates time points based on a specified statistical
distribution. - \`paringits\`: A method for pairing irregular time
series. - \`harmonicfit\`: A method for fitting harmonic models to
data. - \`foldlc\`: A method for folding light curves in time series
analysis.

## Examples

``` r
# Create a utilities object
```
