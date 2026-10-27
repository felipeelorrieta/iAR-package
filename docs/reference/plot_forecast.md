# Plot Forecast Method for unidata Objects

This method visualizes the \`unidata\` and \`multidata\` objects by
plotting both the original time series and the forecasted values. The
\`unidata\` (\`multidata\`) object is converted into a \`zoo\` object
for plotting, and the forecasted values are added as points with a
different color.

## Usage

``` r
plot_forecast(x,...)
```

## Arguments

- x:

  An object of class \`unidata\` or \`multidata\`. This object must
  contain the time series in the \`series\` slot, the associated time
  points in the \`times\` slot, the forecasted values in the
  \`forecast\` slot, and the forecast horizon in the \`tAhead\` slot.

- ...:

  Additional arguments passed to the \`plot.zoo\` function for
  customizing the plot (e.g., titles, colors, etc.).
