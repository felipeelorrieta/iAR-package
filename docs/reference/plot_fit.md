# Plot Fit Method for unidata Objects

This method visualizes the \`unidata\` and \`multidata\` objects by
plotting both the original time series and the fitted values. The
\`unidata\` (\`multidata\`) object is converted into a \`zoo\` object
for plotting, and the fitted values are added as a secondary line with a
different style and color.

## Usage

``` r
plot_fit(x,...)
```

## Arguments

- x:

  An object of class \`unidata\` or \`multidata\`. This object must
  contain the time series in the \`series\` slot, the associated time
  points in the \`times\` slot, and the fitted values in the
  \`fitted_values\` slot.

- ...:

  Additional arguments passed to the \`plot.zoo\` function for
  customizing the plot (e.g., titles, colors, etc.).
