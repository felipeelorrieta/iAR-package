# Summary of iAR package models

This method provides a summary for objects created by the iAR package,
including \`iAR\`, \`CiAR\`, \`BiAR\`, and \`utilities\` objects. For
fitted model objects, the summary reports model-specific information
such as the model family, coefficients, standard errors, t-values, and
p-values. The output is formatted to display the relevant statistics
based on the model class.

## Usage

``` r
summary(object, ...)
```

## Arguments

- object:

  An object of class \`iAR\`, \`CiAR\`, \`BiAR\`, or \`utilities\`. For
  fitted model objects, this contains model-specific parameters such as
  coefficients (\`coef\`), standard errors (\`stderr\`), t-values
  (\`tvalue\`), p-values (\`pvalue\`), and other class-specific
  components.

- ...:

  Additional arguments (unused).
