# Normalize an object to the fitted-margin protocol

`as_margin()` validates fitted margins and adapts the legacy
`list(distr = ...)` representation through
[`stats_margin()`](https://vinecopulib.github.io/rvinecopulib/reference/stats_margin.md).

## Usage

``` r
as_margin(margin)
```

## Arguments

- margin:

  a fitted margin or legacy stats distribution list.

## Value

A validated fitted margin.
