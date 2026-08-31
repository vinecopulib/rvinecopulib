# Normalize an object to the fitted-margin protocol

`as_margin()` validates fitted margins and adapts the legacy `list(distr
= ...)` representation through `stats_margin()`.

## Usage

``` r
as_margin(margin)
```

## Arguments

  - margin:
    
    a fitted margin or legacy stats distribution list.

## Value

A validated fitted margin.
