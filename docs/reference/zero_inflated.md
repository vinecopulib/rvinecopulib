# Declare zero-inflated data

Marks a numeric vector as continuous with an atom at zero. The class is
preserved when stored in or subset from a data frame and is detected by
`vine()`. Alternatively, use `var_types = "zi"`.

## Usage

``` r
zero_inflated(x)
```

## Arguments

  - x:
    
    a numeric vector.

## Value

`x` with an additional `zero_inflated` class.
