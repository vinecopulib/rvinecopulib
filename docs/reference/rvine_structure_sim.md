# Simulate R-vine structures

Simulates from a uniform distribution over all R-vine structures on d
variables. `rvine_structure_sim()` returns an `rvine_structure()`
object, `rvine_matrix_sim()` an `rvine_matrix()`.

## Usage

``` r
rvine_structure_sim(d, natural_order = FALSE)

rvine_matrix_sim(d, natural_order = FALSE)
```

## Arguments

  - d:
    
    the number of variables

  - natural\_order:
    
    boolean; whether the structures should be in natural order
    (counter-diagonal is `1:d`).

## See also

`rvine_structure()`, `rvine_matrix()`, `plot.rvine_structure()`,
`plot.rvine_matrix()`

## Examples

``` r
rvine_structure_sim(10)
#> 10-dimensional R-vine structure ('rvine_structure')
#>  9  9  8 10  9  3  9  3  3  3
#> 10  3  9  9 10  4  3  9  9   
#>  8 10 10  3  3  9  4  4      
#>  5  4  3  4  4 10 10         
#>  3  2  4  2  2  2            
#>  4  1  2  8  8               
#>  2  8  1  1                  
#>  1  5  5                     
#>  7  7                        
#>  6                           

rvine_structure_sim(10, natural_order = TRUE)  # counter-diagonal is 1:d
#> 10-dimensional R-vine structure ('rvine_structure')
#>  7  7 10 10 10  8 10 10 10 10
#> 10 10  8  7  7 10  9  9  9   
#>  9  9  6  9  9  9  8  8      
#>  5  8  9  8  8  7  7         
#>  8  4  7  6  6  6            
#>  6  6  5  5  5               
#>  4  5  4  4                  
#>  3  3  3                     
#>  2  2                        
#>  1                           

rvine_matrix_sim(10)
#> 10-dimensional R-vine matrix ('rvine_matrix')
#>  5  3  6  3  3  6  5  3  3  3
#>  3  6  3  6  6  3  3  5  5   
#>  8  5 10 10  5  5  6  6      
#>  6  8  5  5  8  8  8         
#> 10 10  8  8 10 10            
#>  9  7  4  4  4               
#>  4  4  7  7                  
#>  7  9  9                     
#>  1  1                        
#>  2                           
```
