# Polynomial Utilities

This section describes the core structures and functions used to
represent and manipulate polynomial covers.

## PolynomialCover type

The central object in the package is `PolynomialCover`, which represents
a polynomial map

```math
P : ℂ → ℂ
```
given by its complex coefficients.

```@docs
PolynomialCover
coefficients
degree
```
## Evaluation

Polynomial covers can be evaluated at complex points.

```@docs
value
```
Example:
```julia
P = PolynomialCover([1, 0, -1])   # z² - 1
value(P,2)
```






