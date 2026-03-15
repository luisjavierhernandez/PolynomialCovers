# PolynomialCovers.jl

`PolynomialCovers.jl` provides numerical algorithms for computing fibers of polynomial maps

```math
P(x) = y
```

using robust geometric constructions.

These methods are particularly useful when studying the global
structure of polynomial maps and their monodromy.

---

# Installation

```julia
using Pkg
Pkg.add(url="https://github.com/luisjavierhernandez/PolynomialCovers.jl")
```

## Basic Example

The following example computes the roots of a polynomial equation using
the cyclic lifting solver.

```julia

using PolynomialCovers

# Define polynomial P(x) = x^5 - x + 1
P = PolynomialCover([1, 0, 0, 0, -1, 1])

# Target value
y = 5

# Compute roots of P(x) = y
roots = cyclic_lift_solver(P, y)

println(roots)
```
## Package Structure

The package is organized into several components:
	•	Geometry
Geometric primitives in the complex plane such as disks,
vertices, and polygonal paths.
	•	Polynomial Covers
Representation of polynomial maps and their critical geometry.
	•	Lifting Algorithms
Numerical routines for lifting paths through polynomial maps.
	•	Solvers
Algorithms for computing fibers of polynomial equations.


## Documentation

```@contents
Pages = [
    "geometry.md",
    "api.md"
]
```

## Mathematical Background

The algorithms implemented in PolynomialCovers.jl are based on the
idea that the solutions of

```math
P(x) = y
```

can be computed by lifting paths in the complex plane through the
polynomial map

```math
P \colon \mathbb{C} \to \mathbb{C}.
```

In particular, the package implements numerical algorithms inspired by
monodromy methods and path continuation techniques.

## Examples

```julia

using PolynomialCovers

P = PolynomialCover([1.0, 0.0, -1.0, 0.0, 1.0])   # z^4 - z^2 + 1

F = fiber(P)

transport_to!(F, 5)

println(F.xs)



# P(z) = z² - 1
P = PolynomialCover([1.0, 0.0, -1.0])

y = 3.0

roots = cyclic_lift_solver(P, y)

println("Solutions of P(x) = $y")
println(roots)



```

To compare

```julia

using PolynomialCovers
using Polynomials

P = PolynomialCover([1.0, 0.0, -1.0, 0.0])   # z³ - z

y = 2.0

roots_lift = cyclic_lift_solver(P, y)

# Compare with classical polynomial root solver
poly = Polynomial([ -y, 0, -1, 1 ])  # z³ - z - y
roots_classical = roots(poly)

println("Cyclic lift solver:")
println(roots_lift)

println("\nClassical roots:")
println(roots_classical)
```

## Reference

The package accompanies the article

Root Finding via Polynomial Monodromy at Infinity:
Algorithms and Numerical Perspectives.
---

## Contributing

Contributions, bug reports, and suggestions are welcome.

The source code is available on GitHub.
