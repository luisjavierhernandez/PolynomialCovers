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

$$
P(x) = y
$$

can be computed by lifting paths in the complex plane through the
polynomial map

$$
P : \mathbb{C} \to \mathbb{C}.
$$

In particular, the package implements numerical algorithms inspired by
monodromy methods and path continuation techniques.

## Reference

The package accompanies the article

Root Finding via Polynomial Monodromy at Infinity:
Algorithms and Numerical Perspectives.
---

## Contributing

Contributions, bug reports, and suggestions are welcome.

The source code is available on GitHub.
