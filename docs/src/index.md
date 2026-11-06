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
- Geometry
Geometric primitives in the complex plane such as disks,
vertices, and polygonal paths.
- Polynomial Covers
Representation of polynomial maps and their critical geometry.
- Lifting Algorithms
Numerical routines for lifting paths through polynomial maps.
- Solvers
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

## Citing this package

This package accompanies the article:

> **V. Álvarez-Aparicio** [![ORCID](https://img.shields.io/badge/ORCID-0000--0001--8232--1809-green?logo=orcid)](https://orcid.org/0000-0001-8232-1809),
> **L. J. Hernández-Paricio** [![ORCID](https://img.shields.io/badge/ORCID-0000--0003--4528--7781-green?logo=orcid)](https://orcid.org/0000-0003-4528-7781),
> **M. T. Rivas-Rodríguez** [![ORCID](https://img.shields.io/badge/ORCID-0000--0001--8911--4941-green?logo=orcid)](https://orcid.org/0000-0001-8911-4941).
> *On a polynomial root-finding algorithm based on the monodromy at infinity
> for equations with large independent term*.
> Department of Mathematics and Computing, University of La Rioja, 2026.

### Abstract

We investigate the monodromy structure of complex polynomials at infinity
from a computational viewpoint and introduce a geometric algorithm for
solving polynomial equations of the form $P(x)=b$ when $b$ lies outside
the disk containing the critical values. The proposed method exploits the
cyclic monodromy induced by loops around infinity to generate the entire
fiber from a single initial lift via numerical path lifting.
A pseudocode formulation of the algorithm is provided together with a
comparison with classical root-finding techniques, including companion matrix
eigenvalue methods and homotopy continuation approaches. Numerical experiments
indicate that, for the considered family of sparse polynomials, the method
achieves competitive accuracy while naturally preserving the geometric ordering
of the roots.

**Keywords:** Polynomial coverings · Polynomial root-finding · Monodromy ·
Path-lifting · Numerical continuation · Computational algebraic geometry ·
Homotopy methods

### BibTeX

```bibtex
@article{AlvarezAparicio2026polynomial,
  author    = {{\'{A}}lvarez-Aparicio, V. and
               Hern{\'{a}}ndez-Paricio, L.~J. and
               Rivas-Rodr{\'{i}}guez, M.~T.},
  title     = {On a polynomial root-finding algorithm based on the monodromy
               at infinity for equations with large independent term},
  year      = {2026},
  note      = {Preprint, University of La Rioja},
  url       = {},
}
```

> Once the article is assigned a DOI, replace `note` and `url` with
> `journal`, `volume`, `pages`, and `doi`.

### Funding

This research has been funded by the University of La Rioja, grant REGI25/61.

---

## Contributing

Contributions, bug reports, and suggestions are welcome.

The source code is available on GitHub.
