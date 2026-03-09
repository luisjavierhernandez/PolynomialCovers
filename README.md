# PolynomialCovers.jl

`PolynomialCovers.jl` provides tools for working with **polynomial covering maps**
and performing **path lifting in the complex plane**.

The package implements algorithms for

- computing fibers of polynomial maps
- constructing paths in the codomain that avoid critical values
- lifting paths and transporting fibers
- experimenting with continuation and monodromy of polynomial covers
- solving equations of the form P(x)=b with for values of b sufficiently far 

These tools are intended for research and experimentation in
computational algebra, complex geometry, and numerical continuation.

---

## Installation

Clone the repository and add the package to your Julia environment:

```julia
using Pkg
Pkg.add(path="path/to/PolynomialCovers")

## Quick Start

The following example illustrates the basic workflow when working with
a polynomial covering map.

using PolynomialCovers

# Define a polynomial cover
P = PolynomialCover(1, 0, -1)     # P(z) = z^2 - 1

# Compute a safe basepoint in the codomain
b = basepoint(P)

# Compute the fiber P⁻¹(b)
F = fiber(P)

println("Fiber points:")
println(F)

# Construct a path avoiding critical values
γ = safe_path_arc(b, -b, 2.0)

# Transport the fiber along the path
transport!(F, γ)

println("Fiber after transport:")
println(F)
```

The typical workflow is:

Polynomial → basepoint → fiber → path → transport

## Example: Monodromy Around a Critical Value

A fundamental phenomenon in polynomial covers is that lifting a loop
around a critical value permutes the points of the fiber.

```julia

using PolynomialCovers

# Polynomial P(z) = z^2 - 1
P = PolynomialCover(1, 0, -1)

# Basepoint
b = basepoint(P)

# Fiber over the basepoint
F = fiber(P)

println("Initial fiber:")
println(F)

# Loop around a critical value
γ = basepoint_initial_loop(P)

# Lift the loop
transport!(F, γ)

println("Fiber after lifting the loop:")
println(F)
```

For this polynomial the loop exchanges the two roots,
demonstrating the monodromy permutation of the cover.

⸻

## Cyclic Lift Solver

For values of b sufficiently far from the critical values of the
polynomial, the equation

P(z) = b

an be solved efficiently using a cyclic lifting method.

The package provides the function cyclic_lift_solver, which is
particularly effective when

|b| ≥ 3S / 2

```julia

using PolynomialCovers

P = PolynomialCover(1, 0, -1)

# Radius containing the critical values
S = critical_values_radius(P)

println("Critical values lie in the disk |z| ≤ $S")

# Choose a value far from the critical values
b = 3S

# Solve P(z) = b
z = cyclic_lift_solver(P, b)

println("Solution:", z)
```

In this regime the cyclic lifting method can outperform standard
Newton-based solvers both in runtime and numerical stability.

⸻

## Main Concepts

The package revolves around the following core objects:

Object/Description
PolynomialCover/Polynomial covering map
Fiber/Fiber of the map over a point
Path/Path in the codomain
Disk/Disk used for geometric constructions

Core operations include
	•	fiber(P, b) – compute the fiber over b
	•	transport!(F, γ) – transport a fiber along a path
	•	safe_path_arc – construct paths avoiding critical values
	•	cyclic_lift_solver – solve P(z) = b using cyclic lifting


## Status

This package is currently under active development.

Planned features include
	•	improved lifting algorithms
	•	certified continuation methods
	•	monodromy computation tools
	•	additional numerical solvers for fibers

⸻

## License

MIT License

Copyright (c) 2026 Luis Javier Hernández Paricio

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.







