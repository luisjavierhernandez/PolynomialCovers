# PolynomialCovers.jl

`PolynomialCovers.jl` provides tools for working with **polynomial covering maps**
and performing **path lifting in the complex plane**.

The package implements algorithms for

- computing fibers of polynomial maps
- constructing paths in the codomain that avoid critical values
- lifting paths and transporting fibers
- experimenting with continuation and monodromy of polynomial covers

These tools are intended for research and experimentation in
computational algebra, complex geometry, and numerical continuation.

---

## Installation

Clone the repository and add it to your Julia environment:

```julia
using Pkg
Pkg.add(path="path/to/PolynomialCovers")

## Quick start

using PolynomialCovers

# Define a polynomial cover
P = PolynomialCover([1, 0, -1])     # P(z) = z^2 - 1

# Compute a safe basepoint in the codomain
b = basepoint(P)

# Compute the fiber P⁻¹(b)
F = fiber(P)

println("Fiber points:")
println(F)

# Construct a path avoiding critical values
γ = safe_path_arc(b, -b)

# Transport the fiber along the path
transport!(F, γ)

println("Fiber after transport:")
println(F)

##Example: Monodromy Around a Critical Value

A key phenomenon in polynomial covers is that lifting a loop around
a critical value permutes the points of the fiber.

The following example illustrates this behavior.


using PolynomialCovers

# Polynomial P(z) = z^2 - 1
P = PolynomialCover([1,0,-1])

# Basepoint in the codomain
b = basepoint(P)

# Compute the fiber
F = fiber(P)

println("Initial fiber:")
println(F)

# Construct a loop around the origin (critical value)
γ = basepoint_initial_loop(P)

# Transport the fiber along the loop
transport!(F, γ)

println("Fiber after lifting the loop:")

For this polynomial the loop exchanges the two roots,
demonstrating the monodromy permutation of the cover.

⸻

Status

This package is currently under active development.

Planned features include:
	•	improved lifting algorithms
	•	certified continuation methods
	•	monodromy computation tools
	•	additional numerical solvers for fibers

⸻

License

MIT License