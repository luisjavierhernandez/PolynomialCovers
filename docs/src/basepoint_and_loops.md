# Basepoint and Initial Loops

Many algorithms for polynomial covers rely on lifting paths in the
codomain of the map
```math
P : \mathbb{C}\to \mathbb{C}
```
while avoiding the **critical values** of the polynomial.
In order to perform such lifts safely, one must first choose a
basepoint outside the region containing all critical values and
construct a loop based at this point.

This section describes how `PolynomialCovers.jl` constructs:

1. a safe **basepoint** in the codomain,
2. an **initial polygonal loop** around the origin,
3. a mutable state structure used during lifting.

The construction relies on the exclusion disk computed in the
*Critical Disks* section.

---

## Basepoint Selection

Given a polynomial cover `P`, we first compute a radius `r`
guaranteed to lie outside the disk containing all critical values.
This radius is obtained from the function

critical_values_exclusion_radius(P)

which inflates the critical value disk by a safety factor `s ≥ 3`.

The basepoint `r_0 = r` is then chosen on the positive real axis.

This choice simplifies path construction while guaranteeing that
the basepoint `r_0` and its image `P(r_0)` lie in the regular value set of the map.

## Initial Polygonal Loop

Once a basepoint `r_0` in the domain is chosen and its image `P(r_0)` is computed, we construct a counter-clockwise polygonal loop
centered at the origin and starting at w₀.

The number of edges nl ≥ 3 determines the discretization of the loop.
Because the radius is larger than the exclusion radius,
the entire loop stays outside the critical value disk.

Such loops are used as initial paths for:
- monodromy computations
- cyclic lifting solvers
- fiber exploration algorithms

## API
```@docs
PointState
```
## Geometric Interpretation

The overall construction can be summarized as follows:
	1.	Compute a disk containing all critical values.
	2.	Enlarge this disk by a safety factor.
	3.	Choose a basepoint outside the enlarged disk.
	4.	Construct a loop around the origin based at that point.
	5.	Lift the loop through the covering map.

This guarantees that the lifted paths remain in the regular region
of the polynomial cover, avoiding numerical instability near
critical fibers.
