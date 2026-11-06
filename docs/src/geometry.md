# Geometry

The `Geometry` module provides basic geometric structures and utilities
used throughout `PolynomialCovers`.

These tools operate in the complex plane and support the construction
of polygonal paths used by the lifting algorithms implemented in the
package.

The module introduces three fundamental objects:

- `Disk` — closed disks in the complex plane
- `Vertex` — points in the complex plane
- `Path` — discrete polygonal paths built from vertices

These objects are used to construct safe paths that avoid critical
value regions when performing numerical path lifting.

---

## Disks

A `Disk` represents a closed disk in the complex plane.

Mathematically it corresponds to

```math
\mathbb{D} = \{ z ∈ ℂ : |z - c| ≤ r \}
```
where `c` is the center and `r ≥ 0` the radius.

## Type

```@docs
Disk
```

## Membership and boundary
```@docs
PolynomialCovers.onboundary
```
## Distance

```@docs
distance
```
## Transformations

```@docs
PolynomialCovers.translate
```
## Vertices

Vertices represent points in the complex plane. They are used as the
building blocks of polygonal paths.

## Type

```@docs
Vertex
```
## Polygonal Paths

A Path represents a discrete polygonal curve in the complex plane.

It is defined by a sequence of vertices
```
v₁ → v₂ → ... → vₙ
```
## Type

```@docs
Path
```
## Concatenation

```@docs
PolynomialCovers.concatenate
```
## Polygonal Loops

Regular polygonal loops centered at the origin are frequently used in
monodromy computations.

```@docs
PolynomialCovers.polygonal_loop_origin
```

These loops approximate circular paths while remaining purely
polygonal, which simplifies the lifting algorithms.

---

## Disk Intersections

Many algorithms in PolynomialCovers require constructing paths that
avoid forbidden regions (typically disks containing critical values).

The following utilities detect intersections between segments or paths
and disks.

```@docs
PolynomialCovers.segment_intersects_disk
PolynomialCovers.path_intersects_disk
PolynomialCovers.path_avoids_disk
```
## Safe Path Construction

The module provides algorithms to construct paths that avoid forbidden
regions in the complex plane.

These routines are crucial for robust numerical lifting.

```@docs
safe_path_arc
PolynomialCovers.approach_path
avoid_points_path
```

	•	safe_path_arc constructs a path between two points that detours
around a forbidden disk if necessary.
	•	approach_path constructs a path approaching the boundary of a disk
from the inside.
	•	avoid_points_path constructs a path avoiding a finite set of
forbidden points.

Internally, the avoidance algorithm uses a potential-field strategy
combining attraction toward the target and repulsion from forbidden
points.

---

The geometric utilities (type promotion and pairwise distance estimation) are used internally to ensure numerical robustness.

---

## Example

The following example constructs a polygonal path.


```@example
using PolynomialCovers

v1 = Vertex(0.0)
v2 = Vertex(1.0)
p = Path([v1,v2])
println(p)
```



## Role in PolynomialCovers

The geometry module plays a central role in the implementation of the
lifting algorithms used by the package.

In particular, the path construction routines are used to build
numerically safe paths in the complex plane that avoid the critical
value disks of polynomial maps.