# Critical Disks

This section describes the geometric bounds used to control the
location of critical points and critical values of a polynomial cover.

These bounds are used to construct safe regions in the complex plane
for numerical path lifting.

---

## Root Bounds

Several classical bounds are implemented to estimate the location of
the roots of a polynomial and of its derivative.

```@docs
PolynomialCovers.cb
PolynomialCovers.cbprime
PolynomialCovers.fb
PolynomialCovers.fbprime
PolynomialCovers.bounds_radius
PolynomialCovers.min_critical_radius
```
## Critical Point Bounds

The following functions estimate a disk centered at the origin
containing all critical points of a polynomial cover.

```@docs
PolynomialCovers.critical_points_radius_origin
```
## Critical Value Bounds

Using the bound on critical points we can derive bounds for the
critical values.

```@docs
critical_values_radius_origin
critical_values_disk_origin
```
## Exclusion Disk

For numerical path lifting we require loops that lie strictly outside
the region containing critical values.

The following routines compute a radius guaranteeing this property.


# Mathematical background

This section contains a description of the mathematical results used to construsct disks containing all the critical points and all the critical values.

## Critical Disks

Many algorithms in `PolynomialCovers` require constructing paths in the
complex plane that avoid the critical values of a polynomial map.  
To guarantee numerical stability, the package uses explicit bounds to
enclose the critical points and critical values in disks centered at
the origin.

This section explains the mathematical bounds used to construct those
regions.

---

## Polynomial Covers

A polynomial cover is a map

```math
P : \mathbb{C} \to \mathbb{C}
```
given by

```math
P(z) = a_0 z^n + a_1 z^{n-1} + \cdots + a_n ,
```
with (a_0 \neq 0).

The critical points of the map are the solutions of

```math
P'(z) = 0 ,
```

and the critical values are their images
```math
w_i = P(z_i).
```
The set of critical values plays a central role in the topology of the
covering map: path lifting is well-defined precisely on paths that
avoid this set.

---

## Cauchy Bound

If 
```math
P(z) = a_0 z^n + a_1 z^{n-1} + \cdots + a_n ,
```

then every root satisfies

```math
|z| \le 1 + \max_{1 \le k \le n} \left|\frac{a_k}{a_0}\right|.
```
This gives a disk centered at the origin containing all roots.

The same bound can be applied to the derivative (P’(z)) to enclose
all critical points.

## Fujiwara Bound

A sharper estimate is the Fujiwara bound
```math
|z| \le
2 \max_{1 \le k \le n}
\left|\frac{a_k}{a_0}\right|^{1/k}.
```
This bound is often tighter and therefore produces smaller enclosing
disks.

## Disk Containing Critical Points

Applying the above bounds to the derivative polynomial
```math
P'(z)
```
provides an explicit radius

```math
BF(P')
```
such that all critical points satisfy
```math
|z| \le BF(P').
```
This defines a disk
```math
\mathbb{D}_c = \{ z \in \mathbb{C} : |z| \le BF(P') \}.
```

## Disk Containing Critical Values

Once a bound on the critical points is known, we can bound the critical
values using the triangle inequality.

If
```math
|z| \le R_c ,
```
then
```math
|P(z)|
\le
\sum_{k=0}^{n} |a_k| R_c^{\,n-k}.
```
This produces a radius
```math
S_v
```
such that all critical values lie in the disk

```math
\mathbb{D}_v = \{ w \in \mathbb{C} : |w| \le R_v \}.
```
## Exclusion Radius for Path Lifting

When computing monodromy numerically, paths must remain strictly
outside the disk containing critical values.

To ensure this, the disk (D_v) is inflated by a factor (s>1):

```math
S = s S_v .
```
We then compute the smallest radius (r) such that
```math
|P(z)| > S
\quad \text{for all } |z| = r .
```
This guarantees that the circle of radius (r) in the base space lies
entirely outside the inflated critical value disk.

Polygonal loops constructed on this circle can therefore be lifted
safely without encountering singular fibers.

---

## Role in the Algorithms

These bounds serve several purposes in the algorithms implemented in
PolynomialCovers:
- providing a safe basepoint outside the critical value region,
- constructing loops for monodromy computation, and
- guaranteeing numerical stability of path lifting.

Although these bounds are not optimal, they are inexpensive to compute
and sufficiently sharp for robust numerical algorithms.

