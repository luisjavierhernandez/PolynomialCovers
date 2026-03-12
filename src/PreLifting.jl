"""
    PointState{T}

State of a lifted point during path lifting.

Fields
-------
- `x`: current point in the domain
- `w`: current image in the codomain

This structure is mutated during lifting algorithms in order
to track the current position along the lifted path.
"""
mutable struct PointState{T<:Real}
    x::Complex{T}
    w::Complex{T}
end

#include("PolynomialCriticalDisks.jl")
#include("ApproximateLiftPath.jl")


# -------------------pr-----------------------------------------
# Base point
# ------------------------------------------------------------
"""
    basepoint_radius(P; s=3)

Return a radius guaranteeing that the circle |w| = R lies
outside the disk containing the critical values of `P`.

This radius is computed using `critical_values_exclusion_radius`.
The parameter `s ≥ 3` inflates the critical value disk to ensure
safe lifting paths.
"""
function basepoint_radius(
    P::PolynomialCover{T};
    s::Real = 3,
    tol::Real = 1e-6,
) where {T<:Real}
    s ≥ 3 || throw(ArgumentError("s must be ≥ 3"))

    return critical_values_exclusion_radius(P; s=s, tol=tol)
end


"""
    basepoint(P; s=3)

Return a basepoint in the codomain outside the critical value disk
of the polynomial cover `P`.

The basepoint is chosen on the positive real axis at radius
determined by `basepoint_radius`.
"""
function basepoint(P::PolynomialCover; s=3)
    R = basepoint_radius(P; s=s)
    return complex(R, 0)
end




"""
    basepoint_initial_loop(P; nl=8, s=3)

Construct an initial polygonal loop in the codomain of `P`.

Steps
1. choose a basepoint `z0` in the domain,
2. compute its image `w0 = P(z0)`,
3. construct a regular CCW polygonal loop centered at the origin
   starting and ending at `w0`.

Arguments
---------
- `nl ≥ 3`: number of polygon edges
- `s`: safety parameter used in the basepoint selection

Returns
-------
A `Path` in the codomain suitable for monodromy lifting.
"""
function basepoint_initial_loop(
    P::PolynomialCover;
    nl::Int = 8,
    s::Real = 3
)
    nl ≥ 3 || throw(ArgumentError("nl must be ≥ 3"))

    # 1. basepoint in the domain
    z0 = basepoint(P; s=s)

    # 2. image in the codomain
    w0 = value(P, z0)

    # 3. polygonal loop in the codomain
    return polygonal_loop_origin(w0, nl)
end

#

