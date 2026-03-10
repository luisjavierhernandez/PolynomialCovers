
mutable struct PointState{T<:Real}
    x::Complex{T}          # punto actual en el dominio
    w::Complex{T}          # imagen actual
end

#include("PolynomialCriticalDisks.jl")
#include("ApproximateLiftPath.jl")


# -------------------pr-----------------------------------------
# Base point
# ------------------------------------------------------------

function basepoint_radius(
    P::PolynomialCover{T};
    s::Real = 3,
    tol::Real = 1e-6,
) where {T<:Real}
    s ≥ 3 || throw(ArgumentError("s must be ≥ 3"))

    return critical_values_exclusion_radius(P; s=s, tol=tol)
end


"""
    basepoint(P::PolynomialCover; s=3)

Return a finite basepoint outside the finite critical value locus of `P`,
intended as a starting point for lifting and monodromy constructions.
"""
function basepoint(P::PolynomialCover; s=3)
    R = basepoint_radius(P; s=s)
    return complex(R, 0)
end


"""
    basepoint_initial_loop(P::AbstractBranchedCover; nl=8, s=3)

Construct the initial polygonal loop in the codomain associated
to a branched cover `P`.

Steps:
1. choose a basepoint `z0 = basepoint(P; s=s)` in the domain,
2. compute its image `w0 = value(P, z0)`,
3. construct a regular CCW polygonal loop centered at the origin
   starting and ending at `w0`.

Arguments:
- `nl ≥ 3`: number of sides of the polygon
- `s`: safety parameter passed to `basepoint`

Returns:
- `Path` in the codomain
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

