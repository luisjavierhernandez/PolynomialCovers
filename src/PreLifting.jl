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

### PathState

"""
    PathState{T}

State of a lifted path used for visualization and tracing.

This structure stores the current point of the lifting process
together with the trajectory followed during path lifting.

The invariant relation

    w = P(x)

holds at every accepted step of the algorithm.

Fields
-------
- `x` : current point in the domain
- `w` : current image in the codomain
- `pathx` : recorded trajectory in the domain
- `pathy` : recorded trajectory in the codomain

Notes
-----
This structure is mainly intended for **didactic and visualization
purposes**, since it stores the full lifting trajectory.
For high-performance lifting without recording intermediate
points see the `PointState` based routines.
"""
mutable struct PathState{T<:Real}
    x::Complex{T}
    w::Complex{T}
    pathx::Vector{Complex{T}}
    pathy::Vector{Complex{T}}
end

function PathState(x::Complex{T1}, w::Complex{T2}) where {T1<:Real,T2<:Real}
    T = promote_type(T1,T2)
    PathState{T}(Complex{T}(x), Complex{T}(w), [Complex{T}(x)], [Complex{T}(w)])
end

PathState(x::Complex{T}, w::Complex{T}) where T =
    PathState{T}(x, w, [x], [w])


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

