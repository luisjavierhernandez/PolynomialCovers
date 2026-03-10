"""
Geometry primitives used in PolynomialCovers.

This document defines disks, vertices, and polygonal paths in the
complex plane together with basic geometric utilities used in
path lifting algorithms.
"""

import Base: in

"""
    Disk{T}

Closed disk in the complex plane.

# Fields
- `center::Complex{T}` — center of the disk
- `radius::T` — non–negative radius

# Notes
The disk represents the closed set

    { z ∈ ℂ : |z - center| ≤ radius }.
"""
struct Disk{T<:Real}
    center :: Complex{T}
    radius :: T

    function Disk{T}(c::Complex{T}, r::T) where {T<:Real}
        r ≥ zero(T) || throw(ArgumentError("Radius must be non-negative"))
        new{T}(c, r)
    end
end

# ------------------------------------------------------------
# Constructors
# ------------------------------------------------------------

Disk(c::Real, r::Real) = Disk(complex(c, zero(c)), r)

Disk(c::Complex, r::Real) = begin
    T = promote_type(typeof(real(c)), typeof(r))
    Disk{T}(Complex{T}(c), T(r))
end

# ------------------------------------------------------------
# Membership
# ------------------------------------------------------------
"""
    in(z, d::Disk)

Return `true` if the complex number `z` lies in the closed disk `d`.
"""
in(z::Complex, d::Disk) =
    abs(z - d.center) ≤ d.radius

"""
    onboundary(d::Disk, z::Complex)

Return `true` if the complex number `z` lies on the boundary of the disk `d`.
"""
onboundary(d::Disk, z::Complex) =
    abs(abs(z - d.center) - d.radius) ≤ eps(d.radius)

# Extend Base.contains for Disk
Base.contains(d::Disk, z::Complex) = in(z, d)

# ------------------------------------------------------------
# Distance
# ------------------------------------------------------------
"""
    distance(d::Disk, z)

Distance from the point `z` to the disk `d`.

Returns `0` if `z` lies inside the disk.
"""
function distance(d::Disk, u::Complex)
    δ = abs(u - d.center) - d.radius
    max(zero(d.radius), δ)
end

# ------------------------------------------------------------
# Transformations
# ------------------------------------------------------------
"""
    translate(d::Disk, w)

Translate the disk `d` by the complex displacement `w`.
"""
translate(d::Disk{T}, w::Complex) where {T} =
    Disk{T}(d.center + Complex{T}(w), d.radius)

"""
    scale(d::Disk, λ)

Scale the radius of the disk `d` by the factor `λ ≥ 0`.
The center remains unchanged.
"""

function scale(d::Disk{T}, λ::Real) where {T}
    λ ≥ zero(λ) || throw(ArgumentError("Scale must be non-negative"))
    S = promote_type(T, typeof(λ))
    Disk{S}(Complex{S}(d.center), S(d.radius * λ))
end




#
# Discrete oriented polygonal paths in the complex plane.
#
# This file provides:
# - vertices (`Vertex`) as points in ℂ
# - discrete paths (`Path`) as sequences of vertices
#
# The design is minimal and purely geometric.


# ============================================================
# Vertex
# ============================================================

"""
    Vertex{T}

Vertex (point) in the complex plane.

# Field
- `z::Complex{T}` — complex coordinate of the vertex

Vertices are used to build polygonal paths (`Path`).
"""


struct Vertex{T<:Real}
    z::Complex{T}
end

Vertex(x::T) where {T<:Real} =
    Vertex(complex(x, zero(x)))

Vertex(x::T, y::T) where {T<:Real} =
    Vertex(complex(x, y))

Base.show(io::IO, v::Vertex) = print(io, v.z)

Base.:-(v::Vertex, w::Vertex) = v.z - w.z
Base.:+(v::Vertex, w::Vertex) = v.z + w.z
Base.abs(v::Vertex) = abs(v.z)
import Base: ==
==(v::Vertex, w::Vertex) = v.z == w.z

# ============================================================
# Path
# ============================================================

"""
    Path{T}

Discrete polygonal path represented as a sequence of vertices.

# Field
- `vertices::Vector{Vertex{T}}`

A path with `n` vertices represents `n-1` linear segments.
"""
struct Path{T<:Real}
    vertices::Vector{Vertex{T}}

    function Path{T}(verts::Vector{Vertex{T}}) where {T<:Real}
        length(verts) ≥ 1 ||
            throw(ArgumentError("A path must have at least one vertex"))
        new(verts)
    end
end

# ------------------------------------------------------------
# Constructors
# ------------------------------------------------------------

Path(verts::Vector{Vertex{T}}) where {T<:Real} =
    Path{T}(verts)

Path(points::AbstractVector{<:Complex{T}}) where {T<:Real} =
    Path(Vertex.(points))

Path(points::AbstractVector{<:Real}) =
    Path(complex.(points))

function Base.show(io::IO, p::Path)
    print(io, "Path(")
    for (i, v) in enumerate(p.vertices)
        i > 1 && print(io, " → ")
        print(io, v)
    end
    print(io, ")")
end
# ============================================================
# Basic access
# ============================================================

Base.length(p::Path) = length(p.vertices)

"""
    startpoint(p)

Return the first point of the path `p`.
"""

startpoint(p::Path) = p.vertices[1].z

"""
    endpoint(p)

Return the final point of the path `p`.
"""
endpoint(p::Path)   = p.vertices[end].z



# ============================================================
# Concatenation
# ============================================================
"""
    promote_vertices(z0, zf, S)

Promote vertices to a common floating point type.

This ensures that path construction algorithms operate with a
consistent numeric type.
"""
function concatenate(p::Path{T}, q::Path{T}) where {T}
    endpoint(p) == startpoint(q) ||
        error("Paths are not compatible for concatenation")

    Path(vcat(p.vertices, q.vertices[2:end]))
end


# ============================================================
# Polygonal loops
# ============================================================

"""
    polygonal_loop_origin(z, l)

Construct a regular polygonal loop with:
- center at the origin,
- `l ≥ 3` sides,
- counterclockwise (levogiro) orientation,
- starting and ending at `z`.

Returns a `Path`.
"""
function polygonal_loop_origin(z::Complex{T}, l::Int) where {T<:Real}
    l ≥ 3 || throw(ArgumentError("l must be ≥ 3"))

    ω = cis(2π / l)  # primitive l-th root of unity (CCW)
    verts = Vertex{T}[]

    zk = z
    for _ in 1:l
        push!(verts, Vertex(zk))
        zk *= ω
    end

    push!(verts, Vertex(z))  # close the loop
    Path(verts)
end




"""
    segment_intersects_disk(z1, z2, d::Disk)

Return `true` if the closed segment from `z1` to `z2`
intersects the disk `d`.
"""
function segment_intersects_disk(z1::Complex, z2::Complex, d::Disk)
    R = d.radius

    a = z1 - d.center
    b = z2 - z1

    # Degenerate segment (single point)
    if b == 0
        return abs(a) ≤ R
    end

    # Projection parameter onto the line
    t = -real(conj(b) * a) / abs2(b)
    t = clamp(t, zero(t), one(t))

    # Closest point to the center
    zmin = a + t * b

    return abs(zmin) ≤ R
end

# Convenience wrapper: allow passing a numeric radius (disk centered at origin)
function segment_intersects_disk(z1::Complex, z2::Complex, S::Real)
    d = Disk(0.0 + 0.0im, S)
    return segment_intersects_disk(z1, z2, d)
end
"""
    path_intersects_disk(p::Path, S)

Return `true` if at least one segment of the path `p`
intersects the disk centered at the origin with radius `3S/2`.
"""
function path_intersects_disk(p::Path, S::Real)
    verts = p.vertices

    for i in 1:length(verts)-1
        z1 = verts[i].z
        z2 = verts[i+1].z

        if segment_intersects_disk(z1, z2, S)
            return true
        end
    end

    return false
end

"""
    path_avoids_disk(p::Path, S)

Return `true` if the entire path avoids the disk of radius `3S/2`
centered at the origin.
"""
path_avoids_disk(p::Path, S::Real) =
    !path_intersects_disk(p, S)

# Path-level check for an arbitrary Disk (no shifting)
function path_intersects_disk(p::Path, d::Disk)
    verts = p.vertices

    for i in 1:length(verts)-1
        z1 = verts[i].z
        z2 = verts[i+1].z

        if segment_intersects_disk(z1, z2, d)
            return true
        end
    end

    return false
end

# Convenience wrapper: check whether a path avoids a given Disk (arbitrary center)
function avoids(p::Path, d::Disk)
    !path_intersects_disk(p, d)
end


"""
    safe_path_arc(u, v, ρ; npoints=6)

Construct a polygonal path from `u` to `v` that avoids the disk
of radius `ρ` centered at the origin.

If the straight segment does not intersect the disk, the path
is the segment `[u, v]`. Otherwise, the path detours along a
circular arc of radius slightly larger than `ρ`.

Returns a `Path`.
"""
function safe_path_arc(
    u::Complex{T},
    v::Complex{T},
    rho::Real;
    npoints::Int = 6
) where {T<:Real}

    abs(u) > rho || throw(ArgumentError("u must lie outside the disk"))
    abs(v) > rho || throw(ArgumentError("v must lie outside the disk"))

    # Check if straight segment intersects the disk
    if !segment_intersects_disk(u, v, rho)
        return Path([Vertex(u), Vertex(v)])
    end

    # Angular positions
    θu = angle(u)
    θv = angle(v)

    Δθ = θv - θu
    if Δθ > π
        Δθ -= 2π
    elseif Δθ < -π
        Δθ += 2π
    end

    # Build arc slightly outside the forbidden disk
    r = T(1.05 * rho)
    verts = Vertex{T}[]

    push!(verts, Vertex(u))

    for j in 0:npoints
        θ = θu + j * Δθ / npoints
        z = r * cis(θ)
        push!(verts, Vertex(z))
    end

    push!(verts, Vertex(v))

    return Path(verts)
end
###

"""
    promote_vertices(z0, zf, S)

Promote vertices to a common floating point type.

This ensures that path construction algorithms operate with a
consistent numeric type.
"""

function promote_vertices(
    z0::Vertex{T1},
    zf::Vertex{T2},
    S::Vector{Vertex{T3}},
) where {T1,T2,T3}

    T = promote_type(T1,T2,T3)

    convertV(v) = Vertex{T}(Complex{T}(v.z))

    return convertV(z0),
           convertV(zf),
           convertV.(S)
end


###

"""
    avoid_points_path(z0, zf, S; kwargs...)

Construct a polygonal path from `z0` to `zf` that avoids a set
of forbidden vertices `S`.

The algorithm uses a potential field method combining attraction
towards the target and repulsion from forbidden points.
"""

function avoid_points_path(
    z0::Vertex{T1},
    zf::Vertex{T2},
    S::Vector{Vertex{T3}};
    d_safe = 1.0,
    step   = 0.05,
    λ      = 0.3,
    step_min = 0.15,
    step_max = 0.3,
) where {T1,T2,T3}

    z0T, zfT, ST = promote_vertices(z0, zf, S)

    T = eltype(real(z0T.z))  # ← tipo final REAL

    return _avoid_points_path(
        z0T, zfT, ST;
        d_safe = T(d_safe),
        step   = T(step),
        λ      = T(λ),
        step_min = T(step_min),
        step_max = T(step_max),
    )
end

"""
Internal implementation of `avoid_points_path`.
"""

function _avoid_points_path(
    z0::Vertex{T},
    zf::Vertex{T},
    S::Vector{Vertex{T}};
    step::T = T(0.05),
    λ::T = T(0.3),
    d_safe::T = T(1.0),
    step_min::T = T(0.15),
    step_max::T = T(1.0),
    α::T = T(0.3),
    tol::T = T(1e-6),
    maxiter::Int = 10_000,
) where {T<:AbstractFloat}

    z = z0.z
    target = zf.z

    final_is_forbidden = any(v -> v == zf, S)

    verts = Vertex{T}[z0]

    Vsmooth = zero(z)

    orient = sign(imag(target - z0.z) + eps(T))

    for _ in 1:maxiter

        # ---------------------------
        # Attraction
        # ---------------------------
        V = target - z

        # ---------------------------
        # Repulsion + tangential bypass
        # ---------------------------
        mindist = typemax(T)

        for a in S

            if final_is_forbidden && a == zf
                continue
            end

            d = z - a.z
            r = abs(d)

            mindist = min(mindist, r)

            r < eps(T) &&
                error("Collision with forbidden point")

            if r < d_safe
                radial = d / r

                w = λ * (1/r - 1/d_safe) / r
                V += w * radial

                tangential = orient * im * radial
                V += T(0.5) * w * tangential
            end
        end

        abs(V) < sqrt(eps(T)) && break

        # ---------------------------
        # Curvature smoothing
        # ---------------------------
        if Vsmooth == 0
            Vsmooth = V
        else
            Vsmooth = (one(T)-α)*Vsmooth + α*V
        end

        # ---------------------------
        # Adaptive step
        # ---------------------------
        scale = clamp(
            mindist / d_safe,
            step_min,
            step_max
        )

        step_eff = step * scale

        dir = Vsmooth / abs(Vsmooth)

        # forbid backward motion
        forward = target - z
        if real(conj(dir) * forward) ≤ 0
            dir = forward / abs(forward)
        end

        z_new = z + step_eff * dir

        # ---------------------------
        # Reject cusp creation
        # ---------------------------
        if length(verts) ≥ 2
            z_prev  = verts[end].z
            z_prev2 = verts[end-1].z

            v1 = z_prev - z_prev2
            v2 = z_new - z_prev

            if real(conj(v1)*v2) ≤ 0
                continue
            end
        end

        z = z_new
        push!(verts, Vertex{T}(z))

        abs(z - target) < tol && break
    end

    # ---------------------------
    # SAFE ENDPOINT INSERTION ⭐
    # ---------------------------
    if length(verts) ≥ 2
        z_prev  = verts[end].z
        z_prev2 = verts[end-1].z

        v1 = z_prev - z_prev2
        v2 = zf.z - z_prev

        if real(conj(v1)*v2) > 0
            push!(verts, zf)
        else
            verts[end] = zf
        end
    else
        push!(verts, zf)
    end

    return Path(verts)
end

###

"""
    approach_path(
    u::Complex{T},
    v::Complex{T},
    rho::Real;
    npoints::Int = 6
) where {T<:Real}

TBW
"""
function approach_path(
    u::Complex{T},
    v::Complex{T},
    rho::Real;
    npoints::Int = 6
) where {T<:Real}

    @assert abs(u) < rho "Initial point must satisfy |u| < rho"

    # Caso v = 0 → segmento radial
    if iszero(v)

        dir = u / abs(u)
        w = rho * dir

        pts = Complex{T}[u, w]
        return Path(pts)

    end

    # Caso general → arco seguro hasta la circunferencia
    w = rho * (v / abs(v))

    return safe_path_arc(u, w, rho; npoints=npoints)

end

"""
    min_pairwise_distance(S)

Return the minimal distance between vertices in `S`.
"""


function min_pairwise_distance(S::Vector{Vertex{T}}) where {T<:Real}

    n = length(S)
    n < 2 && return zero(T)

    dmin2 = typemax(T)

    @inbounds for i in 1:n-1
        zi = S[i].z
        for j in i+1:n
            d2 = abs2(zi - S[j].z)
            dmin2 = min(dmin2, d2)
        end
    end

    return sqrt(dmin2)
end





