import Base: in

"""
    Disk{T}

Closed disk in ℂ with center ∈ Complex{T} and radius ∈ T.
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

in(z::Complex, d::Disk) =
    abs(z - d.center) ≤ d.radius

onboundary(d::Disk, z::Complex) =
    abs(abs(z - d.center) - d.radius) ≤ eps(d.radius)

# Extend Base.contains for Disk
Base.contains(d::Disk, z::Complex) = in(z, d)

# ------------------------------------------------------------
# Distance
# ------------------------------------------------------------

function distance(d::Disk, u::Complex)
    δ = abs(u - d.center) - d.radius
    max(zero(d.radius), δ)
end

# ------------------------------------------------------------
# Transformations
# ------------------------------------------------------------

translate(d::Disk{T}, w::Complex) where {T} =
    Disk{T}(d.center + Complex{T}(w), d.radius)

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

Vertex (point) in ℂ.
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

Discrete polygonal path represented as a sequence of vertices (can be a singleton).
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



startpoint(p::Path) = p.vertices[1].z
endpoint(p::Path)   = p.vertices[end].z



# ============================================================
# Concatenation
# ============================================================

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