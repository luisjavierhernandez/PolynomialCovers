

using LinearAlgebra
include("Geometry.jl")
include("PolynomialUtils.jl")










# ============================================================
# No export
    #=
    critical_baricenter,
    critical_radius_baricenter,
    critical_points_disk,
    candidate_disks
    =#
# ============================================================

# ============================================================
# Low-level coefficient bounds (your original ideas)
# ============================================================

"""
    cb(coeffs)

Cauchy bound for roots of a polynomial with leading coefficient `a₀`.
"""
cb(coeffs::AbstractVector{<:Number})::Real = begin
    @assert !isempty(coeffs) "Coefficient vector must be non-empty."
    a0 = coeffs[1]
    @assert a0 != 0 "Leading coefficient must be nonzero."

    if length(coeffs) == 1
        return one(real(a0))
    end

    ratios = abs.(coeffs[2:end] ./ a0)
    return one(real(a0)) + real(maximum(ratios))
end


"""
    cbprime(coeffs)

Cauchy bound for roots of the derivative.
"""
cbprime(coeffs::AbstractVector{<:Number})::Real = begin
    @assert length(coeffs) > 1 "Need degree ≥ 1."
    a0 = coeffs[1]
    @assert a0 != 0 "Leading coefficient must be nonzero."

    n = length(coeffs) - 1
    if n == 1
        return one(real(a0))
    end

    ratios = [
        abs((n - k) * coeffs[k+1] / (n * a0))
        for k in 1:(n-1)
    ]
    return one(real(a0)) + real(maximum(ratios))
end


"""
    fb(coeffs)

Fujiwara bound for roots.
"""
fb(coeffs::AbstractVector{<:Number})::Real = begin
    @assert length(coeffs) > 1 "Need degree ≥ 1."
    a0 = coeffs[1]
    @assert a0 != 0 "Leading coefficient must be nonzero."

    n = length(coeffs) - 1
    values = [
        abs(coeffs[k+1] / a0)^(1/k)
        for k in 1:n
    ]
    return 2 * one(real(a0)) * real(maximum(values))
end


"""
    fbprime(coeffs)

Fujiwara bound for roots of the derivative.
"""
fbprime(coeffs::AbstractVector{<:Number})::Real = begin
    @assert length(coeffs) > 1 "Need degree ≥ 1."
    a0 = coeffs[1]
    @assert a0 != 0 "Leading coefficient must be nonzero."

    n = length(coeffs) - 1
    if n == 1
        return 2 * one(real(a0))
    end

    values = [
        abs((n - k) * coeffs[k+1] / (n * a0))^(1/k)
        for k in 1:(n-1)
    ]
    m = real(maximum(values))
    if m == 0
        return 2 * one(real(a0))
    else
        return 2 * one(real(a0)) * m
    end
end


"""
    bounds_radius(coeffs)

Return all bounds as a named tuple.
"""
function bounds_radius(coeffs::AbstractVector{<:Number})
    return (
        cb      = cb(coeffs),
        fb      = fb(coeffs),
        cbprime = cbprime(coeffs),
        fbprime = fbprime(coeffs),
    )
end


"""
    min_max_critical_radius(coeffs)

Return `(rmin, rmax)` among all bounds.
"""
function min_critical_radius(coeffs::AbstractVector{<:Number})
    b = bounds_radius(coeffs)
    vals = (b.cb, b.fb, b.cbprime, b.fbprime)
    return minimum(vals)
end

# ============================================================
# Geometry of critical points
# ============================================================
#=
    critical_baricenter(P)

Baricenter of the critical points of a polynomial.
"""
function critical_baricenter(P::PolynomialCover)
    n = degree(P)
    if n ≤ 1
        return zero(eltype(P.coeffs))
    end

    a0 = P.coeffs[1]
    a1 = P.coeffs[2]

    return -a1 / (n * a0)
end
=#

"""
    critical_radius_origin(P)

Radius of a disk centered at 0 containing all critical points.
"""
function critical_points_radius_origin(P::PolynomialCover)
    coeffs = P.coeffs
    rmin = min_critical_radius(coeffs)
    return rmin
end
##
function critical_points_disk_origin(P::PolynomialCover)
    c = 0.0+0.0im
    r = critical_points_radius_origin(P)
    return Disk(c, r)
end


##

"""
    bound_values_radius_origin(P::PolynomialCover, R::Real)

Upper bound for |P(z)| on the disk |z| ≤ R.
"""
function bound_values_radius_origin(P::PolynomialCover, R::Real)
    coeffs = coefficients(P)
    n = length(coeffs) - 1
    a0 = coeffs[1]
    a0 != 0 || throw(ArgumentError("Leading coefficient must be nonzero"))

    S = zero(real(R))
    for k in 0:n
        S += abs(coeffs[k+1]) * R^(n - k)
    end
    return S
end

function critical_values_radius_origin(P::PolynomialCover)
    R = critical_points_radius_origin(P)
    bound_values_radius_origin(P, R)
end

function critical_values_disk_origin(P::PolynomialCover)
    r = critical_values_radius_origin(P)
    Disk(zero(r), r)
end

# ============================================================
# Exclusion radius for critical values
# ============================================================



# Internal helper: works on coefficients and a target radius S
function _min_radius(coeffs::AbstractVector{<:Number},
                     S::Real;
                     tol=1e-6,
                     max_iter=1000)

    a0 = coeffs[1]
    a0 != 0 || throw(ArgumentError("Leading coefficient must be nonzero"))

    n = length(coeffs) - 1

    function f(r)
        T = promote_type(typeof(S), typeof(r))
        tail = zero(T)
        for k in 1:n
            tail += abs(coeffs[k+1]) * r^(n-k)
        end
        abs(a0) * r^n - (S + tail)
    end

    low, high = zero(S), one(S)

    while f(high) ≤ 0
        high *= 2
        if high > 1e6
            return (Inf, :no_bracket)
        end
    end

    for _ in 1:max_iter
        mid = (low + high) / 2
        if f(mid) > 0
            high = mid
        else
            low = mid
        end
        if abs(high - low) < tol
            return (high, :ok)
        end
    end

    return (high, :max_iter)
end

"""
    critical_values_exclusion_radius(P; s=3)

Return a radius `r` such that any regular polygonal loop
centered at the origin with vertices on `|z| = r` and
at least 3 sides lies entirely outside the disk of
critical values.

The parameter `s > 1` inflates the critical value disk
to avoid tangencies. The default `s = 3` guarantees that
even triangular paths stay strictly away from critical
values, ensuring safe lifting.
"""


function critical_values_exclusion_radius(
    P::PolynomialCover;
    s::Real = 3,
    tol::Real = 1e-6,
    max_iter::Int = 1000
)
    s > 1 || throw(ArgumentError("s must be > 1"))

    # Disco de valores críticos (centrado en 0)
    Dcrit = critical_values_disk_origin(P)
    R = Dcrit.radius

    # Inflado del disco para evitar tangencias
    S = s * R

    r, status = _min_radius(
        coefficients(P),
        S;
        tol = tol,
        max_iter = max_iter
    )

    status == :ok ||
        @warn "Exclusion radius computation returned status $status"

    return r
end


function critical_values_exclusion_disk(
    P::PolynomialCover;
    s::Real = 3,
    tol::Real = 1e-6,
    max_iter::Int = 1000
)
    r = critical_values_exclusion_radius(
        P;
        s = s,
        tol = tol,
        max_iter = max_iter
    )
    return Disk(0 + 0im, r)
end



