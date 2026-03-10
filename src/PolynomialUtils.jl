
"""
    PolynomialCover{T}

Representation of a polynomial map

    P : ℂ → ℂ

given by its coefficients

    P(z) = a₀ zⁿ + a₁ zⁿ⁻¹ + ⋯ + aₙ.

The coefficients are stored in descending order.

# Fields

- `coeffs`: vector of complex coefficients.
"""
struct PolynomialCover{T<:Real}
    coeffs::Vector{Complex{T}}
end


"""
    coefficients(P)

Return the coefficient vector of the polynomial cover.
"""
coefficients(P::PolynomialCover) = P.coeffs


"""
    degree(P)

Return the degree of the polynomial.
"""
degree(P::PolynomialCover) = length(P.coeffs) - 1


"""
    PolynomialCover(coeffs)

Construct a polynomial cover from a vector of coefficients.
"""
function PolynomialCover(coeffs::AbstractVector{<:Number})
    T = float(real(eltype(coeffs)))
    PolynomialCover{T}(Complex{T}.(coeffs))
end

PolynomialCover(coeffs::Number...) = PolynomialCover(collect(coeffs))


"""
    value(P, z)

Evaluate the polynomial `P` at `z`.
"""
function value(P::PolynomialCover, z::Number)

    T = promote_type(eltype(P.coeffs), typeof(z))
    acc = zero(T)
    zT = T(z)

    for a in P.coeffs
        acc = acc*zT + T(a)
    end

    return acc
end


"""
    P(z)

Evaluate the polynomial cover at `z`.
"""
function (P::PolynomialCover)(z)
    value(P, z)
end


"""
    derivative(P, z)

Evaluate the derivative `P'(z)`.
"""
function derivative(P::PolynomialCover, z::Number)

    n = degree(P)

    T = promote_type(eltype(P.coeffs), typeof(z))
    acc = zero(T)
    zT = T(z)

    for i in 1:n
        acc = acc*zT + T(P.coeffs[i])*(n - i + 1)
    end

    return acc
end


"""
    critical_points(P)

Return the critical points of the polynomial map.
"""
function critical_points(P::PolynomialCover{T}) where T

    n = degree(P)
    n == 0 && return Complex{T}[]

    deriv_coeffs =
        Complex{T}[(P.coeffs[i]*(n - i + 1)) for i in 1:n]

    deriv_poly = PolynomialCover{T}(deriv_coeffs)

    return roots(deriv_poly)
end


"""
    critical_values(P)

Return the critical values `P(c)` where `c` is a critical point.
"""
function critical_values(P::PolynomialCover{T}) where T

    cps = critical_points(P)

    return [value(P, cp) for cp in cps]
end


"""
    fiber_via_companion(P, y)

Compute the solutions of

    P(x) = y

using the companion matrix method.
"""
function fiber_via_companion(
    P::PolynomialCover{T},
    y0::Number
) where T

    n = degree(P)
    n == 0 && return Complex{float(T)}[]

    Tf = float(T)
    CT = Complex{Tf}

    y0c = CT(y0)

    a = CT.(P.coeffs)
    a[end] -= y0c

    C = zeros(CT, n, n)

    @inbounds for i in 2:n
        C[i, i-1] = one(CT)
    end

    @inbounds for j in 1:n
        C[1, j] = -a[j+1] / a[1]
    end

    return Complex{T}.(eigvals(C))
end
#=
# ------------------------------------------------------------
# Tipo concreto
# ------------------------------------------------------------

"""
PolynomialCover(coeffs)

Representation of a polynomial map

    P : ℂ → ℂ

given by its coefficients in descending order.
"""
struct PolynomialCover{T<:Real}

struct PolynomialCover{T<:Real}
    coeffs::Vector{Complex{T}}
end

# ------------------------------------------------------------
# Interfaz básica
# ------------------------------------------------------------

coefficients(P::PolynomialCover) = P.coeffs

import Base: degree

degree(P::PolynomialCover) = length(P.coeffs) - 1

# ------------------------------------------------------------
# Constructores
# ------------------------------------------------------------

function PolynomialCover(coeffs::AbstractVector{<:Number})
    T = float(real(eltype(coeffs)))
    PolynomialCover{T}(Complex{T}.(coeffs))
end

PolynomialCover(coeffs::Number...) = PolynomialCover(collect(coeffs))

#=
# ------------------------------------------------------------
# Tipo concreto
# ------------------------------------------------------------

struct PolynomialCover{T<:Real}
    coeffs::Vector{Complex{T}}
end

# ------------------------------------------------------------
# Interfaz básica
# ------------------------------------------------------------

coefficients(P::PolynomialCover) = P.coeffs

# imports necesarios




function degree(P::PolynomialCover)
    length(P.coeffs) - 1
end



function PolynomialCover(coeffs::AbstractVector{<:Number})
    T = promote_type(real(eltype(coeffs)))
    return PolynomialCover{T}(Complex{T}.(coeffs))
end

=#


function value(P::PolynomialCover, z::Number)
    T = promote_type(eltype(P.coeffs), typeof(z))
    acc = zero(T)
    zT = T(z)

    for a in P.coeffs
        acc = acc*zT + T(a)
    end
    return acc
end

function (P::PolynomialCover)(z)
    value(P, z)
end



function derivative(P::PolynomialCover, z::Number)
    n = degree(P)
    T = promote_type(eltype(P.coeffs), typeof(z))
    acc = zero(T)
    zT = T(z)

    for i in 1:n
        acc = acc*zT + T(P.coeffs[i])*(n - i + 1)
    end
    return acc
end


function critical_points(P::PolynomialCover{T}) where T
    n = degree(P)
    n == 0 && return Complex{T}[]
    deriv_coeffs = Complex{T}[(P.coeffs[i]*(n - i + 1)) for i in 1:n]
    deriv_poly = PolynomialCover{T}(deriv_coeffs)
    return roots(deriv_poly)
end



function critical_values(P::PolynomialCover{T}) where T
    cps = critical_points(P)
    return [value(P, cp) for cp in cps]
end

function fiber_via_companion(
    P::PolynomialCover{T},
    y0::Number
) where T

    n = degree(P)
    n == 0 && return Complex{float(T)}[]

    # tipo de cálculo
    Tf = float(T)
    CT = Complex{Tf}

    # valor objetivo en complejo
    y0c = CT(y0)

    # coeficientes de Q(x) = P(x) - y0
    a = CT.(P.coeffs)
    a[end] -= y0c

    # matriz compañera
    C = zeros(CT, n, n)

    # subdiagonal
    @inbounds for i in 2:n
        C[i, i-1] = one(CT)
    end

    # primera fila
    @inbounds for j in 1:n
        C[1, j] = -a[j+1] / a[1]
    end

    return Complex{T}.(eigvals(C))
end
=#