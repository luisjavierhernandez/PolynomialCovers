
#Me recomiendan

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
