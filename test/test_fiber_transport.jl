############################################################
# test_fiber_transport.jl
#
# Tests for the fiber/transport workflow:
#   PolynomialCover → basepoint → fiber → transport! / transport_to!
############################################################

using Test
using LinearAlgebra
using PolynomialCovers

# ----------------------------------------------------------
# Helpers
# ----------------------------------------------------------

"""
    max_residual(P, xs, y)

Maximum residual |P(x) - y| over a vector of points.
"""
function max_residual(P, xs, y)
    maximum(abs(value(P, x) - y) for x in xs)
end


############################################################
# PolynomialCover basics
############################################################

@testset "PolynomialCover basics" begin

    # P(z) = z^2 - 1  →  coefficients [1, 0, -1]
    P = PolynomialCover([1.0, 0.0, -1.0])

    @test degree(P) == 2
    @test length(coefficients(P)) == 3
    @test value(P, 1.0) ≈ 0.0
    @test value(P, -1.0) ≈ 0.0
    @test value(P, 0.0) ≈ -1.0

    # Degree-3: P(z) = z^3 - z
    Q = PolynomialCover([1.0, 0.0, -1.0, 0.0])
    @test degree(Q) == 3
end


############################################################
# Critical value geometry
############################################################

@testset "Critical value geometry" begin

    # P(z) = z^2 - 1
    P = PolynomialCover([1.0, 0.0, -1.0])

    # Disk containing all critical values
    D = critical_values_disk_origin(P)
    @test D isa Disk
    @test real(D.center) ≈ 0.0
    @test D.radius > 0.0

    # Radius is consistent
    r = critical_values_radius_origin(P)
    @test r ≈ D.radius

    # For P(z) = z^2 - 1 the only critical value is P(0) = -1
    @test r ≥ 1.0
end


############################################################
# Basepoint
############################################################

@testset "Basepoint" begin

    P = PolynomialCover([1.0, 0.0, -1.0])

    b = basepoint(P)

    @test b isa Complex
    @test imag(b) ≈ 0.0   # on the positive real axis

    D = critical_values_disk_origin(P)
    @test abs(b) > D.radius   # lies outside the critical value disk
end


############################################################
# Basepoint initial loop
############################################################

@testset "Basepoint initial loop" begin

    P = PolynomialCover([1.0, 0.0, -1.0])

    loop = basepoint_initial_loop(P)

    @test loop isa Path
    @test length(loop.vertices) ≥ 4           # at least a square

    # loop starts and ends at the same point (closed)
    start = loop.vertices[1].z
    finish = loop.vertices[end].z
    @test abs(start - finish) < 1e-12
end


############################################################
# Fiber computation
############################################################

@testset "Fiber — degree 2" begin

    P = PolynomialCover([1.0, 0.0, -1.0])   # z^2 - 1

    F = fiber(P)

    @test F isa Fiber
    @test length(F.xs) == degree(P)           # exactly n sheets

    # All points satisfy P(x) = w
    @test max_residual(P, F.xs, F.w) < 1e-8
end

@testset "Fiber — degree 3" begin

    P = PolynomialCover([1.0, 0.0, -1.0, 0.0])   # z^3 - z

    F = fiber(P)

    @test length(F.xs) == 3
    @test max_residual(P, F.xs, F.w) < 1e-8

    # All sheets are distinct
    for i in 1:3, j in i+1:3
        @test abs(F.xs[i] - F.xs[j]) > 1e-6
    end
end

@testset "Fiber — degree 5" begin

    P = PolynomialCover([1.0, 0.0, 0.0, 0.0, -1.0, 1.0])   # z^5 - z + 1

    F = fiber(P)

    @test length(F.xs) == 5
    @test max_residual(P, F.xs, F.w) < 1e-8
end


############################################################
# safe_path_arc
############################################################

@testset "safe_path_arc" begin

    P = PolynomialCover([1.0, 0.0, -1.0])
    S = real(critical_values_radius_origin(P))

    u = complex(2S, 0.0)
    v = complex(-2S, 0.0)

    γ = safe_path_arc(u, v, S)

    @test γ isa Path
    @test γ.vertices[1].z ≈ u
    @test γ.vertices[end].z ≈ v

    # No vertex enters the forbidden disk
    for vert in γ.vertices
        @test abs(vert.z) > 0.9 * S
    end
end


############################################################
# transport!
############################################################

@testset "transport! preserves fiber size and residual" begin

    P = PolynomialCover([1.0, 0.0, -1.0])

    F = fiber(P)
    n = length(F.xs)
    S = F.S

    # Path between two safe points
    u = F.w
    v = complex(2 * abs(u), 0.0)
    γ = safe_path_arc(u, v, S)

    transport!(F, γ; parallel=false)

    @test length(F.xs) == n
    @test max_residual(P, F.xs, F.w) < 1e-8
    @test F.w ≈ v
end

@testset "transport! monodromy permutes degree-2 fiber" begin

    # For P(z) = z^2 the monodromy loop exchanges the two sheets
    P = PolynomialCover([1.0, 0.0, 0.0])

    F = fiber(P)
    xs_before = copy(F.xs)

    loop = basepoint_initial_loop(P)
    transport!(F, loop; parallel=false)

    # The set of sheets must be the same (possibly permuted)
    for x in xs_before
        @test any(abs(x - y) < 1e-6 for y in F.xs)
    end
end


############################################################
# transport_to!
############################################################

@testset "transport_to! reaches target" begin

    P = PolynomialCover([1.0, 0.0, -1.0])

    F = fiber(P)
    S = F.S
    y_target = 4 * S   # safely outside the critical value disk

    transport_to!(F, y_target)

    @test F.w ≈ y_target
    @test max_residual(P, F.xs, F.w) < 1e-8
end

@testset "transport_to! rejects unsafe target" begin

    P = PolynomialCover([1.0, 0.0, -1.0])
    F = fiber(P)
    S = F.S

    # Target inside the critical value disk: must throw
    @test_throws ArgumentError transport_to!(F, 0.5 * S)
end
