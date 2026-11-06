############################################################
# test_cyclic_lift_solver.jl
#
# Tests for cyclic_lift_solver — the main algorithm of the paper.
#
# For each test polynomial P and target y with |y| > 3S/2 we verify:
#   1. The solver returns exactly degree(P) roots.
#   2. Every root satisfies |P(x) - y| < tol.
#   3. All roots are distinct.
#   4. The root set matches the companion-matrix roots (up to ordering).
############################################################

using Test
using LinearAlgebra
using PolynomialCovers

# ----------------------------------------------------------
# Helpers
# ----------------------------------------------------------

"""
    companion_roots(P, y)

Roots of P(x) = y computed via the companion matrix (reference solver).
"""
function companion_roots(P, y)
    n = degree(P)
    coeffs = copy(complex.(coefficients(P)))
    coeffs[end] -= complex(y)
    C = zeros(ComplexF64, n, n)
    for i in 2:n
        C[i, i-1] = 1.0
    end
    for j in 1:n
        C[1, j] = -coeffs[j+1] / coeffs[1]
    end
    sort(eigvals(C), by = z -> (real(z), imag(z)))
end

"""
    roots_match(xs, ys; tol=1e-5)

Check whether two root sets agree up to permutation.
"""
function roots_match(xs, ys; tol=1e-5)
    used = falses(length(ys))
    for x in xs
        idx = argmin([used[j] ? Inf : abs(x - ys[j]) for j in 1:length(ys)])
        abs(x - ys[idx]) < tol || return false
        used[idx] = true
    end
    return true
end


# ----------------------------------------------------------
# Common parameters
# ----------------------------------------------------------

RESIDUAL_TOL  = 1e-8   # |P(x) - y| tolerance
DISTINCT_TOL  = 1e-6   # minimum separation between roots
MATCH_TOL     = 1e-5   # tolerance for comparing with companion roots


############################################################
# Degree 2
############################################################

@testset "cyclic_lift_solver — degree 2: z^2 - 1" begin

    P = PolynomialCover([1.0, 0.0, -1.0])
    S = real(critical_values_radius_origin(P))
    y = 4.0 * S   # well outside the critical disk

    xs = cyclic_lift_solver(P, y)

    @test length(xs) == 2
    @test all(abs(value(P, x) - y) < RESIDUAL_TOL for x in xs)
    @test abs(xs[1] - xs[2]) > DISTINCT_TOL
    @test roots_match(xs, companion_roots(P, y); tol=MATCH_TOL)
end

@testset "cyclic_lift_solver — degree 2: z^2" begin

    P = PolynomialCover([1.0, 0.0, 0.0])
    # P(z) = z^2: only critical value is P(0) = 0, so S ≈ 0 and y = 4*S = 0 would
    # fall inside the critical disk. Use a fixed nonzero y instead.
    y = 5.0 + 0im

    xs = cyclic_lift_solver(P, y)

    @test length(xs) == 2
    @test all(abs(value(P, x) - y) < RESIDUAL_TOL for x in xs)
    @test roots_match(xs, companion_roots(P, y); tol=MATCH_TOL)
end


############################################################
# Degree 3
############################################################

@testset "cyclic_lift_solver — degree 3: z^3 - z" begin

    P = PolynomialCover([1.0, 0.0, -1.0, 0.0])
    S = real(critical_values_radius_origin(P))
    y = 4.0 * S

    xs = cyclic_lift_solver(P, y)

    @test length(xs) == 3
    @test all(abs(value(P, x) - y) < RESIDUAL_TOL for x in xs)
    for i in 1:3, j in i+1:3
        @test abs(xs[i] - xs[j]) > DISTINCT_TOL
    end
    @test roots_match(xs, companion_roots(P, y); tol=MATCH_TOL)
end

@testset "cyclic_lift_solver — degree 3: monic sparse" begin

    # P(z) = z^3 + z + 2
    P = PolynomialCover([1.0, 0.0, 1.0, 2.0])
    S = real(critical_values_radius_origin(P))
    y = 4.0 * S

    xs = cyclic_lift_solver(P, y)

    @test length(xs) == 3
    @test all(abs(value(P, x) - y) < RESIDUAL_TOL for x in xs)
    @test roots_match(xs, companion_roots(P, y); tol=MATCH_TOL)
end


############################################################
# Degree 5  (benchmark polynomial from the paper)
############################################################

@testset "cyclic_lift_solver — degree 5: z^5 - z + 1" begin

    P = PolynomialCover([1.0, 0.0, 0.0, 0.0, -1.0, 1.0])
    S = real(critical_values_radius_origin(P))
    y = 4.0 * S

    xs = cyclic_lift_solver(P, y)

    # Degree-5 path lifting accumulates more numerical error than degree 2–3,
    # so we use slightly relaxed tolerances here.
    DEG5_RESIDUAL_TOL = 1e-6
    DEG5_MATCH_TOL    = 1e-4

    @test length(xs) == 5
    @test all(abs(value(P, x) - y) < DEG5_RESIDUAL_TOL for x in xs)
    for i in 1:5, j in i+1:5
        @test abs(xs[i] - xs[j]) > DISTINCT_TOL
    end
    @test roots_match(xs, companion_roots(P, y); tol=DEG5_MATCH_TOL)
end


############################################################
# Large independent term (regime of the paper)
############################################################

@testset "cyclic_lift_solver — large |y| regime" begin

    P = PolynomialCover([1.0, 0.0, -1.0])
    S = real(critical_values_radius_origin(P))

    # For large |y|, absolute residual scales with |y|, so we use a relative
    # tolerance: |P(x) - y| / |y| < REL_TOL.
    REL_TOL = 1e-7

    for y in [10.0*S, 50.0*S, 100.0*S]
        xs = cyclic_lift_solver(P, y)
        @test length(xs) == 2
        @test all(abs(value(P, x) - y) < REL_TOL * abs(y) for x in xs)
    end
end

@testset "cyclic_lift_solver — complex target" begin

    P = PolynomialCover([1.0, 0.0, -1.0, 0.0])
    S = real(critical_values_radius_origin(P))
    y = 4.0 * S * cis(π/4)   # outside disk, complex direction

    xs = cyclic_lift_solver(P, y)

    @test length(xs) == 3
    @test all(abs(value(P, x) - y) < RESIDUAL_TOL for x in xs)
end


############################################################
# Error handling
############################################################

@testset "cyclic_lift_solver — rejects unsafe target" begin

    P = PolynomialCover([1.0, 0.0, -1.0])
    S = real(critical_values_radius_origin(P))

    # y inside the critical disk: must throw
    @test_throws Exception cyclic_lift_solver(P, 0.5 * S)
end

println("\nAll cyclic_lift_solver tests passed ✅")
