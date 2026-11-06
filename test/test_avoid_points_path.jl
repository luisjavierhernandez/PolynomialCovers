############################################################
# test_avoid_points_path.jl
#
# Tests for avoid_points_path
# PolynomialCovers.jl
############################################################

using Test
using LinearAlgebra

# ----------------------------------------------------------
# Load your package
# ----------------------------------------------------------
using PolynomialCovers   # <-- adjust if needed

############################################################
# Auxiliary metrics
############################################################

"""
    min_distance_to_set(path, S)

Minimum distance between path vertices
and forbidden set S.
"""
function min_distance_to_set(p::Path, S)
    minimum(
        abs(v.z - a.z)
        for v in p.vertices
        for a in S
    )
end


"""
    max_turn_angle(path)

Maximum turning angle along polygonal path.
Used to detect violent turns.
"""
function max_turn_angle(p::Path)

    length(p.vertices) ≤ 2 && return 0.0

    θmax = 0.0

    for i in 2:length(p.vertices)-1

        v1 = p.vertices[i].z - p.vertices[i-1].z
        v2 = p.vertices[i+1].z - p.vertices[i].z

        θ = acos(
            clamp(
                real(conj(v1)*v2) /
                (abs(v1)*abs(v2)),
                -1.0,
                1.0
            )
        )

        θmax = max(θmax, θ)
    end

    θmax
end


############################################################
# TEST SETS
############################################################

@testset "avoid_points_path basic" begin

    z0 = Vertex(0.0,0.0)
    zf = Vertex(5.0,0.0)

    S = [
        Vertex(2.5,0.5),
        Vertex(2.5,-0.5)
    ]

    p = avoid_points_path(z0,zf,S)

    @test PolynomialCovers.startpoint(p) == z0.z
    @test PolynomialCovers.endpoint(p)   == zf.z
    @test length(p) > 2
end


@testset "safety distance" begin

    z0 = Vertex(0.0,0.0)
    zf = Vertex(6.0,0.0)

    S = [Vertex(3.0,0.0)]

    p = avoid_points_path(z0,zf,S)

    dmin = min_distance_to_set(p,S)

    @test dmin > 0.2
end


@testset "obstacle bypass" begin

    z0 = Vertex(0.0,0.0)
    zf = Vertex(4.0,0.0)

    S = [Vertex(2.0,0.0)]

    p = avoid_points_path(z0,zf,S)

    ys = [imag(v.z) for v in p.vertices]

    @test maximum(abs.(ys)) > 0.1
end


@testset "smoothness" begin

    p = avoid_points_path(
        Vertex(0,0),
        Vertex(5,0),
        [Vertex(2.5,0.0)]
    )

    @test max_turn_angle(p) < π/2
end


@testset "final allowed collision" begin

    z0 = Vertex(0.0,0.0)
    zf = Vertex(3.0,0.0)

    S = [zf]

    p = avoid_points_path(z0,zf,S)

    hits =
        count(v -> v == zf, p.vertices)

    @test hits == 1
end


@testset "random stability" begin

for k in 1:20

    z0 = Vertex(0.0,0.0)
    zf = Vertex(5.0,0.0)

    S = [
        Vertex(rand()*4, rand()-0.5)
        for _ in 1:5
    ]

    p = avoid_points_path(z0,zf,S)

    @test PolynomialCovers.endpoint(p) == zf.z
end

end


println("\nAll avoid_points_path tests passed ✅")