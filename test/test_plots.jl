using Test
using LinearAlgebra
using Plots
using DelimitedFiles
using PolynomialCovers

include("../MyPlots/MyPlot.jl")

@testset "visual debug" begin
    z0 = Vertex(0.0,0.0)
    zf = Vertex(5.0,0.0)

    S = [Vertex(2.5,0.0)]

    p = avoid_points_path(z0,zf,S)

    plot_path(
        p,S;
        save=true,
        filename="debug_obstacle.png"
    )
end