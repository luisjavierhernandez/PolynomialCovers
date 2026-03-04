

using PolynomialCovers
include("../MyPlots/MyPlot.jl")
T = Float64

R = 3.0
nS = 12
S = Vertex{T}[
    Vertex{T}(R * cis(2π*k/nS)) for k in 0:nS-1
]
typeof(S) == Vector{Vertex{Float64}}
# ------------------------------------------------
# Obstáculos: 12 puntos en un círculo
# ------------------------------------------------


z0 = Vertex{T}(-5.0 + 0.0im)
zf = Vertex{T}( 5.0 + 0.0im)
typeof(z0) == Vertex{Float64}
  # ← claramente fuera del círculo

p_out = avoid_points_path(
    z0,
    zf,
    S;
    d_safe = 0.8,
    step = 0.05,
)

plot_path(
    p_out,
    S;
    save = true,
    filename = joinpath(PolynomialCovers.PLOT_DIR, "avoid_outside_S.png")
)

@assert all(v.z != zf.z for v in S)

z0 = Vertex(-5.0 + 0.0im)
zf = S[4]     # ← destino exactamente en S

p_in = avoid_points_path(
    z0,
    zf,
    S;
    d_safe = 0.8,
    step = 0.05,
)

plot_path(
    p_in,
    S;
    save = true,
    filename = joinpath(PLOT_DIR, "avoid_ending_in_S.png")
)

@assert any(v.z != zf.z for v in S)