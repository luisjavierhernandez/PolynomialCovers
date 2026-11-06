using Plots

function plot_path(
    p::Path,
    S::Vector{<:Vertex};
    z0 = first(p.vertices),
    zf = last(p.vertices),
    save::Bool = false,
    filename::AbstractString = "avoid_points_path.png",
)

    zs = getfield.(p.vertices, :z)
    xs = real.(zs)
    ys = imag.(zs)

    plt = plot(
        xs, ys;
        lw = 2,
        marker = :circle,
        label = "path",
        aspect_ratio = :equal,
    )

    scatter!(
        [real(z0.z)], [imag(z0.z)];
        markersize = 8,
        label = "start",
    )

    scatter!(
        [real(zf.z)], [imag(zf.z)];
        markersize = 8,
        label = "target",
    )

    if !isempty(S)
        zsS = getfield.(S, :z)

        scatter!(
            real.(zsS),
            imag.(zsS);
            markershape = :x,
            markersize = 9,
            label = "obstacles",
        )
    end

    if save
        isdir(PLOT_DIR) || mkpath(PLOT_DIR)

        fullpath = joinpath(PLOT_DIR, filename)
        savefig(plt, fullpath)

        @info "Plot saved at $fullpath"
    end

    return plt
end

const PLOT_DIR = joinpath(@__DIR__, "MyPlots", "PlotExamples")

#####

#Basic plot for debugging and visualization of paths, fibers, etc.


####
#P(x) = x^2

using PolynomialCovers




function example_predictor_newton()

    P = PolynomialCover([1,0,0])  # x^2

    x0 = basepoint(P)
    w0 = value(P, x0)

    state = PathState(x0, w0)

    ypath = Path([
        Vertex(w0),
        Vertex(w0 + 0.5 + 0.2im),
        Vertex(w0 + 1.0 + 0.4im)
    ])

    PolynomialCovers.lift_path_trace!(P, state, ypath; S=1.0)

    fig = plot(
        real.(state.pathx),
        imag.(state.pathx),
        marker=:circle,
        aspect_ratio=:equal,
        label="domain trajectory"
    )

    savefig(fig, joinpath(PLOT_DIR,"predictor_newton.png"))

end

example_predictor_newton()
