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