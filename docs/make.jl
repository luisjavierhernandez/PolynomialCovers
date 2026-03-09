using Documenter
using PolynomialCovers

makedocs(
    sitename = "PolynomialCovers.jl",
    modules = [PolynomialCovers],
    format = Documenter.HTML(),
    checkdocs = :none,
    pages = [
        "Home" => "index.md",
        # "Theory" => "theory.md",
        # "Geometry" => "geometry.md",
        # "Polynomial Covers" => "polynomial_covers.md",
        # "Critical Geometry" => "critical_geometry.md",
        # "Lifting" => "lifting.md",
        # "Cyclic Solver" => "cyclic_solver.md",
        "API" => "api.md",
    ],
)