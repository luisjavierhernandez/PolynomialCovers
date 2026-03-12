

using Documenter
using PolynomialCovers

makedocs(
    sitename = "PolynomialCovers.jl",
    modules = [PolynomialCovers],

    format = Documenter.HTML(
        repolink = "https://github.com/luisjavierhernandez/PolynomialCovers.jl"
    ),

    checkdocs = :none,
    warnonly = true,

    pages = [
        "Home" => "index.md",
        "Geometry" => "geometry.md",
        "Polynomial Covers" => "polynomialUtils.md",
        "Critical Disks" => "polynomialCriticalDisks.md",
        "Basepoint and Loops" => "basepoint_and_loops.md",
        "API" => "api.md",
    ],
)