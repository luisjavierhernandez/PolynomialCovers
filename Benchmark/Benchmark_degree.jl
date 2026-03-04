using BenchmarkTools 
using Polynomials
using HomotopyContinuation
using LinearAlgebra
using Printf

BLAS.set_num_threads(1)





###
using PolynomialCovers



###

############################################################
# MATHEMATICALLY CONSISTENT BENCHMARK
############################################################


############################################################
# Polynomial family
############################################################

coeffs(n,T=Float64) = T[1,zeros(T,n-2)...,1,1]
coeffs(7)

poly_companion(n) = Polynomial(coeffs(n))
poly_cover(n)     = PolynomialCover(coeffs(n))

function poly_homotopy(n,y)
    @polyvar x
    x^n + x + 1 - y
end

############################################################
# Critical radius (usar TU función real aquí)
############################################################

function critical_radius(P)
    R= PolynomialCovers.min_critical_radius(PolynomialCovers.coefficients(P))
    S=PolynomialCovers.bound_values_radius_origin(P, R)
    return S
end

############################################################
# Compute global admissible y
############################################################

function admissible_y(ns)

    Svals = Float64[]

    for n in ns
        P = poly_cover(n)
        push!(Svals, critical_radius(P))
    end

    S = maximum(Svals)
    y = 2S

    println("\nComputed global S = ", S)
    println("Chosen |y| = ", abs(y))
    println("3S/2 = ", 1.5S, "\n")

    return y,S
end

############################################################
# Residual
############################################################
(P::PolynomialCover)(x) = value(P,x)

residual(P,r,y) = maximum(abs(P(x)-y) for x in r)

############################################################
# Companion solver
############################################################

function bench_companion(n,y)

    P = poly_companion(n) - y

    roots(P)               # warm-up
    t = @belapsed roots($P)

    r = roots(P)
    err = maximum(abs.(P.(r)))

    return t,err
end

############################################################
# Homotopy solver
############################################################

function bench_homotopy(n,y)

    F = poly_homotopy(n,y)

    solve([F]; show_progress=false)
    t = @belapsed solve([$F]; show_progress=false)

    sol = solve([F]; show_progress=false)
    r = first.(solutions(sol))

    err = maximum(abs.(r.^n .+ r .+ 1 .- y))

    return t,err
end

############################################################
# Cyclic lifting solver
############################################################

function bench_cyclic(n,y)

    P = poly_cover(n)

    cyclic_lift_solver(P,y)
    t = @belapsed PolynomialCovers.cyclic_lift_solver($P,$y)

    r = cyclic_lift_solver(P,y)
    err = residual(P,r,y)

    return t,err
end

############################################################
# Main benchmark
############################################################

function run_benchmarks(ns)

    y,S = admissible_y(ns)

    println(" n | companion | homotopy | cyclic | err_comp | err_hom | err_cyc")
    println("-"^75)

    for n in ns

        tc,ec = bench_companion(n,y)
        th,eh = bench_homotopy(n,y)
        tl,el = bench_cyclic(n,y)

        @printf("%3d | %8.3e | %8.3e | %8.3e | %8.2e | %8.2e | %8.2e\n",
            n,tc,th,tl,ec,eh,el)
    end

    println()
end

############################################################
# Run experiment
############################################################

ns = 5:5:50
y,S = admissible_y(ns)
run_benchmarks(ns)

