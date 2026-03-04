using BenchmarkTools
using Polynomials
using HomotopyContinuation
using LinearAlgebra
using Printf
using Base.Threads

using PolynomialCovers

Threads.@threads for i in 1:10
    println("thread ", Threads.threadid())
end

############################################################
# Polynomial P25
############################################################

function coeffs25(T=Float64)
    c = zeros(T,26)
    c[1]  = 1
    c[2]  = 1
    c[26] = 1
    return c
end

poly_companion() = Polynomial(coeffs25())
poly_cover()     = PolynomialCover(coeffs25())

function poly_homotopy(y)
    @polyvar x
    x^25 + x + 1 - y
end

############################################################
# Critical radius
############################################################

function critical_radius(P)

    R = PolynomialCovers.min_critical_radius(
        PolynomialCovers.coefficients(P)
    )

    S = PolynomialCovers.bound_values_radius_origin(P,R)

    return S
end

############################################################
# Sequence |y| > 3S/2
############################################################

function y_sequence(P,mmax; α=1.4)

    S = critical_radius(P)

    ys = ComplexF64[]

    y0 = 1.6*S          # > 3S/2

    for m in 0:mmax
        push!(ys, y0 * α^m)
    end

    println()
    println("Critical value bound S = ",S)
    println("Minimum |y| used       = ",abs(first(ys)))
    println("3S/2                   = ",1.5S)
    println()

    return ys,S
end

############################################################
# Residual
############################################################

(P::PolynomialCover)(x) = value(P,x)

residual(P,r,y) = maximum(abs(P(x)-y) for x in r)

############################################################
# Companion
############################################################

function bench_companion(P,y)

    Q = Polynomial(PolynomialCovers.coefficients(P)) - y

    roots(Q)
    t = @belapsed roots($Q)

    r = roots(Q)
    err = maximum(abs.(Q.(r)))

    return t,err
end

############################################################
# HomotopyContinuation
############################################################

function bench_homotopy(y)

    F = poly_homotopy(y)

    solve([F]; show_progress=false)
    t = @belapsed solve([$F]; show_progress=false)

    sol = solve([F]; show_progress=false)
    r = first.(solutions(sol))

    err = maximum(abs.(r.^25 .+ r .+ 1 .- y))

    return t,err
end

############################################################
# Cyclic lifting
############################################################

function bench_cyclic(P,y)

    cyclic_lift_solver(P,y)
    t = @belapsed PolynomialCovers.cyclic_lift_solver($P,$y)

    r = cyclic_lift_solver(P,y)
    err = residual(P,r,y)

    return t,err
end

############################################################
# Main benchmark
############################################################

function run_benchmark(mmax)

    P = poly_cover()

    ys,S = y_sequence(P,mmax)

    println(" m | |y| | companion | homotopy | cyclic | err_comp | err_hom | err_cyc")
    println("-"^95)

    for (m,y) in enumerate(ys)

        tc,ec = bench_companion(P,y)
        th,eh = bench_homotopy(y)
        tl,el = bench_cyclic(P,y)

        @printf("%3d | %8.2e | %8.3e | %8.3e | %8.3e | %8.2e | %8.2e | %8.2e\n",
            m-1,abs(y),tc,th,tl,ec,eh,el)
    end

    println()
end

############################################################
# Run
############################################################

run_benchmark(12)

function stress_test()

    n = 25
    P = poly_cover(n)

    S = critical_radius(P)

    ys = [2S * 2.0^m for m in 0:20]

    println(" |y|  companion_err  cyclic_err")

    for y in ys

        _,ec = bench_companion(n,y)
        _,el = bench_cyclic(n,y)

        @printf("%e  %e  %e\n", abs(y), ec, el)

    end
end