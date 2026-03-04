using BenchmarkTools
using Polynomials
using LinearAlgebra
using Printf
using Base.Threads

using PolynomialCovers

BLAS.set_num_threads(1)

println("Julia threads = ", Threads.nthreads())

############################################################
# Polynomial family
############################################################

coeffs(n,T=Float64) = T[1,zeros(T,n-2)...,1,1]

poly_companion(n) = Polynomial(coeffs(n))
poly_cover(n)     = PolynomialCover(coeffs(n))

############################################################
# Critical radius
############################################################

function critical_radius(P)

    R = PolynomialCovers.min_critical_radius(
        PolynomialCovers.coefficients(P)
    )

    S = PolynomialCovers.bound_values_radius_origin(P, R)

    return S
end

############################################################
# Compute admissible y
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
    println("|y| = ", abs(y))
    println("3S/2 = ", 1.5S, "\n")

    return y,S
end

############################################################
# Residual
############################################################

(P::PolynomialCover)(x) = value(P,x)

residual(P,r,y) = maximum(abs(P(x)-y) for x in r)

############################################################
# Companion
############################################################

function bench_companion(n,y)

    P = poly_companion(n) - y

    roots(P)  # compile

    t = @belapsed roots($P)

    r = roots(P)

    err = maximum(abs.(P.(r)))

    return t,err
end

############################################################
# Cyclic solver
############################################################

function bench_cyclic(n,y)

    P = poly_cover(n)

    PolynomialCovers.cyclic_lift_solver(P,y)

    t = @belapsed PolynomialCovers.cyclic_lift_solver($P,$y)

    r = PolynomialCovers.cyclic_lift_solver(P,y)

    err = residual(P,r,y)

    return t,err
end

############################################################
# Benchmark loop
############################################################

function run_benchmarks(ns)

    y,S = admissible_y(ns)

    println(" n | companion | cyclic | err_comp | err_cyc")
    println("-"^60)

    tc_vec = Vector{Float64}(undef, length(ns))
    ec_vec = Vector{Float64}(undef, length(ns))
    tl_vec = Vector{Float64}(undef, length(ns))
    el_vec = Vector{Float64}(undef, length(ns))

    Threads.@threads for i in 1:length(ns)
        n = ns[i]
        tc,ec = bench_companion(n,y)
        tl,el = bench_cyclic(n,y)
        tc_vec[i] = tc
        ec_vec[i] = ec
        tl_vec[i] = tl
        el_vec[i] = el
    end

    for i in 1:length(ns)
        n = ns[i]
        @printf("%3d | %8.3e | %8.3e | %8.2e | %8.2e\n",
            n, tc_vec[i], tl_vec[i], ec_vec[i], el_vec[i])
    end

    println()
end

############################################################
# Run
############################################################

function main()

    ns = 5:10:100

    run_benchmarks(ns)

end

main()